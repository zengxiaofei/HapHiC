## HapHiC: a fast, reference-independent, allele-aware scaffolding tool based on Hi-C data



![](./images/HapHiC1.png)

<center><strong>[ English | <a href="./README_cn.md)">简体中文</a> ]</strong></center>

HapHiC is an allele-aware scaffolding tool that uses Hi-C data to scaffold haplotype-phased genome assemblies into chromosome-scale pseudomolecules. Unlike [ALLHiC](https://github.com/tangerzhang/ALLHiC), another allele-aware scaffolder, HapHiC can achieve this without the need for reference genomes. Our evaluations indicate that HapHiC outperforms other Hi-C scaffolding tools with higher tolerance to low contig N50, low Hi-C sequencing depth, and various types of assembly errors. Additionally, HapHiC is super-fast and also suitable for haplotype-collapsed diploid and allopolyploid genome assemblies.

**Features:**

- [x] Chromosome-level scaffolding of haplotype-phased assemblies without reference genomes
- [x] Efficient correction of chimeric contigs (misjoins) with little impact on contig N50
- [x] Much higher tolerance to chimeric contigs, collapsed contigs, and switch errors
- [x] Improved performance in chromosome assignment of contigs
- [x] Improved performance in ordering and orienation of contigs
- [x] Super-fast and memory-efficient
- [x] Able to order and orient contigs without prior knowledge of the number of chromosomes (quick view mode)
- [x] Able to utilize phasing information from hifiasm with varying confidence levels
- [x] Extensive compatibility and user-friendly interface: supports chromap; provides a built-in one-command pipeline; able to produce highly customizable vector graphics for contact maps

**Recent updates:**

* Version 1.0.6 (2024.08.26): There is no longer a need to manually set the scale factor in Juicebox. In addition, the saved `.review.assembly` file can now be parsed correctly by Juicebox.
* Version 1.0.5 (2024.07.05): Improved stability in ordering and orientation of contigs through a comparison of fast sorting and ALLHiC optimization.
* Version 1.0.4 (2024.07.03): Add a `haphic refsort` command for [ordering and orienting whole scaffolds according to a reference genome](#refsort).
* Version 1.0.3 (2024.03.21): Add support for the [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) used in [chromap](https://github.com/haowenz/chromap).
* Version 1.0.2 (2023.12.08): We have introduced a `haphic plot` command for [Hi-C contact map visualization](#visualization).
* Version 1.0.1 (2023.11.30): Improved AGP output by incorporating a YaHS-style `scaffolds.raw.agp` for compatibility with the [Juicebox visualization](#juicebox) method suggested by YaHS.

**Terminology:** To ensure conciseness and clarity, we use the term "contigs" to refer to the fragmented genome sequences in the input assembly, although they could be either contigs or scaffolds in actuality.



## Table of contents

- [Installation](#installation)
- [Quick start](#quick_start)
  * [Align Hi-C data to the assembly](#align)
  * [Run HapHiC scaffolding pipeline](#pipeline)
- [Go through the pipeline step by step](#step_by_step)
  * [[Step 1] Clustering](#step1)
  * [[Step 2] Reassignment](#step2)
  * [[Step 3] Ordering and orientation](#step3)
  * [[Step 4] Building scaffolds](#step4)
- [Examples](#examples)
- [Work with hifiasm](#hifiasm)
- [Quick view mode](#quick_view)
- [Juicebox curation](#juicebox)
- [Visualization](#visualization)
- [Order and orient whole scaffolds using a reference genome](#refsort)
- [Frequently asked questions (FAQs)](#faqs)
- [Problems and bug reports](#problems)
- [Citing HapHiC](#citing)
- [Reproducibility](#reproduce)



## <span id="installation">Installation</span>

HapHiC has been tested and validated on servers running Linux, equipped with either Intel Xeon, AMD EPYC, or Hygon C86 CPUs.

```bash
# (1) Download HapHiC from GitHub
$ git clone https://github.com/zengxiaofei/HapHiC.git

# (2) Resolve dependencies
# We strongly recommend using conda to install dependencies. If you prefer manual installation, refer to HapHiC/conda_env/create_conda_env_py310.sh
# We have also included additional environments for Python 3.11 and 3.12 in the directory HapHiC/conda_env/
$ conda env create -f HapHiC/conda_env/environment_py310.yml
# Activate the HapHiC conda environment
$ conda activate haphic # or: source /path/to/conda/bin/activate haphic

# (3) Check whether all dependencies are correctly installed
$ /path/to/HapHiC/haphic check

# (4) Show all available commands and help message
$ /path/to/HapHiC/haphic -h
```



## <span id="quick_start">Quick start</span>

### <span id="align">Align Hi-C data to the assembly</span>

First, you need to prepare a BAM file by aligning Hi-C data to the assembly. Here is the way that we recommend:

```bash
# (1) Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
$ bwa index asm.fa
$ bwa mem -5SP -t 28 asm.fa /path/to/read1_fq.gz /path/to/read2_fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o HiC.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
$ /path/to/HapHiC/utils/filter_bam HiC.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o HiC.filtered.bam
```

> [!NOTE]
>
> * Here, `asm.fa` can be haplotype-collapsed contigs (e.g., `p_ctg` in hifiasm), haplotype-phased unitigs (e.g., `p_utg` in hifiasm), or one or more sets of haplotype-resolved contigs (e.g., `hap*.p_ctg` in hifiasm). In addition, `asm.fa` may also be scaffolds output by other scaffolders.
> * You can prepare the BAM file according to your own preferences or requirements, but **DO NOT** sort it by coordinate. If your BAM file is already sorted by coordinate, you need to resort it by read name (`samtools sort -n`).
> * We **DO NOT** recommend the Juicer pipeline for Hi-C reads alignment, particularly in haplotype-phased assemblies.

### <span id="pipeline">Run HapHiC scaffolding pipeline</span>

**(i) One-line command.** HapHiC provides a one-line command `haphic pipeline` to execute the entire scaffolding pipeline. The required parameters are：

1) `asm.fa`, your genome assembly file in FASTA format.
2) `HiC.filtered.bam`, the BAM file prepared in the previous step (the .pairs file output by chromap is also acceptable since version 1.0.3).
3) `nchrs`, the number of chromosomes present in the assembly, and also the expected number of output scaffolds.

```bash
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs
```

**(ii) Restriction site.** The default restriction site is `GATC` (MboI/DpnII). You can modify this using the `--RE` parameter. If you are unsure or if your Hi-C library was constructed without restriction enzymes (REs), it is acceptable to leave it as the default.

```bash
# For HindIII
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --RE "AAGCTT"
# For Arima two-enzyme chemistry
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --RE "GATC,GANTC"
# For Arima four-enzyme chemistry
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --RE "GATC,GANTC,CTNAG,TTAA"
```

**(iii) Contig correction.** To correct misjoined contigs based on Hi-C linking information, use `--correct_nrounds` to enable assembly correction and set the number of correction rounds. For example:

```bash
# Typically, two rounds of assembly correction are enough
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --correct_nrounds 2
```

**(iv) Switch error.** If your input assembly is haplotype-phased and has a high switch error rate (often introduced by assemblers when the sequence divergence between haplotypes is very low), use `--remove_allelic_links` to remove Hi-C links between allelic contigs, thereby increasing tolerance to such errors. The value should be the ploidy of the assembly. For example:

```bash
# For haplotype-phased assembles of autotetraploids, set the parameter to 4
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --remove_allelic_links 4
```

> [!NOTE]
> If your input assembly is haplotype-phased and the Hi-C reads are aligned using other methods like chromap, we also recommend including this parameter to mitigate the adverse effects of incorrect mapping.

**(v) Performance.** Use `--threads` to set the number of threads for BAM file reading, and `--processes` to create multiple processes for contig ordering and orientation. For example:

```bash
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --threads 8 --processes 8
```

**Parameters**

For more information, run `haphic pipeline --help`.

**Final outputs**

* `01.cluster/corrected_asm.fa`: The corrected contigs in FASTA format. This file is generated only when assembly correction is enabled.
* `04.build/scaffolds.agp`: A [SALSA](https://github.com/marbl/SALSA)-style [AGP file](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) containing information about scaffold assignment, ordering and orientation information for all sequences in `corrected_asm.fa`. If there are any chimeric contigs that are corrected by HapHiC, the broken contigs will be assigned new IDs.
* `04.build/scaffolds.raw.agp`: A [YaHS](https://github.com/c-zhou/yahs)-style AGP file containing information about scaffold assignment, ordering and orientation information for all sequences in `asm.fa`. The broken contigs will not be assigned new IDs, but their starting and ending coordinates in the raw contigs will be displayed in the seventh and eighth columns.
* `04.build/scaffolds.fa`: The final scaffolds in FASTA format.
* `04.build/juicebox.sh`: A shell script for [Juicebox visualization and curation](#juicebox).

> [!NOTE]
>
> Although the one-line command is convenient, the automatic parameter tuning may fail, leading to poor results or even a pipeline interruption in rare cases. If this occurs, we recommend [running each step individually with manual parameter tuning](#step_by_step) or trying the [quick view mode](#quick_view) described below.



## <span id="step_by_step">Go through the pipeline step by step</span>

### <span id="step1">[Step 1] Clustering</span>

Before clustering, HapHiC performs preprocessing to correct assembly misjoins, filter out short, mis-assembled contigs, and remove allelic Hi-C links. After that, the Markov cluster algorithm ([MCL algorithm](https://micans.org/mcl/)) is used to cluster the contigs into groups. Unlike the agglomerative hierarchical clustering (AHC, used in LACHESIS and ALLHiC), which specifies the number of clusters, the MCL Algorithm implicitly controls it with a parameter called "inflation". The higher the inflation, the more the groups are clustered. The main problem with AHC is that even though the number of clusters is specified, contigs from different chromosomes may also be clustered into the same group. This is common in phased diploid or polyploid genome assemblies. To solve this, HapHiC tries a series of inflations to cluster the contigs (controlled by `--min_inflation`, `--max_inflation`, and `--inflation_step`) and recommends a "best" one based on both the expected number of chromosomes `nchrs` provided and the length distribution of the groups.

```bash
$ /path/to/HapHiC/haphic cluster asm.fa HiC.filtered.bam nchrs
```

**Parameters**

For more information, run `haphic cluster --help`.

**Main outputs**

* `corrected_asm.fa`: The corrected assembly in FASTA format. This file is generated only when assembly correction is enabled.

* `corrected_ctgs.txt`: A text file listing the IDs of all corrected contigs.

* `full_links.pkl`: A binary file that stores the number of Hi-C links between each contig pair.

* `HT_links.pkl`: A binary file that records the number of Hi-C links between each half of contig pairs.

* `paired_links.clm`: A text file recording the positional information of paired Hi-C links.

* `inflation_*`: Output directories for respective inflations.

  ├── `group*.txt`: Files containing the contigs and their basic information (including lengths and numbers of restriction sites) for each group, also referred to as `counts_RE.txt` in ALLHiC.

  └── `mcl_inflation_*.clusters.txt`: Markov clustering results.

**The "best" inflation**

You can find the "best" inflation recommendation in the log file `HapHiC_cluster.log`, like:

```
2022-11-07 17:50:08 <HapHiC_cluster.py> [recommend_inflation] You could try inflation from 1.20 (length ratio = 0.75)
```

In some cases, HapHiC cannot get the "best" one. It could be due to inappropriate parameters or extensive assembly errors. Check whether the parameters used are correct/appropriate and then try to tune the parameters for assembly correction, contig/Hi-C link filtration, or Markov Clustering:

```
2022-11-19 13:20:38 <HapHiC_cluster.py> [recommend_inflation] It seems that some chromosomes were grouped together (length ratio = 0.5)...
```

### <span id="step2">[Step 2] Reassignment</span>

In the previous step, some contigs may have been filtered out before clustering or assigned to incorrect groups. Additionally, the number of final clusters output by Markov clustering may exceed the specified number of chromosomes (`nchrs`). To address these issues, we added a reassignment step to rescue contigs that are not in any groups, reassign contigs to the correct groups, and perform an additional agglomerative hierarchical clustering to concatenate groups if necessary. The input files `full_links.pkl`, `mcl_inflation_x.clusters.txt`, and `paired_links.clm` are outputs from the clustering step, where `x` represents the "best" inflation value:

```bash
$ /path/to/HapHiC/haphic reassign asm.fa full_links.pkl mcl_inflation_x.clusters.txt paired_links.clm --nclusters nchrs
```

> [!NOTE]
>
> If assembly correction has been performed, use `corrected_asm.fa` as input FASTA file instead of `asm.fa`.

**Parameters**

For more information, run `haphic reassign --help`.

**Main outputs**

* `final_groups/group*.txt`: Files containing the contigs and their basic information for each final group after reassignment.
* `final_groups/final_cluster.txt`: The final clustering result.
* `split_clms/`: A directory containing group-specific CLM files.

### <span id="step3">[Step 3] Ordering and orientation</span>

The ordering and orientation step in HapHiC is implemented using an integration of algorithms from [3D-DNA](https://github.com/aidenlab/3d-dna) and [ALLHiC](https://github.com/tanghaibao/allhic). First, an efficiency-improved 3D-DNA iterative scaffolding algorithm (refered to as "fast sorting") is used to quickly order and orient the contigs. Then, the order and orientation of contigs are input as an initial configuration and optimized with the ALLHiC program ([a modified version](http://github.com/zengxiaofei/allhic), in which the hot-start optimization has been fixed). The input file `HT_links.pkl` is the output file from the clustering step; the directory `split_clms` and the group files `final_groups/group*.txt` were created in the reassignment step. The optional parameter `--processes` is used to set the number of processes for the ordering and orientation.

```bash
$ /path/to/HapHiC/haphic sort asm.fa HT_links.pkl split_clms final_groups/group*.txt --processes 8
```

> [!NOTE]
>
> If assembly correction has been performed, use `corrected_asm.fa` as input FASTA file instead of `asm.fa`.

**Parameters**

For more information, run `haphic sort --help`.

**Main outputs**

* `group*.tour.sav`: The fast sorting result of contigs within each group.
* `group*.tour`: The contig ordering and orientation result for each group after ALLHiC optimization.
* `final_tours/group*.tour`: The final contig ordering and orientation results.

### <span id="step4">[Step 4] Building scaffolds</span>

The final step is to build the scaffolds (pseudomolecules) using the chromosome assignment, ordering and orientation information of contigs from the `group*.tour` files. By default, the output scaffolds are sorted by scaffold length.

**If assembly correction was not performed:**

```bash
$ /path/to/HapHiC/haphic build asm.fa asm.fa HiC.filtered.bam final_tours/group*.tour
```

**If assembly correction has been performed:**

Use `corrected_asm.fa` as input FASTA file instead of **the first** `asm.fa`. Additionally, specify the corrected contig list `corrected_ctgs.txt` using the `--corrected_ctgs` parameter. Otherwise, the YaHS-style `scaffolds.raw.agp` generated may be incorrect.

```bash
$ /path/to/HapHiC/haphic build corrected_asm.fa asm.fa HiC.filtered.bam final_tours/group*.tour --corrected_ctgs corrected_ctgs.txt
```

> [!NOTE]
>
> The second `asm.fa` (raw uncorrected assembly) and `HiC.filtered.bam` are required since HapHiC version 1.0.1 for generating the shell script for juicebox visualization and curation.

**Parameters**

For more information, run `haphic build --help`.

**Main outputs**

* `scaffolds.agp`: A SALSA-style AGP file containing information about scaffold assignment, ordering and orientation information for all sequences in `corrected_asm.fa`. If there are any chimeric contigs that are corrected by HapHiC, the broken contigs will be assigned new IDs.
* `scaffolds.raw.agp`: A YaHS-style AGP file containing information about scaffold assignment, ordering and orientation information for all sequences in `asm.fa`. The broken contigs will not be assigned new IDs, but their starting and ending coordinates in the raw contigs will be displayed in the seventh and eighth columns.
* `scaffolds.fa`: The final scaffolds in FASTA format.
* `juicebox.sh`: A shell script for [Juicebox visualization and curation](#juicebox).



## <span id="examples">Examples</span>

HapHiC can scaffold most genomes **within 1 hour** using only 8 CPU cores. For large genomes with fragmented contigs, scaffolding typically takes less than half a day. HapHiC has been successfully validated in scaffolding genomes from various taxa, including higher plants, humans, birds, amphibians, fish, insects, mollusks, and annelids. For more examples and detailed information, please refer to the [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41477-024-01755-3/MediaObjects/41477_2024_1755_MOESM1_ESM.pdf) in our [paper](https://doi.org/10.1038/s41477-024-01755-3).

| Species              | Karyotype    | Haplotype-resolved | Assembly size (Gb) | Contig N50 (Mb) | Number of contigs | Hi-C depth after filtering | Wall time (min) | Peak RAM (GiB) |
| :-------------------- | ------------ | ------------------ | ------------------ | --------------- | ----------------- | -------------------------- | --------------- | -------------- |
| Giant Miscanthus     | 2*n*=3*x*=57 | Yes                | 6.11               | 2.19            | 5,761             | 33.58×                     | 115.35          | 17.10          |
| Potato C88           | 2*n*=4*x*=48 | Yes                | 3.16               | 18.78           | 2,490             | 13.4×                      | 20.15           | 5.98           |
| Wild sugarcane Np-X  | 2*n*=4*x*=40 | Yes                | 2.76               | 0.38            | 15,510            | 23.7×                      | 78.97           | 27.02          |
| Alfalfa XinJiangDaYe | 2*n*=4*x*=32 | Yes                | 3.16               | 0.46            | 31,772            | 10.1×                      | 33.13           | 7.68           |
| Tea plant Tieguanyin | 2*n*=2*x*=30 | Yes                | 5.99               | 0.22            | 60,345            | 9.8×                       | 157.53          | 33.68          |
| Human HG002          | 2*n*=2*x*=46 | Yes                | 6.02               | 73.40           | 1,153             | 4.7×                       | 13.42           | 11.50          |
| Wheat         | 2*n*=6*x*=42 | No                 | 14.0               | 2.16            | 12,982            | 1.5×                       | 58.05           | 22.98          |
| Ginkgo tree          | 2*n*=2*x*=24 | No                 | 9.87               | 1.58            | 261,820           | 54.1×                      | 440.78          | 135.83         |
| Northern goshawk     | 2*n*=2*x*=80 | No                 | 1.40               | 17.71           | 638               | 27.2×                      | 16.95           | 2.19           |
| Tropical clawed frog | 2*n*=2*x*=20 | No                 | 1.48               | 0.38            | 9,631             | 47.5×                      | 53.80           | 19.83          |
| Corkwing warsse    | 2*n*=2*x*=46 | No                 | 0.64               | 1.19            | 1,774             | 85.3×                      | 25.73           | 3.13           |
| Chinese oak silkmoth | 2*n*=2*x*=98 | No                 | 0.73               | 0.17            | 9,824             | 70.8×                      | 35.33           | 10.66          |
| Gray topshell        | 2*n*=2*x*=36 | No                 | 1.27               | 6.20            | 843               | 58.3×                      | 27.07           | 5.06           |
| Humus earthworm      | 2*n*=2*x*=36 | No                 | 0.79               | 0.71            | 2,261             | 64.4×                      | 20.32           | 3.23           |



## <span id="hifiasm">Work with hifiasm</span>

When scaffolding a phased [hifiasm](https://github.com/chhylp123/hifiasm) assembly, you can run HapHiC with the GFA file(s) output by hifiasm. Here, the term "phased hifiasm assembly" refers to the haplotype-resolved primary contigs assembled via the trio binning or Hi-C-based algorithm (`*.hap*.p_ctg.gfa`), as well as the phased unitigs (`*.p_utg.gfa`). 

HapHiC uses the read depth information in the GFA file(s) to filter out potential collapsed contigs/unitigs before clustering.  If more than one GFA file is provided, HapHiC assumes these GFA files are haplotype-specific (`*.hap*.p_ctg.gfa`), and artificially removes or reduces the Hi-C links between the haplotypes according to this phasing information during clustering. Note that the contigs/unitigs in GFA file(s) should match those in FASTA file. Either `.gfa` or `noseq.gfa` is acceptable.

```bash
# (1) For hifiasm primary unitigs (`*.p_utg.gfa`), use the read depth informtion in the GFA file to filter out potential collapsed unitigs before clustering
$ /path/to/HapHiC/haphic pipeline p_utg.fa HiC.filtered.bam nchrs --gfa p_utg.gfa

# (2) In addition to read depth filtering, HapHiC can also reduce Hi-C links between haplotypes according to phasing information in GFA files for haplotype-resolved primary contigs (`*.hap*.p_ctg.gfa`)

# By default, all Hi-C links between haplotypes are completely removed and contigs from different haplotypes will not be clustered into the same group
$ /path/to/HapHiC/haphic pipeline allhaps.fa HiC.filtered.bam nchrs --gfa "hap1.p_ctg.gfa,hap2.p_ctg.gfa"

# Use the `--phasing_weight` parameter to control the confidence level of the hifiasm phasing information. The weight can be set to 0 to completely ignore the phasing information. You can also set it between 0 and 1 to run HapHiC (considering both hifiasm and HapHiC phasing results), in which case contigs from different haplotypes might be clustered together
$ /path/to/HapHiC/haphic pipeline allhaps.fa HiC.filtered.bam nchrs --gfa "hap1.p_ctg.gfa,hap2.p_ctg.gfa" --phasing_weight 0
```

> [!NOTE]
>
> This feature is not essential; you can decide whether to use the `--gfa` parameter based on your specific needs. First, even without a GFA file, HapHiC can still filter out collapsed contigs/unitigs based on the Hi-C sequencing depth. Second, artificially removing Hi-C links between haplotypes by using phasing information from the GFA files might lead to a certain degree of reduction in the anchoring rate.



## <span id="quick_view">Quick view mode</span>

You can try the quick view mode in HapHiC when: 
1. The exact number of chromosomes is unknown.
2. HapHiC cannot provide an acceptable clustering result or encounters a pipeline interruption.
3. You need a quick view of your assembly (e.g., to identify the type and approximate proportion of assembly errors).
4. You just want to manually curate your assembly and split chromosomes in Juicebox by yourself. 

In quick view mode, HapHiC simply uses the fast sorting to order and orient all contigs without clustering. The result is similar to `*.0.hic` in [3D-DNA](https://github.com/aidenlab/3d-dna). Most parameters are disabled in this mode, but you can use `--correct_nrounds` to correct input contigs. When scaffolding a haplotype-resolved hifiasm assembly ( `*.hap*.p_ctg.gfa` ), you can still partition contigs into different haplotypes with the haplotype-specific GFA files.

```bash
# HapHiC will ignore the parameter `nchrs`, it can be any integer
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --quick_view
# Correct input contigs in quick view mode
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --quick_view --correct_nrounds 2
# Partition contigs into different haplotypes in quick view mode
$ /path/to/HapHiC/haphic pipeline allhaps.fa HiC.filtered.bam nchrs --quick_view --gfa "hap1.p_ctg.gfa,hap2.p_ctg.gfa"
```



## <span id="juicebox">Juicebox curation</span>

There are two ways of generating `.assembly` and `.hic` files for visualization and manual curation in Juicebox. You can choose one of them according to your preference, **though we recommend the second method**.

#### (1) SALSA-style `scaffolds.agp`

First, install the dependencies, including (1) [3D-DNA](https://github.com/aidenlab/3d-dna), (2) [matlock](https://github.com/phasegenomics/matlock), (3) [Juicebox scripts](https://github.com/phasegenomics/juicebox_scripts). Then, generate the `.assembly` and `.hic` files by following these steps:

```bash
# (1) Generate .mnd file
$ /path/to/matlock bam2 juicer HiC.filtered.bam out.links.mnd
$ sort -k2,2 -k6,6 out.links.mnd > out.sorted.links.mnd

# (2) Generate .assembly file
$ /path/to/juicebox_scripts/agp2assembly.py scaffolds.agp scaffolds.assembly

# (3) Generate .hic file
$ bash /path/to/3d-dna/visualize/run-assembly-visualizer.sh -p false scaffolds.assembly out.sorted.links.mnd
```

> [!NOTE]
>
> If there are any contigs corrected by HapHiC, you need to re-align Hi-C reads to `corrected_asm.fa` and re-filter them instead of using the original `HiC.filtered.bam`. Otherwise, there will not be any Hi-C signals on the corrected contigs in Juicebox. This is because that the IDs of corrected contigs in the SALSA-style `scaffolds.agp` do not match the contig IDs in the original BAM file. This is also the reason why the method is not recommended.

You can recall these steps on the command line:

```bash
$ /path/to/HapHiC/haphic juicer
```

After manual curation in Juicebox, you will obtain the modified assembly file `scaffolds.review.assembly`. To generate the final FASTA file for the scaffolds:

```bash
# Generate the final FASTA file for the scaffolds
$ /path/to/juicebox_scripts/juicebox_assembly_converter.py -a scaffolds.review.assembly -f asm.fa -s
```

#### (2) YaHS-style `scaffolds.raw.agp` (recommended)

To avoid the necessity of re-aligning Hi-C data, we have incorporated a YaHS-style `scaffolds.raw.agp` since HapHiC version 1.0.1. In this AGP file, the broken contigs are not assigned new IDs. Instead, their starting and ending coordinates in the raw contigs are displayed in the seventh and eighth columns. By following [the approach provided by YaHS](https://github.com/c-zhou/yahs#manual-curation-with-juicebox-jbat), you can generate the `.assembly` and `.hic` files without the need for re-aligning.

After constructing the final scaffolds (Step 4), HapHiC automatically generates a shell script for [visualization and curation in Juicebox](#juicebox). Ensure that [Java](https://openjdk.org/install/) and [samtools](https://github.com/samtools/samtools) have been installed and added to `$PATH` on your system. Then, run the following command:

```bash
# bash, not sh
$ bash juicebox.sh
```

> [!NOTE]
>
> * Since HapHiC version 1.0.6, there is no longer a need to manually set the scale factor in Juicebox. In addition, the saved `.review.assembly` file can now be parsed correctly by Juicebox.
> * For large genomes, it is necessary to adjust the memory settings in the juicebox.sh file for Java (e.g. set to `-Xmx64G` or higher) to avoid out-of-memory errors or to improve the execution speed.

After manual curation in Juicebox, you will obtain the modified assembly file `out_JBAT.review.assembly`. To generate the final FASTA file for the scaffolds:

```bash
# Generate the final FASTA file for the scaffolds
$ /path/to/HapHiC/utils/juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp asm.fa
```



## <span id="visualization">Visualization</span>

![](./images/HapHiC2.png)

Since HapHiC version 1.0.2, we have introduced a `haphic plot` command to generate highly customizable Hi-C contact maps. This command requires two input files, a filtered BAM file `HiC.filtered.bam` and a scaffold AGP file containing contig IDs that match those in the BAM file:

```bash
# For HapHiC scaffolding result
$ /path/to/HapHiC/haphic plot scaffolds.raw.agp HiC.filtered.bam
# For the AGP file generated after manual curation in Juicebox
$ /path/to/HapHiC/haphic plot out_JBAT.FINAL.agp HiC.filtered.bam
```

The visualized Hi-C contact map is output as `contact_map.pdf`. This process may be somewhat slow if the BAM file is large, taking several minutes per 10 GiB of the BAM file. Upon completion, the program will produce a binary file named `contact_matrix.pkl`. This file can be utilized in place of `HiC.filtered.bam` for faster visualization (only ~1 minute), facilitating fine-tuning of heatmap-style related parameters:

```bash
# Use previously generated `contact_matrix.pkl` for faster visualization
$ /path/to/HapHiC/haphic plot out_JBAT.FINAL.agp contact_matrix.pkl
```

> [!NOTE]
>
> When using `contact_matrix.pkl` for faster visualization, the input AGP file and the parameters `--bin_size` and `--min_len` must remain consistent throughout.

By default, the bin size is set to 500 Kbp and only scaffolds exceeding 1 Mbp in length will be displayed on the contact map. To modify these parameters:

```bash
# Set the bin size to 1 Mbp and display only scaffolds longer than 5 Mbp
$ /path/to/HapHiC/haphic plot out_JBAT.FINAL.agp HiC.filtered.bam --bin_size 1000 --min_len 5
```

Additionally, you can create `separate_plots.pdf`, which illustrates the contact map for each scaffold individually:

```bash
$ /path/to/HapHiC/haphic plot out_JBAT.FINAL.agp HiC.filtered.bam --separate_plots
```

To change the colormap, origin, border style, and normalization method for the contact maps, refer to the examples provided in the figure above.

**Do you think these contact maps look cool? This function can also visualize results from other scaffolders!** You only need to prepare a BAM file ([by mapping and filtering Hi-C reads](#align-hi-c-data-to-the-assembly)) for your chromosome-level FASTA file and create a corresponding AGP file:

```bash
# This script generates an AGP file for your FASTA file
$ /path/to/HapHiC/utils/mock_agp_file.py chr_asm.fa > chr_asm.agp
# Then, you can visualize your results using `haphic plot` with the BAM file and the AGP file
```


## <span id="refsort">Order and orient whole scaffolds using a reference genome</span>

HapHiC has introduced a separate command, `haphic refsort`, in version 1.0.4 to order and orient whole scaffolds according to a reference genome.

To begin, you should prepare a PAF file by align raw contigs (**not scaffolds**) to a chromosome-level reference genome using [minimap2](https://github.com/lh3/minimap2). The reference genome can be from the same species or a closely related one:

```bash
# The preset can be `asm5` if the reference genome is well-assembled from the same species
$ minimap2 -x asm20 ref.fa asm.fa --secondary=no -t 28 -o asm_to_ref.paf
# `haphic refsort` can also be compatible with other aligners, like wfmash
$ wfmash ref.fa asm.fa -m -n 1 -S 1 -t 28 | cut -f 1-6,8- > asm_to_ref.paf
```

By using the `haphic refsort` command, you can generate new AGP and FASTA files based on the PAF file:

```bash
# By default, scaffolds are output based on the alphabetical order of the chromosome IDs of the reference genome
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf > scaffolds.refsort.agp
# You can specify the order by listing chromosome IDs of the reference genome separated by commas (no spaces)
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf --ref_order "chr1,chr2,chr3,chr4,..." > scaffolds.refsort.agp
# If you want to generate a new FASTA file (default name: `scaffolds.refsort.fa`) as well
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf --fasta asm.fa > scaffolds.refsort.agp
# If you want to run `haphic refsort` on manually curated `out_JBAT.FINAL.agp`
$ haphic refsort out_JBAT.FINAL.agp asm_to_ref.paf > scaffolds.refsort.agp
```

The generated `scaffolds.refsort.agp` file can be directly used for [Juicebox curation](#juicebox) and for `haphic plot` [visualization](#visualization). 

>  [!NOTE]
>
> Please note that **this function is NOT reference-based scaffolding and will NOT alter your scaffolds**, it only changes the way of presentation through overall ordering and orientation of the entire scaffolds. 

Here is an example of the autotetraploid sugarcane Np-X assembly:

![](./images/refsort_example.png)



## <span id="faqs">Frequently asked questions (FAQs)</span>

* **What can I do when the anchoring rate is too low?**

  There are three parameters controlling the anchoring rate through the reassignment step: `--min_RE_sites`, `--min_links`, and `--min_link_density`. By default, these parameters are set to 25, 25, and 0.0001, respectively. However, both the contig contiguity and Hi-C sequencing depth vary across different projects. By checking the `*statistics.txt` files in `01.cluster/inflation_*`, you can find better values for these parameters to get a scaffolding result with a higher anchoring rate.

  For small genomes, the default `--Nx 80` in the clustering step and the default `--min_group_len` in reassignment step may also negatively affect the anchoring rate. To address this, you can increase the value of `--Nx` and decrease `--min_group_len`, or even disable these two functions entirely by using `--Nx 100` and `--min_group_len 0`.

* **How to run HapHiC if I don't know the exact number of chromosomes?**

  You could try [quick view](#quick_view). In this mode, HapHiC ignores the `nchrs` parameter (you can fill in any integer), and orders and orients contigs without clustering (similar to `*.0.hic` in 3D-DNA). After visualizing the results in Juicebox, you can count the number of chromosomes based on the Hi-C contact map, and rerun HapHiC pipeline with this number. Alternatively, you can manually curate the assembly and split chromosomes in Juicebox by yourself.

* **What can I do when I see "It seems that some chromosomes were grouped together" or "The maximum number of clusters is even less than the expected number of chromosomes" in the clustering step?**

  The question is complicated. HapHiC recommends a "best" inflation parameter based on the `nchrs` you specified and the distribution of group lengths. Several factors could cause this problem.

  **(1)** The `--max_inflation` parameter is too low. Try increasing it. The default value of `--max_inflation 3.0` is generally enough. However, sometimes the best inflations are even greater than 7.0 due to a higher background of Hi-C links between chromosomes (We encountered this in the [taro genome](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13239)). This could be due to biological specificity or low quality of the Hi-C library. In rare cases, if a higher `--max_inflation` still doesn't work, try using [quick view](#quick_view) and manually splitting chromosomes in Juicebox.

  **(2)** Some homologous chromosomes may be grouped together due to assembly errors. This is common in scaffolding phased assemblies. Please use more aggressive parameters to correct contigs, filter out contigs, or use `--moreve_allelic_links` to remove Hi-C links between homologous chromosomes. If you are unsure about the type or the proportion of assembly errors, a [quick view](#quick_view) may be helpful.

  **(3)** When the chromosome lengths vary greatly, HapHiC may mistake a large chromosome for two or more chromosomes clustered together. In this case, you can choose a reasonable inflation and [run the remaining steps individually](#step_by_step).



## <span id="problems">Problems and bug reports</span>

* **Issues:** https://github.com/zengxiaofei/HapHiC/issues
* **Before asking questions, please read:** [Important: More information allows me to grasp your issue faster](https://github.com/zengxiaofei/HapHiC/issues/32)
* **We have built a [DeepSeek-Powered HapHiC Knowledge Base](https://github.com/zengxiaofei/HapHiC/issues/106) (WeChat registration is required)**




## <span id="citing">Citing HapHiC</span>

If you have used HapHiC in your work, please cite our paper published in Nature Plants:

> Xiaofei Zeng, Zili Yi, Xingtan Zhang, Yuhui Du, Yu Li, Zhiqing Zhou, Sijie Chen, Huijie Zhao, Sai Yang, Yibin Wang, Guoan Chen. Chromosome-level scaffolding of haplotype-resolved assemblies using Hi-C data without reference genomes. *Nature Plants*, 10:1184-1200. doi: [https://doi.org/10.1038/s41477-024-01755-3](https://doi.org/10.1038/s41477-024-01755-3)

There is also a Research Briefing available in Nature Plants:

> Xiaofei Zeng, Guoan Chen. (2024) Achieving de novo scaffolding of chromosome-level haplotypes using Hi-C data. *Nature Plants*, 10:1157-1158. doi: [https://doi.org/10.1038/s41477-024-01756-2](https://doi.org/10.1038/s41477-024-01756-2)

If you have used the optimization function for contig ordering and orientation, please cite ALLHiC as well:

> Xingtan Zhang, Shengcheng Zhang, Qian Zhao, Ray Ming, Haibao Tang. (2019) Assembly of allele-aware, chromosomal-scale autopolyploid genomes based on Hi-C data. *Nature Plants*, 5:833-845. doi: [https://doi.org/10.1038/s41477-019-0487-8](https://doi.org/10.1038/s41477-019-0487-8)




## <span id="reproduce">Reproducibility</span>

To reproduce the results in our paper, please use the HapHiC code in commit [431b7b6](https://github.com/zengxiaofei/HapHiC/tree/431b7b6f3e0471c960985ecdbbab4d18b452a22f).



## Star history

Love HapHiC? A star ⭐️ would be a huge encouragement for us and help others discover the project.

[![Star History Chart](https://api.star-history.com/svg?repos=zengxiaofei/HapHiC&type=Date)](https://star-history.com/#zengxiaofei/HapHiC&Date)
