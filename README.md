## HapHiC: a fast, reference-independent, allele-aware scaffolding tool based on Hi-C data



![](./images/HapHiC1.png)

HapHiC is an allele-aware scaffolding tool that uses Hi-C data to scaffold haplotype-phased genome assemblies into chromosome-scale pseudomolecules. Unlike [ALLHiC](https://github.com/tangerzhang/ALLHiC), another allele-aware scaffolder, HapHiC can achieve this without the need for reference genomes. Our evaluations indicate that HapHiC outperforms other Hi-C scaffolding tools with higher tolerance to low contig N50, low Hi-C sequencing depth, and various types of assembly errors. Additionally, HapHiC is super-fast and also suitable for haplotype-collapsed diploid and allopolyploid genome assemblies.

**Features:**

- [x] Chromosome-level scaffolding of haplotype-phased assemblies without reference genomes
- [x] Efficient correction of chimeric contigs (misjoins) with little impact on contig N50
- [x] Much higher tolerance to chimeric contigs, collapsed contigs, and switch errors
- [x] Improved performance in chromosome assignment of contigs
- [x] Improved performance in ordering and orienation of contigs
- [x] Super-fast and memory-efficient
- [x] Able to order and orient contigs without prior knowledge of the number of chromosomes
- [x] Able to utilize phasing information from hifiasm with varying confidence levels 

**Recent updates:**

* Version 1.0.1 (2023.11.30): Improved AGP output by incorporating a YaHS-style `scaffolds.raw.agp`  for compatibility with the Juicebox visualization method suggested by YaHS.

**Terminology:** To ensure concision and clarity, we use the term "contigs" to refer to the fragmented genome sequences in the input assembly, although they could be either contigs or scaffolds in actuality.

## Table of contents

- [Installation](#installation)
- [Quick start](#quick_start)
  * [Align Hi-C data to the assembly](#align)
  * [Run HapHiC scaffolding pipeline](#pipeline)
- [Go through the pipeline step by step](#step_by_step)
  * [[Step 1] Clustering](#step1)
  * [[Step 2] Reassignment](#step2)
  * [[Step 3] Ordering and orientation](#step3)
  * [[Step 4] Building pseudomolecules](#step4)
- [Work with hifiasm (experimental)](#hifiasm)
- [Quick view mode](#quick_view)
- [Juicebox visualization and curation](#juicebox)
- [Frequently asked questions (FAQs)](#faqs)
- [Problems and bug reports](#problems)
- [Citing HapHiC](#citing)
- [Reproducibility](#reproduce)

## <span id="installation">Installation</span>

HapHiC has been tested and validated on servers running Linux, equipped with either Intel Xeon or AMD EPYC CPUs.

```bash
# (1) Download HapHiC from GitHub
$ git clone https://github.com/zengxiaofei/HapHiC.git

# (2) Resolve dependencies
# We strongly recommend using conda to install dependencies. If you prefer manual installation, refer to HapHiC/conda_env/create_conda_env_py310.sh
$ conda env create -f HapHiC/conda_env/environment_py310.yml
# Activate the HapHiC conda environment
$ conda activate haphic # or: source /path/to/conda/bin/activate haphic

# (3) Show all available commands and help message
$ /path/to/HapHiC/haphic -h
```



## <span id="quick_start">Quick start</span>

### <span id="align">Align Hi-C data to the assembly</span>

First, you need to prepare a BAM file by aligning Hi-C data to the assembly. Here is the way that we recommend:

```bash
# (1) Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
$ bwa index asm.fa
$ bwa mem -5SP asm.fa /path/to/read1_fq.gz /path/to/read2_fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o HiC.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
$ /path/to/HapHiC/utils/filter_bam.py HiC.bam 1 --NM 3 --threads 14 | samtools view - -b -@ 14 -o HiC.filtered.bam
```

**Notes:**

* Here, `asm.fa` can be haplotype-collapsed contigs (e.g., `p_ctg` in hifiasm), haplotype-phased unitigs (e.g., `p_utg` in hifiasm), or one or more sets of haplotype-resolved contigs (e.g., `hap*.p_utg` in hifiasm). In addition, `asm.fa` may also be scaffolds output by other scaffolders.
* You can prepare the BAM file according to your own preferences or requirements, but **DO NOT** sort it by coordinate. If your BAM file is already sorted by coordinate, you need to resort it by read name ( `samtools sort -n` ).
* We **DO NOT** recommend the Juicer pipeline for Hi-C reads alignment, particularly in haplotype-phased assemblies.

### <span id="pipeline">Run HapHiC scaffolding pipeline</span>

**(i) One-line command**. HapHiC provides a one-line command `haphic pipeline` to execute the entire scaffolding pipeline. The required parameters are 1) `asm.fa` , your genome assembly file in FASTA format; 2) `HiC.filtered.bam` , the BAM file prepared in the previous step; 3) `nchrs` , the number of chromosomes present in the assembly, and also the expected number of output scaffolds.

```bash
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs
```

**(ii) Restriction site.** The default restriction site is `GATC`  (MboI/DpnII). You can modify this  using the `--RE` parameter. If you are unsure or if your Hi-C library was constructed without restriction enzymes (REs), it is acceptable to leave it as the default.

```bash
# For HindIII
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --RE "AAGCTT"
# For Arima two-enzyme chemistry
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --RE "GATC,GANTC"
# For Arima four-enzyme chemistry
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --RE "GATC,GANTC,CTNAG,TTAA"
```

**(iii) Contig correction.** To correct input contigs based on Hi-C linking information, use `--correct_nrounds` to enable assembly correction and set the number of correction rounds. For example:

```bash
# Typically, two rounds of assembly correction are enough
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --correct_nrounds 2
```

**(iv) Switch error.** If your input assembly is haplotype-phased and has a high switch error rate (often introduced by assemblers when the sequence divergence between haplotypes is very low), use `--remove_allelic_links` to remove Hi-C links between allelic contigs. The value should be the ploidy of the assembly. For example:

```bash
# For haplotype-phased assembles of autotetraploids, set the parameter to 4
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --remove_allelic_links 4
```

**(v) Performance.** Use `--threads` to set the number of threads for BAM file reading, and `--processes` to create multiple processes for contig ordering and orientation. For example:

```bash
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --threads 8 --processes 8
```

**Parameters**

For more information, run `haphic pipeline --help` .

**Final outputs**

* `01.cluster/corrected_asm.fa` : The corrected assembly in FASTA format. This file is generated only when assembly correction is enabled.
* `04.build/scaffolds.agp` : A [SALSA](https://github.com/marbl/SALSA)-style [AGP file](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) containing information about scaffold assignment, ordering and orientation information for all sequences in `corrected_asm.fa`. If there are any chimeric contigs that are corrected by HapHiC, the broken contigs will be assigned new IDs.
* `04.build/scaffolds.raw.agp` : A [YaHS](https://github.com/c-zhou/yahs)-style AGP file containing information about scaffold assignment, ordering and orientation information for all sequences in `asm.fa`. The broken contigs will not be assigned new IDs, but their starting and ending coordinates in the raw contigs will be displayed in the seventh and eighth columns.
* `04.build/scaffolds.fa` : The final scaffolds in FASTA format.
* `04.build/juicebox.sh` : A shell script for [Juicebox visualization and curation](#juicebox).

**Note:** Although the one-line command is convenient, the automatic parameter tuning may fail, leading to poor results or even a pipeline interruption in rare cases. If this occurs, we recommend [running each step individually with manual parameter tuning](#step_by_step) or trying the [quick view mode](#quick_view) described below.



## <span id="step_by_step">Go through the pipeline step by step</span>

### <span id="step1">[Step 1]. Clustering</span>

Before clustering, HapHiC performs preprocessing to correct assembly misjoins, filter out short, mis-assembled contigs, and remove allelic Hi-C links. After that, a Markov cluster algorithm (MCL algorithm) is used to cluster the contigs into groups. Unlike agglomerative hierarchical clustering (AHC, used in LACHESIS and ALLHiC), which specifies the number of clusters, the MCL Algorithm implicitly controls it with a parameter called "inflation". The higher the inflation, the more the groups are clustered. The main problem with AHC is that even though the number of clusters is specified, contigs from different chromosomes may also be clustered into the same group. This is common in phased diploid or polyploid genome assemblies. To solve this, HapHiC tries a series of inflations to cluster the contigs (controlled by `min_inflation` and `max_inflation` ) and recommends a "best" one based on both the expected number of chromosomes `nchrs` provided and the length distribution of the groups.

```bash
$ /path/to/HapHiC/haphic cluster asm.fa HiC.filtered.bam nchrs
```

**Parameters**

For more information, run `haphic cluster --help` .

**Main outputs**

* `corrected_asm.fa` : The corrected assembly in FASTA format. This file is generated only when assembly correction is enabled.

* `corrected_ctgs.txt` : A text file listing the IDs of all corrected contigs.

* `full_links.pkl` : A binary file that stores the number of Hi-C links between each contig pair.

* `HT_links.pkl` : A binary file that records the number of Hi-C links between each half of contig pairs.

* `paired_links.clm` : A text file recording the positional information of paired Hi-C links.

* `inflation_*` : Output directories for respective inflations.

  ├── `group*.txt` : Files containing the contigs and their basic information (including lengths and numbers of restriction sites) for each group, also reffered to as `counts_RE.txt` in ALLHiC.

  └── `mcl_inflation_*.clusters.txt` : Markov clustering results.

**The "best" inflation**

You can find the "best" inflation recommendation in the log file `HapHiC_cluster.log` , like:

```
2022-11-07 17:50:08 <HapHiC_cluster.py> [recommend_inflation] You could try inflation from 1.20 (length ratio = 0.75)
```

In some cases, HapHiC cannot get the "best" one. It could be due to inappropriate parameters or extensive assembly errors. Check whether the parameters used are correct / appropriate and then try to tune the parameters for assembly correction, contig / Hi-C link filtration, or Markov Clustering:

``````
2022-11-19 13:20:38 <HapHiC_cluster.py> [recommend_inflation] It seems that some chromosomes were grouped together (length ratio = 0.5)...
``````

### <span id="step2">[Step 2] Reassignment</span>

In the previous step, some contigs may have been filtered out before clustering or assigned to incorrect groups. Additionally, the number of final clusters output by Markov clustering may exceed the specified number of chromosomes ( `nchrs` ). To address these issues, we added a reassignment step to rescue contigs that are not in any groups, reassign contigs to the correct groups, and perform an additional agglomerative hierarchical clustering to concatenate groups if necessary. The input files `full_links.pkl` , `mcl_inflation_x.clusters.txt ` , and `paired_links.clm` are outputs from the clustering step, where `x` represents the "best" inflation value:

```bash
$ /path/to/HapHiC/haphic reassign asm.fa full_links.pkl mcl_inflation_x.clusters.txt paired_links.clm --nclusters nchrs
```

**Note:**  If assembly correction has been performed, use `corrected_asm.fa` as input FASTA file instead of `asm.fa`.

**Parameters**

For more information, run `haphic reassign --help` .

**Main outputs**

* `final_groups/group*.txt` : Files containing the contigs and their basic information for each final group after reassignment.
* `final_groups/final_cluster.txt` : The final clustering result.
* `split_clms/` : A directory containing group-specific CLM files.

### <span id="step3">[Step 3] Ordering and orientation</span>

The ordering and orientation step in HapHiC is implemented using an integration of algorithms from [3D-DNA](https://github.com/aidenlab/3d-dna) and [ALLHiC](https://github.com/tanghaibao/allhic). First, an efficiency-improved 3D-DNA iterative scaffolding algorithm (refered to as "fast sorting") is used to quickly order and orient the contigs. Then, the ordering and orientation of contigs are input as an initial configuration and optimized with the ALLHiC program ([a modified version](http://github.com/zengxiaofei/allhic), in which the hot-start optimization has been fixed). The input file `HT_links.pkl` is the output file from the clustering step; the directory `split_clms` and the group files `final_groups/group*.txt` were created in the reassignment step. The optional parameter `--processes` is used to set the number of processes for the ordering and orientation.

```bash
$ /path/to/HapHiC/haphic sort asm.fa HT_links.pkl split_clms final_groups/group*.txt --processes 8
```

**Note:**  If assembly correction has been performed, use `corrected_asm.fa` as input FASTA file instead of `asm.fa`.

**Parameters**

For more information, run `haphic sort --help` .

**Main outputs**

* `group*.tour.sav` : The fast sorting result of contigs within each group.
* `group*.tour` : The final contig ordering and orientation result for each group after ALLHiC optimization.

### <span id="step4">[Step 4] Building pseudomolecules</span>

The final step is to build the scaffolds (pseudomolecules) using the chromosome assignment, ordering and orientation information of contigs from the `group*.tour` files. By default, the output scaffolds are sorted by scaffold length.

**If assembly correction was not performed**:

```bash
$ /path/to/HapHiC/haphic build asm.fa asm.fa HiC.filtered.bam group*.tour
```

**If assembly correction has been performed**, use `corrected_asm.fa` as input FASTA file instead of **the first** `asm.fa` . Additionally, specify the corrected contig list `corrected_ctgs.txt` using the `--corrected_ctgs` parameter. Otherwise, the YaHS-style `scaffolds.raw.agp` generated may be incorrect.

```bash
$ /path/to/HapHiC/haphic build corrected_asm.fa asm.fa HiC.filtered.bam group*.tour --corrected_ctgs corrected_ctgs.txt
```

**Note:**  

* The second `asm.fa` (raw uncorrected assembly) and `HiC.filtered.bam` are required since HapHiC version 1.0.1 for generating the script for juicebox visualization and curation.

**Parameters**

For more information, run `haphic build --help` .

**Main outputs**

* `scaffolds.agp` : A SALSA-style AGP file containing information about scaffold assignment, ordering and orientation information for all sequences in `corrected_asm.fa`. If there are any chimeric contigs that are corrected by HapHiC, the broken contigs will be assigned new IDs.
* `scaffolds.raw.agp` : A YaHS-style AGP file containing information about scaffold assignment, ordering and orientation information for all sequences in `asm.fa`. The broken contigs will not be assigned new IDs, but their starting and ending coordinates in the raw contigs will be displayed in the seventh and eighth columns.
* `scaffolds.fa` : The final scaffolds in FASTA format.
* `juicebox.sh` : A shell script for [Juicebox visualization and curation](#juicebox).



## <span id="hifiasm">Work with hifiasm (experimental)</span>

When scaffolding a phased [hifiasm](https://github.com/chhylp123/hifiasm) assembly, you can run HapHiC with the GFA file(s) output by hifiasm. Here, the term "phased hifiasm assembly" refers to the haplotype-resolved primary contigs assembled via the trio binning or Hi-C-based algorithm ( `*.hap*.p_ctg.gfa` ), as well as the phased unitigs ( `*.p_utg.gfa` ). 

HapHiC uses the read depth information in the GFA file(s) to filter out potential collapsed contigs/unitigs before clustering.  If more than one GFA file is provided, HapHiC assumes these GFA files are haplotype-specific ( `*.hap*.p_ctg.gfa` ), and artificially removes or reduces the Hi-C links between the haplotypes according to this phasing information. Note that the contigs/unitigs in GFA file(s) should match those in FASTA file. Either `.gfa`  or `noseq.gfa` is acceptable.

```shell
# (1) For hifiasm primary unitigs, use the GFA file to filter out potential collapsed unitigs before clustering
$ /path/to/HapHiC/haphic pipeline p_utg.fa HiC.filtered.bam nchrs --gfa p_utg.gfa

# (2) In addition to read depth filtering, HapHiC can also reduce Hi-C links between haplotypes according to phasing information in GFA files for haplotype-resolved primary contigs

# By default, all Hi-C links between haplotypes are completely removed and contigs from different haplotypes will not be clustered into the same group
$ /path/to/HapHiC/haphic pipeline allhaps.fa HiC.filtered.bam nchrs --gfa "hap1.p_ctg.gfa,hap2.p_ctg.gfa"

# The weight can be set to 0 to ignore the phasing information. You can also set it between 0 and 1 to run HapHiC as a double check. In the latter case, contigs from different haplotypes might be clustered together
$ /path/to/HapHiC/haphic pipeline allhaps.fa HiC.filtered.bam nchrs --gfa "hap1.p_ctg.gfa,hap2.p_ctg.gfa" ---phasing_weight 0
```



## <span id="quick_view">Quick view mode</span>

You can try the quick view mode in HapHiC when: 
1. The exact number of chromosomes is unknown.
2. HapHiC cannot provide an acceptable clustering result or encounters a pipeline interruption.
3. You need a quick view of your assembly (e.g., to identify the type and approximate proportion of assembly errors).
4. You just want to manually curate your assembly and split chromosomes in Juicebox by yourself. 

In quick view mode, HapHiC simply uses the fast sorting to order and orient all contigs without clustering. The result is similar to `*.0.hic` in [3D-DNA](https://github.com/aidenlab/3d-dna). Most parameters are disabled in this mode, but you can use `--correct_nrounds` to correct input contigs. When scaffolding a haplotype-resolved hifiasm assembly ( `*.hap*.p_ctg.gfa` ), you can still partition contigs into different haplotypes with the haplotype-specific GFA files.

```shell
# Autohic will ignore the parameter "nchrs", it can be any integer
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --quick_view
# Correct input contigs before a quick view
$ /path/to/HapHiC/haphic pipeline asm.fa HiC.filtered.bam nchrs --quick_view --correct_nrounds 2
# Partition contigs into different haplotypes in quick view mode
$ /path/to/HapHiC/haphic pipeline allhaps.fa HiC.filtered.bam nchrs --quick_view --gfa "hap1.p_ctg.gfa,hap2.p_ctg.gfa"
```



## <span id="juicebox">Juicebox visualization and curation</span>

There are two ways of generating `.assembly` and `.hic` files for visualization and manual curation in Juicebox. You can choose one of them according to your preference.

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

**Note:** If there are any contigs corrected by HapHiC, you need to re-align Hi-C reads to `corrected_asm.fa` and re-filter them instead of using the original `HiC.filtered.bam` . Otherwise, there will not be any Hi-C signals on the corrected contigs in Juicebox. This is because that the IDs of corrected contigs in the SALSA-style `scaffolds.agp` do not match the contig IDs in the original BAM file.

You can recall these steps on the command line:

```bash
$ /path/to/HapHiC/haphic juicer
```

#### (2) YaHS-style `scaffolds.raw.agp` (recommended)

To avoid the necessity of re-aligning Hi-C data, we have incorporated a YaHS-style `scaffolds.raw.agp` since HapHiC version 1.0.1. In this AGP file, the broken contigs are not assigned new IDs. Instead, their starting and ending coordinates in the raw contigs are displayed in the seventh and eighth columns. By following [the approach provided by YaHS](https://github.com/c-zhou/yahs#manual-curation-with-juicebox-jbat), you can generate the `.assembly` and `.hic` files without the need for re-aligning.

After constructing the final scaffolds, HapHiC automatically generates a shell script for visualization and curation in Juicebox. Ensure that [Java](https://openjdk.org/install/) and [samtools](https://github.com/samtools/samtools) have been installed and added to `$PATH` on your system. Then, run the following command:

```bash
$ bash juicebox.sh
```

**Note:** In the output log file `out_JBAT.log` , you can find the corresponding scale factor, e.g., `[I::main_pre] scale factor: 2` . To ensure proper alignment of Hi-C contact maps with the boundaries of scaffolds and superscaffolds in Juicebox, please set your own scale factor in Juicebox through the menu `Assembly > Set Scale` .

## <span id="faqs">Frequently asked questions (FAQs)</span>

* **How can I do when the anchoring rate is too low?**

  There are three parameters controlling the anchoring rate through the reassignment step: `--min_RE_sites` , `--min_links`  , and `--min_link_density` . By default, these parameters are set to 25, 25, and 0.0001, respectively. However, both the contig contiguity and Hi-C sequencing depth vary across different projects. By checking the `*statistics.txt` files in `01.cluster/inflation_*` , you can find better values for these parameters to get a scaffolding result with a higher anchoring rate.

  For small genomes, the default `--Nx 80` in the clustering step and the default `--min_group_len` in reassignment step may also negatively affect the anchoring rate. To address this, you can increase the value of `--Nx` and decrease `--min_group_len`, or even disable these two functions entirely by using `--Nx 100` and `--min_group_len 0`.

* **How to run HapHiC if I don't know the exact number of chromosomes?**

  You could try [quick view](#quick_view). In this mode, HapHiC ignores the `nchrs` parameter (you can fill in any integer), and scaffold contigs without clustering (similar to `*.0.hic` in 3D-DNA). After visualizing the results in Juicebox, you can count the number of chromosomes based on the Hi-C contact map, and rerun HapHiC pipeline with this number. Alternatively, you can manually curate the assembly and split chromosomes in Juicebox by yourself.

* **How can I do when I see "It seems that some chromosomes were grouped together" in the clustering step?**

  The question is complicated. HapHiC recommends a "best" inflation parameter based on the `nchrs` you specified and the distribution of group lengths. Several factors could cause this problem.

  **(1)** The `--max_inflation` parameter is too low. Try increasing it. The default value of `--max_inflation 3.0` is generally enough. However, sometimes the best inflations are even greater than 7.0 due to a higher background of Hi-C links between chromosomes (We encountered this in the [taro genome](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13239)). This could be due to biological specificity or low quality of the Hi-C library. In rare cases, if a higher `--max_inflation` still doesn't work, try using [quick view](#quick_view) and manually splitting chromosomes in Juicebox.

  **(2)** Some homologous chromosomes may be grouped together due to assembly errors. This is common in scaffolding phased assemblies. Please use more aggressive parameters to correct contigs, filter out contigs, or remove Hi-C links between homologous chromosomes. If you are unsure about the type or the proportion of assembly errors, a [quick view](#quick_view) may be helpful.

  **(3)** When the chromosome lengths vary greatly, HapHiC may mistake a large chromosome for two or more chromosomes clustered together. In this case, you can choose a reasonable inflation and [run the remaining steps individually](#step_by_step).



## <span id="problems">Problems and bug reports</span>

* **Issues:** https://github.com/zengxiaofei/HapHiC/issues




## <span id="citing">Citing HapHiC</span>

If you use HapHiC in your work, please cite our preprint on bioRxiv:

> Xiaofei Zeng, Zili Yi, Xingtan Zhang, Yuhui Du, Yu Li, Zhiqing Zhou, Sijie Chen, Huijie Zhao, Sai Yang, Yibin Wang, Guoan Chen. (2023) Chromosome-level scaffolding of haplotype-resolved assemblies using Hi-C data without reference genomes. *bioRxiv*, 2023.11.18.567668. doi: [https://doi.org/10.1101/2023.11.18.567668](https://doi.org/10.1101/2023.11.18.567668)

If you use the optimization function for contig ordering and orientation (by default), please also cite ALLHiC:

> Xingtan Zhang, Shengcheng Zhang, Qian Zhao, Ray Ming, Haibao Tang. (2019) Assembly of allele-aware, chromosomal-scale autopolyploid genomes based on Hi-C data. *Nature Plants*, 5:833-845. doi: [https://doi.org/10.1038/s41477-019-0487-8](https://doi.org/10.1038/s41477-019-0487-8)




## <span id="reproduce">Reproducibility</span>

To reproduce the results in our paper, please use the HapHiC code in commit [431b7b6](https://github.com/zengxiaofei/HapHiC/tree/431b7b6f3e0471c960985ecdbbab4d18b452a22f).
