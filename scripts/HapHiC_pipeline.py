#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-10-24 16:27

import os
import sys
import subprocess
import argparse
import logging
import re
import glob
import time

import HapHiC_cluster
import HapHiC_reassign
# import HapHiC_sort
import HapHiC_build

from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

abspath = os.path.abspath
pathjoin = os.path.join
pathsplit = os.path.split

def parse_argument():

    parser = argparse.ArgumentParser(prog='haphic pipeline')

    # Parameters for parsing input files
    input_group = parser.add_argument_group('>>> Parameters for parsing input files')
    input_group.add_argument(
            'fasta', help='draft genome in FASTA format')
    input_group.add_argument(
            'alignments', help='filtered Hi-C read alignments in BAM/pairs format (DO NOT sort it by coordinate)')
    input_group.add_argument(
            'nchrs', type=int, help='expected number of chromosomes')
    input_group.add_argument(
            '--aln_format', choices={'bam', 'pairs', 'bgzipped_pairs', 'auto'}, default='auto', help='file format for Hi-C read alignments, default: %(default)s')
    input_group.add_argument(
            '--gfa', default=None,
            help='(experimental) matched GFA file(s) of the phased hifiasm assembly, separated with commas (e.g., `--gfa "gfa1,gfa2"`), default: %(default)s. '
            'Here, the "phased hifiasm assemblies" include the haplotype-resolved primary contigs assembled using the trio binning or Hi-C-based algorithm '
            '(`*.hap*.p_ctg.gfa`) and the phased unitigs (`*.p_utg.gfa`). HapHiC uses the read depth information in each GFA file to filter out contigs '
            '(see `--read_depth_upper`). If multiple GFA files are provided, HapHiC assumes these GFA files are haplotype-specific, and then artificially '
            'reduces the Hi-C links between the haplotypes according to this phasing information (see `--phasing_weight`). This parameter can also work with '
            '`--quick_view` to separate contigs into different haplotype groups')
    input_group.add_argument(
            '--ul', default=None,
            help='ultra-long read alignments in BAM format, default: %(default)s')

    # Paramters for pipeline control and shared parameters
    pipe_group = parser.add_argument_group('>>> Paramters for pipeline control and shared parameters')
    pipe_group.add_argument(
            '--RE', default='GATC',
            help='restriction enzyme site(s) (e.g., GATC for MboI & AAGCTT for HindIII), default: %(default)s. If more than one enzyme was used '
            'in the Hi-C library construction, such as with the Arima genomics kit, separate the RE sites with commas (e.g., `--RE "GATC,GANTC"` for Arima '
            'two-enzyme chemistry and `--RE "GATC,GANTC,CTNAG,TTAA"` for Arima four-enzyme chemistry)')
    pipe_group.add_argument(
            '--steps', default='1,2,3,4',
            help='run specified HapHiC pipeline steps, separated with commas, default: %(default)s. `steps` should be continuous and start with "1" (e.g., '
            '--steps "1,2,3"). (1: cluster, 2: reassign, 3: sort, 4: build)')
    pipe_group.add_argument(
            '--quick_view', default=False, action='store_true',
            help='in quick view mode, HapHiC will skip the clustering and reassignment steps, and order and orient all contigs with fast sorting, default: %(default)s. '
            'This is helpful when you encounter difficulties in generating ideal clusters or when you are unsure of the exact number of chromosomes')
    pipe_group.add_argument(
            '--normalize_by_nlinks', default=False, action='store_true',
            help='normalize inter-contig Hi-C links by the number of links to other contigs or groups, default: %(default)s. This parameter is shared by '
            'both the clustering and reassignment steps') 
    pipe_group.add_argument(
            '--outdir', default=None,
            help='output directory, default: %(default)s (current directory)')

    # Parameters for assembly correction
    correct_group = parser.add_argument_group('>>> Parameters for assembly correction')
    correct_group.add_argument(
            '--correct_nrounds', type=int, default=0,
            help='maximum number of rounds for assembly correction, default: %(default)s. If no breakpoints are detected, the iteration will be '
            'interrupted. Set the value to 0 to skip the assembly correction step')
    correct_group.add_argument(
            '--correct_resolution', type=int, default=500,
            help='resolution for breakpoint detection and assembly correction, default: %(default)s (bp)')
    correct_group.add_argument(
            '--median_cov_ratio', type=float, default=0.2,
            help='for each contig, the coverage cutoff is calculated by multiplying the median coverage by this value, default: %(default)s. '
            'A continuous region with coverages >= the cutoff is considered as a "high-coverage region". Candidate breakpoints are determined in low-coverage regions '
            'bounded by two "large" high-coverage regions')
    correct_group.add_argument(
            '--region_len_ratio', type=float, default=0.1,
            help='a "large" high-coverage region should span a region longer than a threshold, which is the larger one of either the contig length * `region_len_ratio`, '
            'or `min_region_cutoff`, default: %(default)s')
    correct_group.add_argument(
            '--min_region_cutoff', type=int, default=5000,
            help='minimum length cutoff for a "large" high-coverage region, default: %(default)s (bp)')

    # Parameters for preprocessing (contig / link filtration) before Markov clustering
    filter_group = parser.add_argument_group('>>> Parameters for preprocessing (contig / Hi-C link filtration) before Markov clustering')
    filter_group.add_argument(
            '--Nx', type=int, default=80,
            help='keep contigs longer than Nx, default: %(default)s. This parameter speeds up the process and reduces memory consumption by reducing "n" in the '
            'Markov cluster algorithm, which has a time and space complexity of O(n^3) in its straightforward implementation. Additionally, fragmented '
            'contigs often fail to elongate due to a high proportion of repeats and assembly errors, so ignoring them in the clustering step could be better. '
            'The contigs filtered out by `--Nx` can be rescued in the following reassignment steps. To disable it, set this parameter to 100')
    filter_group.add_argument(
            '--RE_site_cutoff', type=int, default=5,
            help='minimum restriction enzyme sites, default: %(default)s. Contigs with fewer RE sites than this value will be filtered out')
    filter_group.add_argument(
            '--density_lower', default='0.2X',
            help='lower limit for contig Hi-C link density, used to filter out contigs that provide minimal linking information, default: %(default)s. '
            'This parameter works in two modes. [Fraction mode] If `density_lower` is a numeric value within the range [0, 1) , contigs with a link density '
            'less than (1 - `density_lower`) * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `density_lower` is a string '
            'ending with "X", contigs with a link density lower than `density_lower` times the median density will be filtered out. To disable it, set this '
            'parameter to 0')
    filter_group.add_argument(
            '--density_upper', default='1.9X',
            help='upper limit for contig Hi-C link density, used to filter out potential collapsed contigs, default: %(default)s. This parameter works '
            'in two modes. [Fraction mode] If `density_upper` is a numeric value within the range (0, 1] , contigs with a link density greater than '
            '`density_upper` * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `density_upper` is a string ending with "X", '
            'contigs with a link density greater than `density_upper` times the median density will be filtered out. To disable it, set this parameter to 1')
    filter_group.add_argument(
            '--read_depth_upper', default='1.5X',
            help='upper limit for contig read depth, used to filter out potential collapsed contigs, default: %(default)s. This parameter works '
            'in two modes. [Fraction mode] If `read_depth_upper` is a numeric value within the range (0, 1], contigs with a read depth greater than '
            '`read_depth_upper` * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `read_depth_upper` is a string ending with "X", '
            'contigs with a read depth greater than Q3 + `read_depth_upper` times the IQR (interquartile range) will be filtered out. To disable it, set this '
            'parameter to 1')
    filter_group.add_argument(
            '--topN', type=int, default=10,
            help='top N nearest neighbors (contigs) used for rank-sum calculation, default: %(default)s')
    filter_group.add_argument(
            '--rank_sum_hard_cutoff', type=int, default=0,
            help='hard filtering cutoff for rank-sum value, default: %(default)s (disabled). This parameter may be helpful if you are experienced in '
            'estimating the range of rank-sum values, but it can also be dangerous. Be aware that the rank-sum values are affected by variables such as `topN` '
            'and contig length')
    filter_group.add_argument(
            '--rank_sum_upper', default='1.5X',
            help='upper limit for rank-sum value, used to filter out potential chimeric and collapsed contigs, default: %(default)s. This parameter works '
            'in two modes. [Fraction mode] If `rank_sum_upper` is a numeric value within the range (0, 1], contigs with a rank-sum value greater than '
            '`rank_sum_upper` * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `rank_sum_upper` is a string ending with "X", '
            'contigs with a rank-sum value greater than Q3 + `rank_sum_upper` times the IQR (interquartile range) will be filtered out. To disable it, set this '
            'parameter to 1')
    filter_group.add_argument(
            '--remove_allelic_links', type=int, default=0,
            help='This parameter identifies allelic contig pairs and removes the Hi-C links between them. The value should be the ploidy and must be >= 2. '
            'By default, this function is disabled. It is designed to handle haplotype-phased assemblies with high switch error rates. Allelic contig pairs are '
            'identified based on the pattern of Hi-C links')
    filter_group.add_argument(
            '--concordance_ratio_cutoff', type=float, default=0.2,
            help='concordance ratio cutoff for allelic contig pair identification, default: %(default)s. In most cases, the default cutoff works well. Increasing '
            'the cutoff will increase both the true and false positive rates (TPR and FPR). This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--nwindows', type=int, default=50,
            help='this parameter is used in the concordance ratio calculation to eliminate the interference of contig length. The range (2 times the distance) '
            'for counting the "concordant Hi-C links" is dynamically defined by dividing the length of the shorter contig from a pair by `nwindows`, default: %(default)s')
    filter_group.add_argument(
            '--remove_concentrated_links', default=False, action='store_true',
            help='remove concentrated links that are concentrated in a tiny region of the genome with a remarkably higher link density, default: %(default)s')
    filter_group.add_argument(
            '--max_read_pairs', type=int, default=200,
            help='maximum Hi-C read pairs for allelic contig pair identification, default: %(default)s. This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--min_read_pairs', type=int, default=20,
            help='minimum Hi-C read pairs for allelic contig pair identification, default: %(default)s. This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--phasing_weight', type=float, default=1.0,
            help='weight of phasing information, default: %(default)s. When this parameter is set to 1.0, all Hi-C links between haplotypes will be removed. Set this'
            'parameter to 0 to completely ignore the phasing information in the GFA files. This parameter only works if `--gfa` is set and mulitple gfa files are input')

    # Parameters for parsing ultra-long reads
    ul_group = parser.add_argument_group('>>> Parameters for parsing ultra-long reads')
    ul_group.add_argument(
            '--min_ul_mapq', type=int, default=30,
            help='MAPQ cutoff, alignments with MAPQ lower than this value will be removed, default: %(default)s')
    ul_group.add_argument(
            '--min_ul_alignment_length', type=int, default=10000,
            help='alignment length cutoff, alignments shorter than this value will be removed, default: %(default)s (bp)')
    ul_group.add_argument(
            '--max_distance_to_end', type=int, default=100,
            help='the distance to the ends of a contig should be less than this value, default: %(default)s (bp)')
    ul_group.add_argument(
            '--max_overlap_ratio', type=float, default=0.5,
            help='the length ratio of overlap between primary and supplementary alignments should be less than this value, default: %(default)s')
    ul_group.add_argument(
            '--max_gap_len', type=int, default=10000,
            help='maximum gap length between primary and supplementary alignments, default: %(default)s (bp)')
    ul_group.add_argument(
            '--min_ul_support', type=int, default=2,
            help='a linkage between two contigs should be supported by more than this number of UL reads, default: %(default)s')

    # Parameters for adjacency matrix construction and Markov clustering
    mcl_group = parser.add_argument_group('>>> Parameters for adjacency matrix construction and Markov clustering')
    mcl_group.add_argument(
            '--bin_size', type=int, default=-1,
            help="if a contig's length is greater than this value (unit: Kbp), the contig will be split into bins of this size "
            "to construct an adjacency matrix. By default, `bin_size` is the smaller of either average chromosome length divided by 30, or 2000 kbp, "
            "but it will not be shorter than 100 Kbp. You can manually designate this value or set it to 0 to disable the function")
    mcl_group.add_argument(
            '--flank', type=int, default=500,
            help="consider only the Hi-C links from the flanking regions (ends) of each contig to construct the adjacency matrix, default: %(default)s (Kbp). "
            "If a contig's length is less than or equal to 2 times `flank`, all links will be used. Set this parameter to 0 to use the whole contigs")
    mcl_group.add_argument(
            '--expansion', type=int, default=2,
            help='expansion for Markov clustering, default: %(default)s')
    mcl_group.add_argument(
            '--min_inflation', type=float, default=1.1,
            help='minimum inflation for parameter tuning, default: %(default)s')
    mcl_group.add_argument(
            '--max_inflation', type=float, default=3.0,
            help='maximum inflation for parameter tuning, default: %(default)s')
    mcl_group.add_argument(
            '--inflation_step', type=float, default=0.1,
            help='step length of inflation for parameter tuning, default: %(default)s')
    mcl_group.add_argument(
            '--max_iter', type=int, default=200,
            help='maximum number of iterations of Markov clustering for each inflation, default: %(default)s')
    mcl_group.add_argument(
            '--pruning', type=float, default=0.0001,
            help='pruning threshold for Markov clustering, default: %(default)s')
    mcl_group.add_argument(
            '--skip_clustering', default=False, action='store_true',
            help='skip Markov clustering, default: %(default)s')

    # Parameters for reassignment and rescue
    reassign_group = parser.add_argument_group('>>> Parameters for reassignment and rescue')
    reassign_group.add_argument(
            '--min_group_len', type=float, default=5,
            help='minimum group length, default: %(default)s (Mbp). Groups smaller than this value will be dismissed and the contigs inside will be rescued and reassigned')
    reassign_group.add_argument(
            '--max_ctg_len', type=float, default=10000,
            help='maximum contig length for reassignment, default: %(default)s (Kbp). Contigs longer than this value will NOT be reassigned as they are error-prone. '
            'Set this parameter to 0 to rescue only, and no contig will be reassigned')
    reassign_group.add_argument(
            '--min_RE_sites', type=int, default=25,
            help='minimum restriction enzyme sites, default: %(default)s. Contigs with fewer RE sites than this value will be moved to "ungrouped"')
    reassign_group.add_argument(
            '--min_links', type=int, default=25,
            help='minimum Hi-C links, default: %(default)s. Contigs with fewer Hi-C links than this value will be moved to "ungrouped"')
    reassign_group.add_argument(
            '--min_link_density', type=float, default=0.0001,
            help='minimum Hi-C link density (= links / RE_sites), default: %(default)s. Contigs with Hi-C link density less than this value will be moved to "ungrouped"')
    reassign_group.add_argument(
            '--min_density_ratio', type=float, default=4,
            help='minimum Hi-C link density ratio (= link density to the best group / to the average of other groups), default: %(default)s. Contigs with link density ratio '
            'less than this value will NOT be reassigned and rescued')
    reassign_group.add_argument(
            '--ambiguous_cutoff', type=float, default=0.6,
            help='contigs with Hi-C links to the second-best group / to the best group >= this value are defined as ambiguous contigs, default: %(default)s. '
            'These contigs will NOT be reassigned, but only rescued in an additional round of rescue.')
    reassign_group.add_argument(
            '--reassign_nrounds', type=int, default=5,
            help='maximum number of rounds for reassignment, default: %(default)s. If the result converges, the iterations will be interrupted')
    reassign_group.add_argument(
            '--no_additional_rescue', default=False, action='store_true',
            help='do not run the additional round of rescue, default: %(default)s')

    # Parameters for fast sorting
    fast_sort_group = parser.add_argument_group('>>> Parameters for fast sorting')
    fast_sort_group.add_argument(
            '--skip_fast_sort', default=False, action='store_true',
            help='skip fast sorting and only run ALLHiC optimization, default: %(default)s')
    fast_sort_group.add_argument(
            '--flanking_region', type=int, default=0,
            help='consider only the Hi-C links from the flanking regions (ends) of each contig or scaffold during fast sorting (unit: kbp). By default, this parameter is set to '
            '0 to use the whole contigs or scaffolds')
    fast_sort_group.add_argument(
            '--density_cal_method', choices={'multiplication', 'sum', 'geometric_mean'}, default='multiplication',
            help='method for Hi-C link density calculation during fast sorting, default: %(default)s')
    fast_sort_group.add_argument(
            '--confidence_cutoff', type=float, default=1,
            help='cutoff for confidence filtering, default: %(default)s')

    # Parameters for ALLHiC optimization
    allhic_group = parser.add_argument_group('>>> Parameters for ALLHiC optimization')
    allhic_group.add_argument(
            '--skip_allhic', default=False, action='store_true',
            help='skip the entire ALLHiC optimization step, default: %(default)s')
    allhic_group.add_argument(
            '--skipGA', default=False, action='store_true',
            help='skip the genetic algorithm optimization step in ALLHiC, default: %(default)s')
    allhic_group.add_argument(
            '--mutprob', type=float, default=0.2,
            help='mutation probability in the genetic algorithm, default: %(default)s')
    allhic_group.add_argument(
            '--ngen', type=int, default=5000,
            help='number of generations for convergence, default: %(default)s')
    allhic_group.add_argument(
            '--npop', type=int, default=100,
            help='population size, default: %(default)s')
    allhic_group.add_argument(
            '--seed', type=int, default=42,
            help='random seed, default: %(default)s')

    # Parameters for building final scaffolds (pseudomolecules)
    build_group = parser.add_argument_group('>>> Parameters for building final scaffolds (pseudomolecules)')
    build_group.add_argument(
            '--Ns', type=int, default=100,
            help='number of Ns representing gaps, default: %(default)s')
    build_group.add_argument(
            '--sort_by_input',
            default=False, action='store_true',
            help='sort output scaffolds by the order of input `*.tour` files, otherwise by group length, default: %(default)s')
    build_group.add_argument(
            '--keep_letter_case',
            default=False, action='store_true',
            help='keep the letter case of bases in the draft genome (usually lowercase represents softmaksing), default: %(default)s')
    build_group.add_argument(
            '--prefix', default='scaffolds',
            help='prefix of output FASTA and AGP files, default: %(default)s')

    # Parameters for performance
    performance_group = parser.add_argument_group('>>> Parameters for performance')
    performance_group.add_argument(
            '--threads', type=int, default=8,
            help='threads for reading BAM file, default: %(default)s')
    performance_group.add_argument(
            '--dense_matrix', default=False, action='store_true',
            help='the Markov clustering is optimized with Scipy sparse matrix and Intel Math Kernal Library (MKL) by default. Set this parameter to '
            'disable it using dense matrix, default: %(default)s')
    performance_group.add_argument(
            '--processes', type=int, default=8,
            help='processes for fast sorting and ALLHiC optimization, always less than or equal to the number '
            'of input group files, default: %(default)s. Be aware that multiprocessing will increase RAM consumption')

    # Parameters for logging
    logging_group = parser.add_argument_group('>>> Parameters for logging')
    logging_group.add_argument(
            '--verbose', default=False, action='store_true',
            help='verbose logging, default: %(default)s')

    args = parser.parse_args()

    return args


def haphic_cluster(args):

    logger.info('Step1: Execute preprocessing and Markov clustering for contigs...')

    DIR = '01.cluster'
    LOG_FILE = 'HapHiC_cluster.log'
    os.mkdir(DIR)
    os.chdir(DIR)

    HapHiC_cluster.run(args, log_file=LOG_FILE)

    # modify args for the next step
    # use the corrected assembly instead the original one
    if args.correct_nrounds:
        args.fasta = abspath('corrected_asm.fa')
        args.corrected_ctgs = abspath('corrected_ctgs.txt')
        if args.quick_view and args.gfa and len(args.gfa.split(',')) >= 2:
            args.gfa = ','.join([pathjoin(pathsplit(gfa)[0], DIR, 'corrected_' + pathsplit(gfa)[1]) for gfa in args.gfa.split(',')])
    else:
        args.corrected_ctgs = None
    os.chdir('..')

    args.links = abspath(pathjoin(DIR, 'full_links.pkl'))
    args.HT_links = abspath(pathjoin(DIR, 'HT_links.pkl'))
    args.clm = abspath(pathjoin(DIR, 'paired_links.clm'))
    args.nclusters = args.nchrs

    # quick view
    if args.quick_view:
        return None

    # get recommended inflation
    with open(pathjoin(DIR, LOG_FILE)) as f:

        for line in f:
            match = re.match(r'.+You could try inflation from ([\d.]+)', line)

            if match:
                inflation = match.groups(0)[0]

                # modify args for the next step
                args.clusters = abspath(pathjoin(
                    DIR, 'inflation_{}'.format(inflation),
                    'mcl_inflation_{}.clusters.txt'.format(inflation)))

                return None

    raise RuntimeError(
            'Pipeline Aborted: Inflation recommendation failed. It seems that '
            'some chromosomes were grouped together, or the maximum number of '
            'clusters is even less than the expected number of chromosomes. For '
            'more details, please check out the logs.')


def haphic_reassign(args):

    logger.info('Step2: Reassign and rescue contigs...')

    DIR = '02.reassign'
    LOG_FILE = 'HapHiC_reassign.log'
    os.mkdir(DIR)
    os.chdir(DIR)

    HapHiC_reassign.run(args, log_file=LOG_FILE)

    os.chdir('..')

    # modify args for the next step
    args.clm_dir = abspath(pathjoin(DIR, 'split_clms'))
    args.groups = [f for f in glob.glob(abspath(pathjoin(DIR, 'final_groups/group*.txt')))]


def haphic_sort(args):

    logger.info('Step3: Order and orient contigs within each group...')

    DIR = '03.sort'
    # LOG_FILE = 'HapHiC_sort.log'
    os.mkdir(DIR)
    os.chdir(DIR)

    # To improve memory usage of multiprocessing, the HapHiC_sort.py is called via subprocess
    # instead of import HapHiC_sort
    # HapHiC_sort.run(args, log_file=LOG_FILE)
    commands = [pathjoin(os.path.dirname(os.path.realpath(__file__)), 'HapHiC_sort.py')]
    commands.extend([args.fasta, args.HT_links, args.clm_dir])
    commands.extend(args.groups)

    if args.quick_view:
        commands.append('--quick_view')
    else:
        if args.skip_fast_sort:
            commands.append('--skip_fast_sort')
        if args.skip_allhic:
            commands.append('--skip_allhic')
        if args.skipGA:
            commands.append('--skipGA')
        commands.extend([
            '--mutprob', str(args.mutprob), '--ngen', str(args.ngen), 
            '--npop', str(args.npop), '--seed', str(args.seed)])
        commands.extend(['--processes', str(args.processes)])

    commands.extend(['--flanking_region', str(args.flanking_region)])
    commands.extend(['--density_cal_method', args.density_cal_method])
    commands.extend(['--confidence_cutoff', str(args.confidence_cutoff)])

    if args.verbose:
        commands.append('--verbose')

    subprocess.run(commands, check=True)

    os.chdir('..')

    # modify args for the next step
    args.tours = [f for f in glob.glob(abspath(pathjoin(DIR, 'final_tours', '*.tour')))]


def haphic_build(args):

    logger.info('Step4: Build final scaffolds (pseudomolecules)...')

    DIR = '04.build'
    LOG_FILE = 'HapHiC_build.log'
    os.mkdir(DIR)
    os.chdir(DIR)

    if args.alignments.endswith('.pairs') or args.alignments.endswith('.pairs.gz'):
        args.alignments = '../01.cluster/alignments.bed'

    HapHiC_build.run(args, log_file=LOG_FILE)

    os.chdir('..')


def main():

    # parse arguments
    args = parse_argument()

    start_time = time.time()
    logger.info('Pipeline started, HapHiC version: {} (update: {})'.format(__version__, __update_time__))
    logger.info('Python version: {}'.format(sys.version.replace('\n', '')))
    logger.info('Command: {}'.format(' '.join(sys.argv)))

    # check steps
    steps = {int(step) for step in args.steps.split(',')}
    if steps not in ({1}, {1, 2}, {1, 2, 3}, {1, 2, 3, 4}):
        raise Exception('Illegal steps: {}'.format(args.steps))

    # modify args to make them compatible across scripts
    args.fasta = abspath(args.fasta)
    args.raw_fasta = args.fasta
    args.alignments = abspath(args.alignments)
    if args.ul:
        args.ul = abspath(args.ul)
    if args.gfa:
        args.gfa = ','.join([abspath(gfa) for gfa in args.gfa.split(',')])

    # prepare output directory
    if args.outdir:
        try:
            os.mkdir(args.outdir)
        except FileExistsError:
            logger.warning('The directory {} already exists'.format(args.outdir))
        os.chdir(args.outdir)

    # Step1: Cluster contigs into groups
    haphic_cluster(args)

    # Step2: Reassign and rescue contigs
    if 2 in steps:
        haphic_reassign(args)

    # Step3: Order and orient contigs for each group
    if 3 in steps:
        haphic_sort(args)

    # Step4: Build final scaffolds (pseudomolecules)
    if 4 in steps:
        haphic_build(args)

    finished_time = time.time()
    logger.info('HapHiC pipeline finished in {}s'.format(finished_time-start_time))


if __name__ == '__main__':
    main()

