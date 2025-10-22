#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-06-06 10:01


import sys
import time
import argparse
import logging

from collections import defaultdict
from portion import closed
from HapHiC_cluster import parse_fasta
from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_agp(agp, min_ctg_len):

    logger.info('Parsing input AGP file...')

    ctg_group_dict = defaultdict(list)
    group_ctg_dict = defaultdict(list)
    group_len_dict = defaultdict(int)
    group_agp_lines = defaultdict(list)

    with open(agp) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            group, group_start, group_end = cols[0], int(cols[1]), int(cols[2])
            if group_end > group_len_dict[group]:
                group_len_dict[group] = group_end
            group_agp_lines[group].append(line)
            if cols[4] != 'W':
                continue
            ctg, ctg_start, ctg_end = cols[5], int(cols[6]), int(cols[7])
            if cols[8] == '+':
                ctg_orient = 1
            else:
                ctg_orient = -1
            ctg_group_dict[ctg].append((group, ctg_start, ctg_end, group_start, group_end, ctg_orient))
            group_ctg_dict[group].append((ctg, ctg_end - ctg_start + 1))

    one_ctg_groups = set()
    for group, ctg_len_list in group_ctg_dict.items():
        if len(ctg_len_list) == 1 and sum([l for _, l in ctg_len_list]) < min_ctg_len * 1000000:
            one_ctg_groups.add(group)
            for i, ctg_info in enumerate(ctg_group_dict[ctg_len_list[0][0]][::-1]):
                if ctg_info[0] == group:
                    ctg_group_dict[ctg_len_list[0][0]].pop(i)

    return ctg_group_dict, group_agp_lines, group_len_dict, one_ctg_groups


def get_max_ovl_group(groups, ctg_aln_start, ctg_aln_end):

    max_ovl_group, max_ovl = None, -1

    for group, ctg_start, ctg_end, _, __, ___ in groups:
        ovl = closed(ctg_start, ctg_end) & closed(ctg_aln_start, ctg_aln_end)
        ovl_len = 0 if ovl.empty else ovl.upper - ovl.lower + 1
        if ovl_len > max_ovl:
            max_ovl = ovl_len
            max_ovl_group = group

    return max_ovl_group


def parse_paf(paf, ctg_group_dict, aln_len_cutoff):

    logger.info('Parsing input PAF file...')

    group_ref_dict = defaultdict(dict)
    with open(paf) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if int(cols[11]) < 1:
                continue
            ctg, ctg_aln_start, ctg_aln_end, ref = cols[0], int(cols[2]), int(cols[3]), cols[5]
            if cols[4] == '+':
                orient = 1
            else:
                orient = -1
            if ctg_aln_end - ctg_aln_start < aln_len_cutoff:
                continue

            if ctg in ctg_group_dict:

                # solo contig
                if len(ctg_group_dict[ctg]) == 0:
                    continue

                ref_aln_start, ref_aln_end = int(cols[7]), int(cols[8])

                if len(ctg_group_dict[ctg]) == 1:
                    group = ctg_group_dict[ctg][0][0]
                # make it compatible with scaffolds.raw.agp after assembly correction
                else:
                    group = get_max_ovl_group(ctg_group_dict[ctg], ctg_aln_start, ctg_aln_end)

                if ref not in group_ref_dict[group]:
                    # [[(ctg, aln_len, aln_mid, ref_mid, orient)], aln_len_sum]
                    group_ref_dict[group][ref] = [
                            [(
                                ctg,
                                ctg_aln_end - ctg_aln_start + 1,
                                (ctg_aln_end - ctg_aln_start) / 2 + ctg_aln_start, 
                                (ref_aln_end - ref_aln_start) / 2 + ref_aln_start, 
                                orient)], 
                            ctg_aln_end - ctg_aln_start + 1]
                else:
                    group_ref_dict[group][ref][0].append((
                        ctg,
                        ctg_aln_end - ctg_aln_start + 1,
                        (ctg_aln_end - ctg_aln_start) / 2 + ctg_aln_start, 
                        (ref_aln_end - ref_aln_start) / 2 + ref_aln_start, 
                        orient))
                    group_ref_dict[group][ref][1] += ctg_aln_end - ctg_aln_start + 1

    return group_ref_dict


def alignment_check(group_len_dict, group_ref_dict, one_ctg_groups, aln_len_cutoff):
    
    group_list = []

    for group in group_len_dict:
        if group not in group_ref_dict and group not in one_ctg_groups:
            group_list.append(group)

    if group_list:
        err_message = 'Alignment check failed. Cannot find any alignment >= {} bp in the following group(s): {}'.format(aln_len_cutoff, ','.join(group_list))
        logger.error(err_message)
        raise Exception(err_message)


def order_and_orient_groups(ctg_group_dict, group_ref_dict, group_agp_lines, group_len_dict, one_ctg_groups, args, fa_dict=None, fout=None):

    logger.info('Ordering and orienting scaffolds based on alignments...')

    base_plus = 'ATCGNatcgn'
    base_minus = 'TAGCNtagcn'
    com_tab = str.maketrans(base_plus, base_minus)

    def revcom(seq):
        return seq.translate(com_tab)[::-1]

    def orient_seq(seq, orient):
        if orient == '+':
            return seq
        else:
            assert orient == '-'
            return revcom(seq)

    def get_reversed_orientation(orient):
        if orient == '+':
            return '-'
        assert orient == '-'
        return '+'

    def find_lis(aln_order_list, aln_len_list, forward=True):

        order_list = []
        order_aln_dict = dict()
        order_len_dict = dict()

        for i, (aln, order) in enumerate(aln_order_list):
            if forward and order < 0:
                continue
            if not forward and order > 0:
                continue
            if order in order_aln_dict:
                continue
            order_list.append(order)
            order_aln_dict[order] = aln
            order_len_dict[order] = aln_len_list[i]

        if not order_list:
            return 0, []

        dp = [0] * len(order_list)
        sequence = [None] * len(order_list)
        max_sum_idx = 0

        for i in range(len(order_list)):
            dp[i] = order_len_dict[order_list[i]]
            for j in range(i):
                if order_list[i] > order_list[j] and dp[i] < dp[j] + order_len_dict[order_list[i]]:
                    dp[i] = dp[j] + order_len_dict[order_list[i]]
                    sequence[i] = j
            if dp[i] >= dp[max_sum_idx]:
                max_sum_idx = i
        max_sum = dp[max_sum_idx]
        max_sequence = []
        while max_sum_idx is not None:
            max_sequence.append(order_aln_dict[order_list[max_sum_idx]])
            max_sum_idx = sequence[max_sum_idx]
        max_sequence.reverse()

        return max_sum, max_sequence


    ref_groups_dict = defaultdict(list)

    for group, ref_aln_dict in group_ref_dict.items():

        max_ref = max(ref_aln_dict, key=lambda x: ref_aln_dict[x][1])

        max_aln_list = ref_aln_dict[max_ref][0]

        aln_list = []
        for aln in max_aln_list:
            ctg, aln_len, aln_mid, ref_mid, orient = aln
            for _, ctg_start, ctg_end, group_start, group_end, ctg_orient in ctg_group_dict[ctg]:
                if aln_mid not in closed(ctg_start, ctg_end):
                    continue
                if orient * ctg_orient == 1:
                    order = group_start + aln_mid
                else:
                    assert orient * ctg_orient == -1
                    order = - (group_start + aln_mid)
                aln_list.append((aln, order, aln_len, ref_mid))

        aln_order_list, aln_len_list = [], []

        aln_list.sort(key=lambda x: x[-1])

        for aln, order, aln_len, _ in aln_list:
            aln_order_list.append([aln, order])
            aln_len_list.append(aln_len)

        max_sum_f, lis_aln_list_f = find_lis(aln_order_list, aln_len_list, forward=True)
        max_sum_r, lis_aln_list_r = find_lis(aln_order_list, aln_len_list, forward=False)

        logger.info('group: {}\tforward LIS: {}\treverse LIS: {}'.format(group, max_sum_f, max_sum_r))

        if max_sum_f > max_sum_r:
            ref_groups_dict[max_ref].append((group, 1, max_sum_f))
        else:
            ref_groups_dict[max_ref].append((group, -1, max_sum_r))

    if not args.ref_order:
        order_list = sorted(ref_groups_dict.keys())
    else:
        order_list = args.ref_order.split(',')

    seq_list = []
    output_groups = set()
    for ref in order_list:
        ref_groups_dict[ref].sort(key=lambda x: x[-1], reverse=True)
        for group, max_orient, _ in ref_groups_dict[ref]:
            if group in one_ctg_groups or group is None:
                continue
            output_groups.add(group)
            if fout:
                if args.keep_original_ids:
                    fout.write('>{}\n'.format(group))
                else:
                    if max_orient == 1:
                        fout.write('>{}:{}:+\n'.format(group, ref))
                    else:
                        fout.write('>{}:{}:-\n'.format(group, ref))

            if max_orient == 1:
                logger.info('{}: {} +'.format(ref, group))
                for line in group_agp_lines[group]:
                    if args.keep_original_ids:
                        print(line, end='')
                    else:
                        rest = line.split(maxsplit=1)[-1]
                        print('{}:{}:+\t{}'.format(group, ref, rest), end='')
                    if fout:
                        cols = line.split()
                        if cols[4] == 'W':
                            ctg, start, end, orient = cols[5:9]
                            seq_list.append(orient_seq(fa_dict[ctg][0][int(start)-1:int(end)], orient))
                        else:
                            Ns = int(cols[5])
                            seq_list.append('N'*Ns)

            else:
                logger.info('{}: {} -'.format(ref, group))
                for n, line in enumerate(group_agp_lines[group][::-1], 1):
                    cols = line.split()
                    group, start, end = cols[0], int(cols[1]), int(cols[2])
                    group_len = group_len_dict[group]
                    reversed_start, reversed_end = group_len - end + 1, group_len - start + 1
                    if cols[4] == 'W':
                        last_col = get_reversed_orientation(cols[-1])
                        if fout:
                            seq_list.append(orient_seq(fa_dict[cols[5]][0][int(cols[6])-1:int(cols[7])], last_col))
                    else:
                        last_col = cols[-1]
                        if fout:
                            Ns = int(cols[5])
                            seq_list.append('N'*Ns)
                    if args.keep_original_ids:
                        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            group, reversed_start, reversed_end, n, cols[4], cols[5], cols[6], cols[7], last_col))
                    else:
                        print('{}:{}:-\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            group, ref, reversed_start, reversed_end, n, cols[4], cols[5], cols[6], cols[7], last_col))
            if fout:
                seq = ''.join(seq_list)
                for i in range(0, len(seq), args.max_width):
                    fout.write('{}\n'.format(seq[i:i+args.max_width]))
                seq_list.clear()

    for group, lines in group_agp_lines.items():
        # the remaining groups (mostly unanchored contigs)
        if group not in output_groups:
            if fout:
                fout.write('>{}\n'.format(group))
            for line in lines:
                print(line, end='')
                if fout:
                    cols = line.split()
                    if cols[4] == 'W':
                        ctg, start, end, orient = cols[5:9]
                        seq_list.append(orient_seq(fa_dict[ctg][0][int(start)-1:int(end)], orient))
                    else:
                        Ns = int(cols[5])
                        seq_list.append('N'*Ns)
            if fout:
                seq = ''.join(seq_list)
                for i in range(0, len(seq), args.max_width):
                    fout.write('{}\n'.format(seq[i:i+args.max_width]))
                seq_list.clear()


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'agp', help='final scaffolds (`scaffolds.raw.agp`) or unsorted contigs in AGP format')
    parser.add_argument(
            'paf', help='alignments between the draft assembly (query) and the reference genome (target) in PAF format generated by minimap2')
    parser.add_argument(
            '--min_ctg_len', type=int, default=10,
            help='minimum length for solo contigs to be ordered and oriented, default: %(default)s (Mbp). Scaffolds with multiple contigs are always processed')
    parser.add_argument(
            '--aln_len_cutoff', type=int, default=5000,
            help='alignment length cutoff, default: %(default)s. Alignments shorter than this value are excluded from refsort')
    parser.add_argument(
            '--skip_aln_check', default=False, action='store_true',
            help='skip scaffold alignment check, default: %(default)s. By default, the program will terminate if any scaffold, or any contig longer than '
                 '`--min_ctg_len` lacks alignments meeting the length cutoff (`--aln_len_cutoff`), which may indicate inappropriate reference genome '
                 'or alignment issues. Use this flag to bypass the check')
    parser.add_argument(
            '--keep_original_ids', default=False, action='store_true',
            help='keep original sequence IDs in output files without appending reference information, default: %(default)s. By default, sequence IDs are '
                 'formatted as "original_id:ref_chrom:orientation" (e.g., group1:chr1:+, group4:chr2:-) to indicate mapping relationships. When this flag is set, '
                 'original IDs (e.g., group1, group4) are preserved instead')
    parser.add_argument(
            '--fasta', default=None,
            help='raw (uncorrected) draft genome in FASTA format. When provided, the program will output sorted scaffolds in FASTA format')
    parser.add_argument(
            '--max_width', type=int, default=60,
            help='maximum number of bases per line in the output FASTA file, default: %(default)s (bp)')
    parser.add_argument(
            '--fout', default='scaffolds.refsort.fa',
            help='file name for output FASTA file, default: %(default)s. This parameter takes effect only when `--fasta` is provided')
    parser.add_argument(
            '--ref_order', default=None,
            help='order of the reference chromosomes for outputting scaffolds, default: %(default)s (sorted by chromosome ID of the reference genome). '
                 'You can specify the order by listing chromosome IDs separated by commas, e.g., "chr1,chr2,chr3"')
    return parser.parse_args()


def run(args, log_file=None):

    if log_file:
        file_handler = logging.FileHandler(log_file, 'w')
        formatter=logging.Formatter(
                fmt='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = time.time()
    logger.info('Program started, HapHiC version: {} (update: {})'.format(__version__, __update_time__))
    logger.info('Python version: {}'.format(sys.version.replace('\n', '')))
    logger.info('Command: {}'.format(' '.join(sys.argv)))

    ctg_group_dict, group_agp_lines, group_len_dict, one_ctg_groups = parse_agp(args.agp, args.min_ctg_len)
    
    group_ref_dict = parse_paf(args.paf, ctg_group_dict, args.aln_len_cutoff)
    if not args.skip_aln_check:
        alignment_check(group_len_dict, group_ref_dict, one_ctg_groups, args.aln_len_cutoff)
    
    if args.fasta:
        fa_dict = parse_fasta(args.fasta)
        with open(args.fout, 'w') as fout:
            order_and_orient_groups(ctg_group_dict, group_ref_dict, group_agp_lines, group_len_dict, one_ctg_groups, args, fa_dict, fout)
    else:
        order_and_orient_groups(ctg_group_dict, group_ref_dict, group_agp_lines, group_len_dict, one_ctg_groups, args)

    end_time = time.time()
    logger.info('Program finished in {}s'.format(end_time-start_time))


def main():

    # get arguments
    args = parse_arguments()

    run(args, log_file='HapHiC_refsort.log')


if __name__ == '__main__':
    main()

