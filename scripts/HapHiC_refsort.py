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
from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_agp(agp):

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
            group_agp_lines[cols[0]].append(line)
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
        if len(ctg_len_list) == 1 and sum([l for _, l in ctg_len_list]) < 10000000:
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
            best_group = group

    return max_ovl_group


def parse_paf(paf, ctg_group_dict):

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
            if ctg_aln_end - ctg_aln_start < 5000:
                continue

            if ctg in ctg_group_dict:

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


def order_and_orient_groups(ctg_group_dict, group_ref_dict, group_agp_lines, group_len_dict, one_ctg_groups, ref_order):

    logger.info('Ordering and orienting scaffolds based on alignments...')

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

    if not ref_order:
        order_list = sorted(ref_groups_dict.keys())
    else:
        order_list = ref_order.split(',')

    output_groups = set()
    for ref in order_list:
        ref_groups_dict[ref].sort(key=lambda x: x[-1], reverse=True)
        for group, max_orient, _ in ref_groups_dict[ref]:
            if group in one_ctg_groups or group is None:
                continue
            output_groups.add(group)
            if max_orient == 1:
                logger.info('{}: {} +'.format(ref, group))
                for line in group_agp_lines[group]:
                    print(line, end='')
            else:
                logger.info('{}: {} -'.format(ref, group))
                for n, line in enumerate(group_agp_lines[group][::-1], 1):
                    cols = line.split()
                    group, start, end = cols[0], int(cols[1]), int(cols[2])
                    group_len = group_len_dict[group]
                    reversed_start, reversed_end = group_len - end + 1, group_len - start + 1
                    if cols[4] == 'W':
                        last_col = get_reversed_orientation(cols[-1])
                    else:
                        last_col = cols[-1]
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        group, reversed_start, reversed_end, n, cols[4], cols[5], cols[6], cols[7], last_col))

    for group, lines in group_agp_lines.items():
        if group not in output_groups:
            for line in lines:
                print(line, end='')


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'agp', help='final scaffolds (`scaffolds.raw.agp) or unsorted contigs in AGP format')
    parser.add_argument(
            'paf', help='alignments between the draft assembly (query) and the reference genome (target) in PAF format generated by minimap2')
    parser.add_argument(
            '--ref_order', default=None,
            help='the order of the reference chromosomes for outputting scaffolds, default: %(default)s (sorted by chromosome ID of the reference genome). '
                 'You can manually specify the order by listing them, separated with commas, e.g., chr1,chr2,chr3')
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

    ctg_group_dict, group_agp_lines, group_len_dict, one_ctg_groups = parse_agp(args.agp)
    group_ref_dict = parse_paf(args.paf, ctg_group_dict)
    order_and_orient_groups(ctg_group_dict, group_ref_dict, group_agp_lines, group_len_dict, one_ctg_groups, args.ref_order)


def main():

    # get arguments
    args = parse_arguments()

    run(args, log_file='HapHiC_refsort.log')


if __name__ == '__main__':
    main()

