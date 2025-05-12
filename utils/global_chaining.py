#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2025-04-01 16:26


import argparse
import sys
import collections
import re
import networkx as nx
from itertools import combinations


def parse_paf(paf, mapq, min_len, min_aln_len, div):

    """Parse PAF file, extract necessary information"""

    qry_aln_dict = collections.defaultdict(dict)
    qry_len_dict, ref_len_dict = dict(), dict()
    div_pattern = re.compile(r'.+{}:f:([0-9.]+)'.format(div))
    
    with open(paf) as f:
        for n, line in enumerate(f):
            # skip empty lines
            if not line.strip():
                continue
            
            cols = line.split()
            # MAPQ filtering
            if int(cols[11]) < mapq:
                continue
            
            # extract useful information from each line
            qry, qry_len = cols[0], int(cols[1])
            ref, ref_len = cols[5], int(cols[6])
            if min(qry_len, ref_len) < min_len:
                continue
            qry_start, qry_end = int(cols[2]), int(cols[3])
            ref_start, ref_end = int(cols[7]), int(cols[8])
            orient = 1 if cols[4] == '+' else -1
            residue_matches, block_len = int(cols[9]), int(cols[10]) 
            
            if ref_end - ref_start < min_aln_len:
                continue
            
            div = div_pattern.match(line).group(1) if div_pattern.match(line) else None
            if div is None:
                continue

            # record query and reference lengths
            if qry not in qry_len_dict:
                qry_len_dict[qry] = qry_len
            if ref not in ref_len_dict:
                ref_len_dict[ref] = ref_len
            
            # line number, alignment length, orientation * query mid, ref mid, identity, div
            info_tuple = (
                    n,
                    ref_end - ref_start + 1,
                    orient * ((qry_end - qry_start) / 2 + qry_start),
                    (ref_end - ref_start) / 2 + ref_start,
                    residue_matches,
                    block_len,
                    float(div)
                    )
            
            # record alignment information for each query-reference pair
            if ref in qry_aln_dict[qry]:
                qry_aln_dict[qry][ref].append(info_tuple)
            else:
                qry_aln_dict[qry][ref] = [info_tuple]

    return qry_aln_dict, qry_len_dict, ref_len_dict


def find_lis(aln_list, forward=True):

    """Find longest increasing subsequences, using alignment lengths as weights"""

    order_list = []
    order_aln_dict = dict()
    order_len_dict = dict()

    for i, aln in enumerate(aln_list):
        aln_len, order = aln[1:3]
        if forward and order < 0:
            continue
        if not forward and order > 0:
            continue
        if order in order_aln_dict:
            if order_len_dict[order] < aln_len:
                order_aln_dict[order] = aln
                order_len_dict[order] = aln_len
                order_list.remove(order)
                order_list.append(order)
            else:
                continue
        else:
             order_list.append(order)
             order_aln_dict[order] = aln
             order_len_dict[order] = aln_len

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


def chain_alignments(qry_aln_dict, qry_len_dict, ref_len_dict, args):
   
    print('Query\tQuery_len\tReference\tReference_len\tOrientation\tAln_len\tAln_num\tPercent_identity\tGap_compressed_Percent_identity', file=sys.stderr)
    
    G = nx.Graph()
    chained_lines = set()
    qry_ref_chained_lines = collections.defaultdict(set)
    for qry, ref_aln_list in qry_aln_dict.items():
        all_lis = []
        for ref, aln_list in ref_aln_list.items():
            # cov ratio filtering
            aln_len = sum([aln[1] for aln in aln_list])
            if aln_len / min(qry_len_dict[qry], ref_len_dict[ref]) < args.min_cov_ratio:
                continue
            # sort by ref mid
            aln_list.sort(key=lambda x: x[3])
            # find longest increasing subsequences, using alignment lengths as weights
            max_sum_f, lis_aln_list_f = find_lis(aln_list, forward=True)
            max_sum_r, lis_aln_list_r = find_lis(aln_list, forward=False)
            if max_sum_f >= max_sum_r:
                max_sum = max_sum_f
                lis_aln_list = lis_aln_list_f
                lis_orient = '+'
            else:
                max_sum = max_sum_r
                lis_aln_list = lis_aln_list_r
                lis_orient = '-'
            # record LISs for each query
            lis_info = (max_sum, lis_aln_list, ref, lis_orient)
            if not all_lis:
                all_lis.append(lis_info)
            elif max_sum > all_lis[0][0]:
                all_lis.insert(0, lis_info)
            else:
                all_lis.append(lis_info)
        
        if not all_lis:
            continue
        
        filtered_lis = [all_lis[0]]
        if len(all_lis) > 1:
            filtered_lis.extend([lis_info for lis_info in all_lis[1:] if lis_info[0] >= args.min_sb_ratio * all_lis[0][0]])

        for lis_info in filtered_lis:
            aln_len, aln_num, ref, orient = lis_info[0], len(lis_info[1]), lis_info[2], lis_info[3]
            # cov ratio filtering
            if aln_len / min(qry_len_dict[qry], ref_len_dict[ref]) < args.min_cov_ratio:
                continue
            residue_matches, block_len, div_sum = 0, 0, 0
            for aln in lis_info[1]:
                chained_lines.add(aln[0])
                qry_ref_chained_lines[frozenset({qry, ref})].add(aln[0])
                residue_matches += aln[4]
                block_len += aln[5]
                div_sum += aln[4] * aln[6]
            
            # identity filtering
            gap_compressed_identity = (1-div_sum/residue_matches)*100
            if gap_compressed_identity < args.min_identity:
                continue
            identity = residue_matches/block_len*100
            
            G.add_edge(qry, ref)
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    qry, qry_len_dict[qry], ref, ref_len_dict[ref], 
                    orient, aln_len, aln_num, identity, gap_compressed_identity 
                ), file=sys.stderr)

    return chained_lines, qry_ref_chained_lines, G


def filter_paf(paf, chained_lines, prefix='all'):

    with open(paf) as f, open('{}_chained.paf'.format(prefix), 'w') as fchain:
        for n, line in enumerate(f):
            if n in chained_lines:
                fchain.write(line)


def perform_clustering(paf, qry_ref_chained_lines, G):
    
    for n, subgraph in enumerate(nx.connected_components(G), 1):
        chained_lines = set()
        for ctg1, ctg2 in combinations(subgraph, 2):
            chained_lines |= qry_ref_chained_lines[frozenset({ctg1, ctg2})]
        filter_paf(paf, chained_lines, prefix='cluster{}'.format(n)) 


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
            'paf',
            help='aligments in PAF format (output by minimap2)')
    parser.add_argument(
            '--mapq', type=int, default=0,
            help='MAPQ filtering threshold, default: %(default)s. Alignments with a MAPQ value '
            'less than this value will be filtered out')
    parser.add_argument(
            '--min_len', type=int, default=100000,
            help='length threshold for query and reference sequences, default: %(default)s')
    parser.add_argument(
            '--min_aln_len', type=int, default=10000,
            help='alignment length threshold, default: %(default)s. Alignments with a length '
            'less than this value will be filtered out')
    parser.add_argument(
            '--div', choices={'de', 'dv'}, default='de',
            help='tag for sequence divergence, default: %(default)s')
    parser.add_argument(
            '--min_identity', type=float, default=90,
            help='minimum gap-compressed identity for the output chained alignments, default: %(default)s')
    parser.add_argument(
            '--min_cov_ratio', type=float, default=0,
            help='minimum ratio of (chained alignment length) / min(query length, reference length) to retain a query-reference combination, default: %(default)s')
    parser.add_argument(
            '--min_sb_ratio', type=float, default=0.2,
            help='minimum ratio of (secondary chained alignment length) / (best chained alignment length) to retain a query-reference combination, default: %(default)s')
    parser.add_argument(
            '--perform_clustering', action='store_true', default=False,
            help='cluster query-reference combinations, default: %(default)s')
    args = parser.parse_args()
    
    qry_aln_dict, qry_len_dict, ref_len_dict = parse_paf(args.paf, args.mapq, args.min_len, args.min_aln_len, args.div)
    chained_lines, qry_ref_chained_lines, G = chain_alignments(qry_aln_dict, qry_len_dict, ref_len_dict, args)
    
    # for all chained alignments
    filter_paf(args.paf, chained_lines)
    
    # cluster query-reference combinations
    if args.perform_clustering:
        perform_clustering(args.paf, qry_ref_chained_lines, G)


if __name__ == '__main__':
    main()

