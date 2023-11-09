#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-08-06 16:44


import argparse
from portion import closed
from itertools import combinations
import random
import collections
import numpy as np
from decimal import Decimal

def parse_allele_table(allele_table):
    
    synteny_allele_set = set()

    with open(allele_table) as f:
        for line in f:
            cols = line.split()
            if len(cols[2:]) > 1:
                for ctg_pair in combinations(cols[2:], 2):
                    synteny_allele_set.add(tuple(sorted(ctg_pair)))

    return synteny_allele_set


def parse_ctg_name(ctg):

    elements = ctg.split('_')
    homo_group = elements[0]
    source_chr = '_'.join(elements[:2])
    ctg_pos_range = closed(int(elements[3]), int(elements[4]))
    ctg_len = int(elements[6])

    return homo_group, source_chr, ctg_pos_range, ctg_len


def parse_result(result, synteny_allele_set, overlap_len_cutoff, link_cutoff, tag):
    
    # use a dict to cache ctg info
    ctg_info_dict = dict()

    # statistics for ROC
    total_ctg_pair_set = set()
    allelic_ctg_pair_set = set()
    allhic_positive_ctg_pair_set = set()
    concordance_ratio_list = list()

    with open(result) as fin, open('concordance_ratio.txt', 'w') as fout:
        
        for line in fin:
            if 'concordance_ratio' not in line:
                continue
            
            cols = line.split()
            ctg_1, ctg_2 = cols[3:5]
            links = cols[5].split('=')[-1]
            concordance_ratio = cols[6].split('=')[-1]
            sorted_ctg_pair = tuple(sorted(cols[3:5]))

            if sorted_ctg_pair in synteny_allele_set:
                is_synteny_allele = True
            else:
                is_synteny_allele = False

            if ctg_1 in ctg_info_dict:
                homo_group_1, source_chr_1, ctg_pos_range_1, ctg_len_1 = ctg_info_dict[ctg_1]
            else:
                ctg_info_dict[ctg_1] = homo_group_1, source_chr_1, ctg_pos_range_1, ctg_len_1 = parse_ctg_name(ctg_1)
            
            if ctg_2 in ctg_info_dict:
                homo_group_2, source_chr_2, ctg_pos_range_2, ctg_len_2 = ctg_info_dict[ctg_2]
            else:
                ctg_info_dict[ctg_2] = homo_group_2, source_chr_2, ctg_pos_range_2, ctg_len_2 = parse_ctg_name(ctg_2)

            ovl_len, ovl_percentage = 0, 0
            if homo_group_1 == homo_group_2 and source_chr_1 != source_chr_2:
                type_ = 'Inter_homo'
                ovl_range = ctg_pos_range_1 & ctg_pos_range_2
                if not ovl_range.empty:
                    ovl_len = ovl_range.upper - ovl_range.lower + 1
                    ovl_percentage = ovl_len * 2 / (ctg_len_1 + ctg_len_2) * 100
            elif homo_group_1 == homo_group_2:
                type_ = 'Intra_chrom'
            
            if homo_group_1 == homo_group_2:
            
                if ovl_percentage == 0:
                    ovl_percentage_interval = '0'
                elif ovl_percentage <= 10:
                    ovl_percentage_interval = '(0, 10]'
                elif ovl_percentage <= 20:
                    ovl_percentage_interval = '(10, 20]'
                elif ovl_percentage <= 30:
                    ovl_percentage_interval = '(20, 30]'
                elif ovl_percentage <= 40:
                    ovl_percentage_interval = '(30, 40]'
                elif ovl_percentage <= 50:
                    ovl_percentage_interval = '(40, 50]'
                elif ovl_percentage <= 60:
                    ovl_percentage_interval = '(50, 60]'
                elif ovl_percentage <= 70:
                    ovl_percentage_interval = '(60, 70]'
                elif ovl_percentage <= 80:
                    ovl_percentage_interval = '(70, 80]'
                elif ovl_percentage <= 90:
                    ovl_percentage_interval = '(80, 90]'
                elif ovl_percentage <= 100:
                    ovl_percentage_interval = '(90, 100]'
                else:
                    assert False
                
                # statistics for ROC
                # firstly, links should >= link_cutoff
                if int(links) >= link_cutoff:
                    # positive: overlap length >= overlap_len_cutoff
                    if ovl_len >= overlap_len_cutoff:
                        total_ctg_pair_set.add(sorted_ctg_pair)
                        allelic_ctg_pair_set.add(sorted_ctg_pair)
                        if is_synteny_allele:
                            allhic_positive_ctg_pair_set.add(sorted_ctg_pair)
                        concordance_ratio_list.append((sorted_ctg_pair, float(concordance_ratio)))
                    # negative: intra-chromosome contig pairs
                    elif type_ == 'Intra_chrom':
                        total_ctg_pair_set.add(sorted_ctg_pair)
                        if is_synteny_allele:
                            allhic_positive_ctg_pair_set.add(sorted_ctg_pair)
                        concordance_ratio_list.append((sorted_ctg_pair, float(concordance_ratio)))

                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    tag, ctg_1, ctg_2, type_, ovl_percentage, ovl_percentage_interval, 
                    concordance_ratio, links, is_synteny_allele))

    return total_ctg_pair_set, allelic_ctg_pair_set, concordance_ratio_list, allhic_positive_ctg_pair_set


def output_ROC(total_ctg_pair_set, allelic_ctg_pair_set, concordance_ratio_list, allhic_positive_ctg_pair_set, N50):
    
    # sort concordance_ratio_list
    # concordance ratio is dispersed (because number of coord pairs are limited)
    # shuffle the list before sort could be more robust
    random.seed(12345)
    random.shuffle(concordance_ratio_list)
    concordance_ratio_list.sort(key=lambda x: x[1], reverse=True)
    total_nctg_pairs = len(concordance_ratio_list)

    TP_dict = collections.defaultdict(set)
    FP_dict = collections.defaultdict(set)

    for param in np.arange(Decimal('0'), Decimal('1.05'), Decimal('0.05')):

        for ctg_pair, concordance_ratio in concordance_ratio_list:

            if concordance_ratio <= param:
                break

            if ctg_pair in allelic_ctg_pair_set:
                TP_dict[param].add(ctg_pair)
            else:
                FP_dict[param].add(ctg_pair)

    with open('ROC.txt', 'w') as fout:
        
        # for HapHiC
        fout.write('HapHiC\t{}\t{}\t{}\t{}\n'.format(N50, 1, 0, 0))
        for param, TP_ctg_set in TP_dict.items():
            FP_ctg_set = FP_dict[param]
            FN_ctg_set = allelic_ctg_pair_set - TP_ctg_set
            TN_ctg_set = (total_ctg_pair_set - allelic_ctg_pair_set) - (TP_ctg_set | FP_ctg_set)

            TPR = len(TP_ctg_set)/(len(TP_ctg_set)+len(FN_ctg_set))
            FPR = len(FP_ctg_set)/(len(FP_ctg_set)+len(TN_ctg_set))

            fout.write('HapHiC\t{}\t{}\t{}\t{}\n'.format(N50, param, TPR, FPR))
        
        # for ALLHiC
        allhic_TP_ctg_set = set()
        allhic_FP_ctg_set = set()

        for ctg_pair in allhic_positive_ctg_pair_set:
            if ctg_pair in allelic_ctg_pair_set:
                allhic_TP_ctg_set.add(ctg_pair)
            else:
                allhic_FP_ctg_set.add(ctg_pair)

        allhic_FN_ctg_set = allelic_ctg_pair_set - allhic_TP_ctg_set
        allhic_TN_ctg_set = (total_ctg_pair_set - allelic_ctg_pair_set) - (allhic_TP_ctg_set | allhic_FP_ctg_set)

        allhic_TPR = len(allhic_TP_ctg_set)/(len(allhic_TP_ctg_set)+len(allhic_FN_ctg_set))
        allhic_FPR = len(allhic_FP_ctg_set)/(len(allhic_FP_ctg_set)+len(allhic_TN_ctg_set))

        fout.write('ALLHiC\t{}\t1\t{}\t{}\n'.format(N50, allhic_TPR, allhic_FPR))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'result', help='input standard error of HapHiC_cluster.py --verbose')
    parser.add_argument(
            'allele_table', help='input "Allele.ctg.table" generated by ALLHiC')
    parser.add_argument(
            'tag', help='tag of input result (e.g. N50)')
    parser.add_argument(
            '--overlap_len_cutoff', default=10000, type=int,
            help='contig pairs with a overlap length >= this value will be considered as "allelic" in the ROC statistics, default: %(default)s')
    parser.add_argument(
            '--link_cutoff', default=20, type=int,
            help='contig pairs with Hi-C links less than this value will be filtered out in the ROC statistics, default: %(default)s')

    args = parser.parse_args()

    synteny_allele_set = parse_allele_table(args.allele_table)

    total_ctg_pair_set, allelic_ctg_pair_set, concordance_ratio_list, allhic_positive_ctg_pair_set = parse_result(
            args.result, synteny_allele_set, args.overlap_len_cutoff, args.link_cutoff, args.tag)

    output_ROC(total_ctg_pair_set, allelic_ctg_pair_set, concordance_ratio_list, allhic_positive_ctg_pair_set, args.tag)


if __name__ == '__main__':
    main()

