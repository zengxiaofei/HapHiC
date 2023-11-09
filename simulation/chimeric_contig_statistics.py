#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-04-09 20:24

import argparse
import collections
import numpy as np
from decimal import Decimal
import re

def parse_fasta(fasta):
    
    total_ctg_set = set()
    inter_homo_ctg_set = set()
    inter_nonhomo_ctg_set = set()
    inner_chrom_ctg_set = set()

    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                total_ctg_set.add(ID)
                if 'inter_homo' in line:
                    inter_homo_ctg_set.add(ID)
                elif 'inter_nonhomo' in line:
                    inter_nonhomo_ctg_set.add(ID)
                elif 'inner' in line:
                    inner_chrom_ctg_set.add(ID)

    return total_ctg_set, inter_homo_ctg_set, inter_nonhomo_ctg_set, inner_chrom_ctg_set


def parse_result(result, total_ctg_set, N50):
    
    rank_sum_list = list()
    
    with open(result) as f, open('rank_sum.txt', 'w') as fout:
        for line in f:
            match = re.match(r'.+rank sum filtering] Fragment ([\w+-]+).+rank sum=(\d+)', line)
            if match:
                ctg, rank_sum = match.groups()
                rank_sum_list.append(ctg)
                fout.write('{}\t{}\t{}\n'.format(N50, ctg, rank_sum))
            
    assert len(total_ctg_set) == len(rank_sum_list)
    
    return rank_sum_list


def output_ROC(fout, total_ctg_set, rank_sum_list, chimeric_set, chimeric_type, N50):
    
    TP_dict = collections.defaultdict(set)
    FP_dict = collections.defaultdict(set)

    for rank_sum_upper in np.arange(Decimal('0'), Decimal('1.025'), Decimal('0.025')):
        
        param = rank_sum_upper
        rank_sum_upper = int(len(rank_sum_list) * rank_sum_upper)

        for ctg in rank_sum_list[rank_sum_upper:]:
            if ctg in chimeric_set:
                TP_dict[param].add(ctg)
            else:
                FP_dict[param].add(ctg)
   

    fout.write('HapHiC\t{}\t{}\t{}\t{}\t{}\n'.format(N50, chimeric_type, 1, 0, 0))
    
    for param, TP_ctg_set in TP_dict.items():
        FP_ctg_set = FP_dict[param]
        FN_ctg_set = chimeric_set - TP_ctg_set
        TN_ctg_set = (total_ctg_set - chimeric_set) - (TP_ctg_set | FP_ctg_set)
        
        TPR = len(TP_ctg_set)/(len(TP_ctg_set)+len(FN_ctg_set))
        FPR = len(FP_ctg_set)/(len(FP_ctg_set)+len(TN_ctg_set))

        fout.write('HapHiC\t{}\t{}\t{}\t{}\t{}\n'.format(N50, chimeric_type, param, TPR, FPR))


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file (to get all contig names)')
    parser.add_argument('result', help='input standard error of HapHiC_cluster.py')
    parser.add_argument('program', help='program used for chimeric contig detection')
    parser.add_argument('N50', help='contig N50 of input fasta')
    
    args = parser.parse_args()

    total_ctg_set, inter_homo_ctg_set, inter_nonhomo_ctg_set, inner_chrom_ctg_set = parse_fasta(args.fasta)

    rank_sum_list = parse_result(args.result, total_ctg_set, args.N50)

    with open('ROC.txt', 'w') as fout:
        output_ROC(fout, total_ctg_set, rank_sum_list, inter_homo_ctg_set, 'inter_homo', args.N50)
        output_ROC(fout, total_ctg_set, rank_sum_list, inter_nonhomo_ctg_set, 'inter_nonhomo', args.N50)
        output_ROC(fout, total_ctg_set, rank_sum_list, inner_chrom_ctg_set, 'intra_chrom', args.N50)


if __name__ == '__main__':
    main()
