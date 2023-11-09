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
    
    two_haps_ctg_set = set()
    three_haps_ctg_set = set()
    four_haps_ctg_set = set()

    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                total_ctg_set.add(ID)
                haps = ID.split('_')[1]
                if len(haps) == 2:
                    two_haps_ctg_set.add(ID)
                elif len(haps) == 3:
                    three_haps_ctg_set.add(ID)
                elif len(haps) == 4:
                    four_haps_ctg_set.add(ID)
                else:
                    assert len(haps) == 1
                    assert 'collapsed' not in ID

    return total_ctg_set, two_haps_ctg_set, three_haps_ctg_set, four_haps_ctg_set


def parse_result(result, total_ctg_set, method, tag):
    
    link_density_list = list()
    
    with open(result) as f, open('link_density.txt', 'w') as fout:
        for line in f:
            if method == 'link_density':
                match = re.match(r'.+link density filtering] Fragment ([\w+-]+).+density=([\d.]+)', line)
            else:
                match = re.match(r'.+rank sum filtering] Fragment ([\w+-]+).+rank sum=(\d+)', line)
            if match:
                ctg, link_density = match.groups()
                link_density_list.append(ctg)
                fout.write('{}\t{}\t{}\n'.format(tag, ctg, link_density))
            
    assert len(total_ctg_set) == len(link_density_list)
    
    return link_density_list


def output_ROC(fout, total_ctg_set, link_density_list, chimeric_set, collapse_type, N50):
    
    TP_dict = collections.defaultdict(set)
    FP_dict = collections.defaultdict(set)

    for link_density_upper in np.arange(Decimal('0'), Decimal('1.025'), Decimal('0.025')):
        
        param = link_density_upper
        link_density_upper = int(len(link_density_list) * link_density_upper)

        for ctg in link_density_list[link_density_upper:]:
            if ctg in chimeric_set:
                TP_dict[param].add(ctg)
            else:
                FP_dict[param].add(ctg)
   

    fout.write('HapHiC\t{}\t{}\t{}\t{}\t{}\n'.format(N50, collapse_type, 1, 0, 0))
    
    for param, TP_ctg_set in TP_dict.items():
        FP_ctg_set = FP_dict[param]
        FN_ctg_set = chimeric_set - TP_ctg_set
        TN_ctg_set = (total_ctg_set - chimeric_set) - (TP_ctg_set | FP_ctg_set)
        
        TPR = len(TP_ctg_set)/(len(TP_ctg_set)+len(FN_ctg_set))
        FPR = len(FP_ctg_set)/(len(FP_ctg_set)+len(TN_ctg_set))

        fout.write('HapHiC\t{}\t{}\t{}\t{}\t{}\n'.format(N50, collapse_type, param, TPR, FPR))


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file (to get all contig names)')
    parser.add_argument('result', help='input standard error of HapHiC_cluster.py')
    parser.add_argument('program', help='program used for chimeric contig detection')
    parser.add_argument('tag', help='tag of input fasta')
    parser.add_argument('--method', default='link_density', choices={'link_density', 'rank_sum'},
            help='method used for collapsed contig identification')

    args = parser.parse_args()

    total_ctg_set, two_haps_ctg_set, three_haps_ctg_set, four_haps_ctg_set = parse_fasta(args.fasta)

    link_density_list = parse_result(args.result, total_ctg_set, args.method, args.tag)

    with open('ROC.txt', 'w') as fout:
        output_ROC(fout, total_ctg_set, link_density_list, two_haps_ctg_set, 'two_hap_collapsed', args.tag)
        output_ROC(fout, total_ctg_set, link_density_list, three_haps_ctg_set, 'three_hap_collapsed', args.tag)
        output_ROC(fout, total_ctg_set, link_density_list, four_haps_ctg_set, 'four_hap_collapsed', args.tag)


if __name__ == '__main__':
    main()
