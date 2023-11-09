#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-05-22 15:59

import argparse
from itertools import combinations

def parse_allele(allele_table):
    
    allele_ctg_set = set()
    with open(allele_table) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            ctgs = cols[2:]
            if len(ctgs) <= 1:
                continue
            for ctg_pair in combinations(ctgs, 2):
                # not from a same chromosome
                if ctg_pair[0].split('_')[:2] != ctg_pair[1].split('_')[:2]:
                    allele_ctg_set.add(tuple(sorted(ctg_pair)))
    
    return allele_ctg_set


def parse_correlation_file(cor_file, allele_ctg_set):
    
    with open(cor_file) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if tuple(sorted(cols[:2])) in allele_ctg_set:
                is_allele = 'True'
            else:
                is_allele = 'False'
            print('{}\t{}\t{}'.format(line.strip(), min(float(cols[4]), float(cols[6])), is_allele))


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cor_file', help='input correlation.txt')
    parser.add_argument('allele_table', help='input Allele.ctg.table')
    args = parser.parse_args()

    allele_ctg_set = parse_allele(args.allele_table)
    parse_correlation_file(args.cor_file, allele_ctg_set)


if __name__ == '__main__':
    main()

