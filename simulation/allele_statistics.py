#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-09-26 15:03


import argparse


def get_allele_type(base1, base2):

    if base1 == base2:
        return 'Identical'
    elif '-' in (base1, base2) or len(base1) != len(base2):
        return 'InDel'
    elif base1 in 'AG' and base2 in 'TC':
        return 'SNP_transversion'
    elif base2 in 'AG' and base1 in 'TC':
        return 'SNP_transversion'
    else:
        return 'SNP_transition'


def parse_allele_info(allele_info, prefix):

    with open(allele_info) as fin, open(prefix+'.txt', 'w') as fout:
        # skip the first line
        fin.readline()
        for line in fin:
            n, chrom, ref_coord, ref_base, hap_1, hap_2, hap_3, hap_4 = line.split()
            if chrom != 'Chr1_1':
                break
            hap12_type = get_allele_type(hap_1, hap_2)
            hap13_type = get_allele_type(hap_1, hap_3)
            hap14_type = get_allele_type(hap_1, hap_4)
            hap23_type = get_allele_type(hap_2, hap_3)
            hap24_type = get_allele_type(hap_2, hap_4)
            hap34_type = get_allele_type(hap_3, hap_4)
        
            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                line.strip(), hap12_type, hap13_type, hap14_type,
                hap23_type, hap24_type, hap34_type))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('allele_info', help='input allele_info.txt')
    parser.add_argument('prefix', help='prefix for output file')

    args = parser.parse_args()
    
    parse_allele_info(args.allele_info, args.prefix)

if __name__ == '__main__':
    main()
