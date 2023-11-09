#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-09-18 16:16


import argparse


def parse_allele_info(allele_info):

    allele_chrom_pos_list = list()
    allele_hap_list = list()

    with open(allele_info) as f:
        f.readline()
        for line in f:
            cols = line.split()
            chrom, pos = cols[1:3]
            hap1, hap2, hap3, hap4 = cols[4:]
            allele_chrom_pos_list.append((chrom, pos))
            allele_hap_list.append((hap1, hap2, hap3, hap4))

    return allele_chrom_pos_list, allele_hap_list


def interpret(allele_chrom_pos_list, allele_hap_list, new_chrom_pos_list, new_hap_list):

    assert len(allele_chrom_pos_list) == len(new_chrom_pos_list)

    for n, chrom_pos in enumerate(allele_chrom_pos_list):
        assert chrom_pos == new_chrom_pos_list[n]
        if allele_hap_list[n] != new_hap_list[n]:
            print('<->'.join(['hap{}'.format(m+1) for m in range(4) if allele_hap_list[n][m] != new_hap_list[n][m]]))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('allele_info', help='allele_info.txt')
    parser.add_argument('new_allele_info', help='new_allele_info.txt')
    args = parser.parse_args()

    allele_chrom_pos_list, allele_hap_list = parse_allele_info(args.allele_info)
    new_chrom_pos_list, new_hap_list = parse_allele_info(args.new_allele_info)

    interpret(allele_chrom_pos_list, allele_hap_list, new_chrom_pos_list, new_hap_list)


if __name__ == '__main__':
    main()
