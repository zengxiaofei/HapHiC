#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-05-17 16:34

import argparse


def parse_bed(bed):

    gene_to_chr_dict = dict()
    gene_order_dict = dict()
    with open(bed) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            gene_to_chr_dict[cols[3]] = cols[0]
            if cols[0] in gene_order_dict:
                gene_order_dict[cols[0]].append(cols[3])
            else:
                gene_order_dict[cols[0]] = [cols[3]]

    return gene_to_chr_dict, gene_order_dict


def parse_simple(anchors_simple, gene_to_chr_dict1, gene_to_chr_dict2, gene_order_dict1, gene_order_dict2, chrs1, chrs2):

    chr_list1 = chrs1.split(',')
    chr_list2 = chrs2.split(',')
    
    former_genes1, former_genes2 = list(), list()
    with open(anchors_simple) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            chr1 = gene_to_chr_dict1[cols[0]]
            chr2 = gene_to_chr_dict2[cols[2]]
            if chr_list1.index(chr1) != chr_list2.index(chr2):
                continue
            if former_genes1 and gene_to_chr_dict1[former_genes1[-1]] == chr1:
                assert max([gene_order_dict1[chr1].index(g) for g in former_genes1]) < gene_order_dict1[chr1].index(cols[0])
            else:
                if cols[-1] == '-':
                    print('g*'+line, end='')
                else:
                    print(line, end='')
                former_genes1, former_genes2 = [cols[0]], [cols[2]]
                continue

            assert former_genes2 and gene_to_chr_dict2[former_genes2[-1]] == chr2
            
            if max([gene_order_dict2[chr2].index(g) for g in former_genes2]) > gene_order_dict2[chr2].index(cols[2]):
                if cols[-1] == '-':
                    print('r*'+line, end='')
                else:
                    print('b*'+line, end='')
                continue
            if cols[-1] == '-':
                print('g*'+line, end='')
                former_genes2.append(cols[2])
                continue
            print(line, end='')
            former_genes1.append(cols[0])
            former_genes2.append(cols[2])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('anchors_simple', help='*.anchors.simple file')
    parser.add_argument('bed1', help='bed file of the first genome')
    parser.add_argument('bed2', help='bed file of the second genome')
    parser.add_argument('chrs1', help='chromosomes in the first genome, separated with commas')
    parser.add_argument('chrs2', help='chromosomes in the second genome (paired), separated with commas')

    args = parser.parse_args()

    gene_to_chr_dict1, gene_order_dict1 = parse_bed(args.bed1)
    gene_to_chr_dict2, gene_order_dict2 = parse_bed(args.bed2)

    parse_simple(args.anchors_simple, gene_to_chr_dict1, gene_to_chr_dict2, gene_order_dict1, gene_order_dict2, args.chrs1, args.chrs2)


if __name__ == '__main__':
    main()
