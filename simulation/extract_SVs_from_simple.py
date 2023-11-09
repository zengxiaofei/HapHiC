#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-05-19 18:07

import argparse


def parse_bed(bed):

    gene_pos_dict = dict()
    with open(bed) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if cols[0] == 'Chr02':
                gene_pos_dict[cols[3]] = [int(cols[1]), int(cols[2])]
    return gene_pos_dict


def parse_simple(simple, bed, gene_pos_dict):
    
    prefix = bed.rsplit('.')[0]
    g, b, r = 0, 0, 0
    with open(simple) as f, open('{}_inv.bed'.format(prefix), 'w') as finv, open('{}_trans.bed'.format(prefix), 'w') as ftrans, open('{}_invtr.bed'.format(prefix), 'w') as finvtr:
        for line in f:
            if not line.strip() or "*" not in line:
                continue
            sv_type = line[0]
            cols = line[2:].split()
            if cols[0] in gene_pos_dict:
                assert cols[1] in gene_pos_dict
                pos_list = gene_pos_dict[cols[0]] + gene_pos_dict[cols[1]]
                start, end = min(pos_list), max(pos_list)
                if sv_type == 'g':
                    g += 1
                    finv.write('Chr02\t{}\t{}\tINV{}\t0\t+\n'.format(start, end, g))
                elif sv_type == 'b':
                    b += 1
                    ftrans.write('Chr02\t{}\t{}\tTRANS{}\t0\t+\n'.format(start, end, b))
                else:
                    assert sv_type == 'r'
                    r += 1
                    finvtr.write('Chr02\t{}\t{}\tINVTR{}\t0\t+\n'.format(start, end, r))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('simple', help='.simple file')
    parser.add_argument('gene_bed', help='.bed file of genes')
    args = parser.parse_args()

    gene_pos_dict = parse_bed(args.gene_bed)
    parse_simple(args.simple, args.gene_bed, gene_pos_dict)


if __name__ == '__main__':
    main()
