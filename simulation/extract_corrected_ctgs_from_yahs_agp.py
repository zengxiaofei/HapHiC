#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-04-03 17:25


import argparse


def parse_fasta(fasta):
    
    fa_dict = dict()

    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                fa_dict[ID] = list()
            else:
                fa_dict[ID].append(line.strip())

    for ID, seq_list in fa_dict.items():
        fa_dict[ID] = ''.join(seq_list)

    return fa_dict


def parse_agp(agp, fa_dict):
    
    corrected_ctg_dict = dict()

    with open(agp) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if cols[4] == 'W':
                ID, start, end = cols[5], int(cols[6])-1, int(cols[7])
                corrected_ctg_dict['{}_{}_{}'.format(ID, cols[6], cols[7])] = fa_dict[cols[5]][start:end]

    return corrected_ctg_dict


def output_seqs(corrected_ctg_dict):

    for new_ID, seq in corrected_ctg_dict.items():
        print('>{}\n{}'.format(new_ID, seq))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='final AGP file from yahs outputs')
    parser.add_argument('fasta', help='raw assembly file in FASTA format')
    args = parser.parse_args()

    fa_dict = parse_fasta(args.fasta)

    corrected_ctg_dict = parse_agp(args.agp, fa_dict)

    output_seqs(corrected_ctg_dict)

if __name__ == '__main__':
    main()
