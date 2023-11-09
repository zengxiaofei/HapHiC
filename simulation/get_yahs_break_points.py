#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-28 15:14

import argparse

def parse_fasta(raw_fasta):

    raw_ID_dict = dict()
    with open(raw_fasta) as f:
        for line in f:
            if line.startswith('>'):
                raw_ID = line.split()[0][1:]
                raw_ID_dict[raw_ID] = list()

    return raw_ID_dict


def parse_agp(raw_ID_dict, agp):

    with open(agp) as f:
        for line in f:
            cols = line.split()
            if cols[4] == 'W':
                ID = cols[5]
                assert ID in raw_ID_dict
                raw_ID_dict[ID].append(int(cols[7]))


def output_statistics(raw_ID_dict, N50):
    for raw_ID, break_points in raw_ID_dict.items():
        if 'chimeric' in raw_ID:
            if 'inter_homo' in raw_ID:
                contig_type = 'Inter_homo'
            elif 'inner_chrom' in raw_ID:
                contig_type = 'Intra_chrom'
            else:
                assert 'inter_nonhomo' in raw_ID
                contig_type = 'Inter_nonhomo'
        else:
            contig_type = 'Non_chimeric'
        if len(break_points) > 1:
            sorted_break_points = sorted(break_points)[:-1]
            print('YaHS\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, len(sorted_break_points), ' '.join([str(point) for point in sorted_break_points])))
        else:
            print('YaHS\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, 0, 'NA'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_fasta', help='input raw assembly in FASTA format')
    parser.add_argument('agp', help='input final agp file')
    parser.add_argument('N50', help='N50')
    args = parser.parse_args()
    
    raw_ID_dict = parse_fasta(args.raw_fasta)
    parse_agp(raw_ID_dict, args.agp)
    output_statistics(raw_ID_dict, args.N50)

if __name__ == '__main__':
    main()
