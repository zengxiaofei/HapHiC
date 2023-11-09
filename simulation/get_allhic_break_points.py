#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-28 15:14

import argparse

def parse_fastas(raw_fasta, corrected_fasta):

    raw_ID_dict = dict()
    with open(raw_fasta) as f:
        for line in f:
            if line.startswith('>'):
                raw_ID = line.split()[0][1:]
                raw_ID_dict[raw_ID] = list()

    with open(corrected_fasta) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                # broken contigs
                if ID not in raw_ID_dict:
                    splits = ID.rsplit('_', 2)
                    raw_ID = splits[0]
                    raw_ID_dict[raw_ID].append(int(splits[-1]))

    return raw_ID_dict


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
        if break_points:
            print('ALLHiC\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, len(break_points)-1, ' '.join([str(point) for point in sorted(break_points)[:-1]])))
        else:
            print('ALLHiC\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, 0, 'NA'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_fasta', help='input raw assembly in FASTA format')
    parser.add_argument('corrected_fasta', help='input corrected assembly in FASTA format')
    parser.add_argument('N50', help='N50')
    args = parser.parse_args()
    
    raw_ID_dict = parse_fastas(args.raw_fasta, args.corrected_fasta)
    output_statistics(raw_ID_dict, args.N50)

if __name__ == '__main__':
    main()
