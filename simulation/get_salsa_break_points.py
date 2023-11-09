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
                # broken contigs
                if ID not in raw_ID_dict:
                    raw_ID, n = ID.rsplit('_', 1)
                    raw_ID_dict[raw_ID].append((int(n), int(cols[7])))


def output_statistics(raw_ID_dict, N50):
    for raw_ID, frag_len_list in raw_ID_dict.items():
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
        if frag_len_list:
            p = 0
            break_points = list()
            for n, frag_len in sorted(frag_len_list)[:-1]:
                break_points.append(frag_len + p)
                p += frag_len
            print('SALSA2\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, len(break_points), ' '.join([str(point) for point in break_points])))
        else:
            print('SALSA2\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, 0, 'NA'))


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
