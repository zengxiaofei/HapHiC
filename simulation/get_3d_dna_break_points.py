#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-28 15:14

import argparse
from portion import closed, empty

def parse_fasta(raw_fasta):
    
    len_dict = dict()
    raw_ID_dict = dict()
    with open(raw_fasta) as f:
        for line in f:
            if line.startswith('>'):
                raw_ID = line.split()[0][1:]
                raw_ID_dict[raw_ID] = empty()
                len_dict[raw_ID] = 0
            else:
                len_dict[raw_ID] += len(line.strip())

    return len_dict, raw_ID_dict


def parse_assembly(raw_ID_dict, assembly):

    with open(assembly) as f:
        last_ID = ''
        for line in f:
            if line.startswith('>'):
                cols = line.split()
                ID = cols[0][1:]
                # broken contigs
                if ID not in raw_ID_dict:
                    frag_len = int(cols[2])
                    raw_ID = ID.split(':::')[0]
                    if raw_ID != last_ID:
                        last_ID = raw_ID
                        start = 0
                    if ID.endswith('debris'):
                        raw_ID_dict[raw_ID] |= closed(start+1, start+frag_len+1)
                    start += frag_len


def output_statistics(len_dict, raw_ID_dict, N50):
    
    for raw_ID, debris_regions in raw_ID_dict.items():
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
        if not debris_regions.empty:
            break_points = [(r.upper-r.lower-1)//2+r.lower for r in debris_regions if 0 not in debris_regions and r.upper-1 != len_dict[raw_ID]]
            if len(break_points):
                print('3D-DNA\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, len(break_points), ' '.join([str(point) for point in break_points])))
            else:
                print('3D-DNA\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, 0, 'NA'))
        else:
            print('3D-DNA\t{}\t{}\t{}\t{}\t{}'.format(N50, raw_ID, contig_type, 0, 'NA'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_fasta', help='input raw assembly in FASTA format')
    parser.add_argument('assembly', help='input assembly file')
    parser.add_argument('N50', help='N50')
    args = parser.parse_args()
    
    len_dict, raw_ID_dict = parse_fasta(args.raw_fasta)
    parse_assembly(raw_ID_dict, args.assembly)
    output_statistics(len_dict, raw_ID_dict, args.N50)

if __name__ == '__main__':
    main()
