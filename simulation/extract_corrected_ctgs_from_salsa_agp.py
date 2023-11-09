#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-04-03 17:25


import argparse
import collections

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
    frag_dict = collections.defaultdict(list)
    with open(agp) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if cols[4] == 'W':
                if cols[5] not in fa_dict:
                    ID, n = cols[5].rsplit('_', 1)
                    frag_dict[ID].append((int(n), int(cols[6]), int(cols[7])))
                else:
                    ID, start, end = cols[5], int(cols[6])-1, int(cols[7])
                    corrected_ctg_dict['{}_{}_{}'.format(ID, cols[6], cols[7])] = fa_dict[cols[5]][start:end]
    
    for ID, frag_list in frag_dict.items():
        accumulated_pos = 0
        for n, start, end in sorted(frag_list):
            real_start, real_end = start+accumulated_pos-1, end+accumulated_pos
            corrected_ctg_dict['{}_{}_{}'.format(ID, real_start+1, real_end)] = fa_dict[ID][real_start:real_end]
            accumulated_pos += end

    return corrected_ctg_dict


def output_seqs(corrected_ctg_dict):

    for new_ID, seq in corrected_ctg_dict.items():
        print('>{}\n{}'.format(new_ID, seq))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='final AGP file from salsa outputs')
    parser.add_argument('fasta', help='raw assembly file in FASTA format')
    args = parser.parse_args()

    fa_dict = parse_fasta(args.fasta)

    corrected_ctg_dict = parse_agp(args.agp, fa_dict)

    output_seqs(corrected_ctg_dict)

if __name__ == '__main__':
    main()
