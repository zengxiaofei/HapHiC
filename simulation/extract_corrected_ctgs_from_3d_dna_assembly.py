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


def parse_assembly(assembly, fa_dict):
    
    corrected_ctg_dict = dict()
    frag_dict = dict() 
    with open(assembly) as f:
        for line in f:
            if not line.startswith('>') or line.startswith('>hic_gap_'):
                continue
            cols = line.split()
            ID = cols[0].split(':::fragment')[0][1:]
            frag_len = int(cols[2])
            if ID in frag_dict:
                frag_dict[ID].append((accumulated_pos+1, accumulated_pos+frag_len))
                accumulated_pos += frag_len
            else:
                frag_dict[ID] = [(1, frag_len)]
                accumulated_pos = frag_len

    for ID, frag_list in frag_dict.items():
        for start, end in frag_list:
            corrected_ctg_dict['{}_{}_{}'.format(ID, start, end)] = fa_dict[ID][start-1:end]

    return corrected_ctg_dict


def output_seqs(corrected_ctg_dict):

    for new_ID, seq in corrected_ctg_dict.items():
        print('>{}\n{}'.format(new_ID, seq))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('assembly', help='final .assembly file from 3D-DNA outputs')
    parser.add_argument('fasta', help='raw assembly file in FASTA format')
    args = parser.parse_args()

    fa_dict = parse_fasta(args.fasta)

    corrected_ctg_dict = parse_assembly(args.assembly, fa_dict)

    output_seqs(corrected_ctg_dict)

if __name__ == '__main__':
    main()
