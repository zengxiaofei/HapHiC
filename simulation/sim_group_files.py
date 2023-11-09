#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-07-24 17:05


import argparse
import collections

def parse_fasta(fasta_file):

    group_ctg_dict = collections.defaultdict(list)
    fa_dict = collections.defaultdict(list)    
    
    with open(fasta_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                group = ctg.split('_')[0]
                assert ctg not in group_ctg_dict[group]
                group_ctg_dict[group].append(ctg)
            else:
                fa_dict[ctg].append(line.strip().upper())
    
    for ctg, seq_list in fa_dict.items():
        seq = ''.join(seq_list)
        # Using fixed 'GATC' here is ok. Becuase in the sorting step,
        # we don't acre about the restriction sites
        RE_sites = seq.count('GATC')
        fa_dict[ctg] = [seq, len(seq), RE_sites]

    return fa_dict, group_ctg_dict


def generate_group_files(fa_dict, group_ctg_dict):
    
    # sort by length
    for group, ctgs in group_ctg_dict.items():
        with open('group_{}.txt'.format(group), 'w') as fout:
            fout.write('#Contig\tRECounts\tLength\n')
            sorted_ctgs = sorted(ctgs, key=lambda x: fa_dict[x][1], reverse=True)
            for ctg in sorted_ctgs:
                fout.write('{}\t{}\t{}\n'.format(ctg, fa_dict[ctg][2], fa_dict[ctg][1]))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='contig-level genome file in FASTA format')
    args = parser.parse_args()

    fa_dict, group_ctg_dict = parse_fasta(args.fasta)

    generate_group_files(fa_dict, group_ctg_dict)


if __name__ == '__main__':
    main()

