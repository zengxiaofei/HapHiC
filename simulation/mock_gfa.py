#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-04-02 15:48


import argparse
import collections
import re

def is_chr(ID, chr_pattern):

    if re.match(chr_pattern, ID):
        return True
    else:
        return False

def parse_fasta(fasta, chr_pattern):
    
    hap_ctg_dict = collections.defaultdict(list)
    ctg_len_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                if is_chr(ID, chr_pattern):
                    hap = ID.split('_')[1]
                    hap_ctg_dict[hap].append(ID)
                else:
                    hap_ctg_dict['other'].append(ID)
                ctg_len_dict[ID] = 0
            else:
                ctg_len_dict[ID] += len(line.strip())

    return hap_ctg_dict, ctg_len_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='contigs in FASTA format')
    parser.add_argument('chr_pattern', help='pattern for chromosome IDs')
    args = parser.parse_args()

    hap_ctg_dict, ctg_len_dict = parse_fasta(args.fasta, args.chr_pattern)
    for hap, ctg_list in hap_ctg_dict.items():
        with open('hap{}.gfa'.format(hap), 'w') as f:
            for ctg in ctg_list:
                f.write('S\t{}\t*\tLN:i:{}\trd:i:30\n'.format(ctg, ctg_len_dict[ctg]))


if __name__ == '__main__':
    main()

