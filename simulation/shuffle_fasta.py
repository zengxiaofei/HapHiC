#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-30 14:56

import os
import argparse
import collections
import random


def parse_fasta(fasta):
    
    seq_dict = collections.OrderedDict()
    
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                seq_dict[ID] = list()
            else:
                seq_dict[ID].append(line)

    return seq_dict


def shuffle(fasta, seq_dict, seed, offset):
    
    ID_list = list(seq_dict.keys())
    random.seed(seed+offset)
    random.shuffle(ID_list)
    
    with open('shuffled_'+os.path.basename(fasta), 'w') as f:
        for ID in ID_list:
            f.write('>{}\n'.format(ID))
            for line in seq_dict[ID]:
                f.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('--seed', type=int, default=12345, help='ramdom seed, default: %(default)s')
    parser.add_argument('--offset', type=int, default=0, help='seed offset, default: %(default)s')
    args = parser.parse_args()
    
    seq_dict = parse_fasta(args.fasta)
    shuffle(args.fasta, seq_dict, args.seed, args.offset)

if __name__ == '__main__':
    main()
