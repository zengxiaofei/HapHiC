#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-10-22 11:20


import argparse
import collections

def parse_fasta(fasta):
    len_dict = collections.OrderedDict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                len_dict[ID] = 0
            else:
                len_dict[ID] += len(line.strip())
    return len_dict


def mock_agp(len_dict):
    for ID, length in len_dict.items():
        print('{0}\t1\t{1}\t1\tW\t{0}\t1\t{1}\t+'.format(ID, length))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    args = parser.parse_args()

    len_dict = parse_fasta(args.fasta)
    mock_agp(len_dict)

if __name__ == '__main__':
    main()
