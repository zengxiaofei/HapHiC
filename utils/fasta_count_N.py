#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-05-10 09:34

import argparse
import collections
import re

def parse_fasta(fasta):
    count_dict = collections.defaultdict(int)
    with open(fasta) as f:
        for line in f:
            lstr = line.strip()
            if not line.startswith('>') and lstr:
                for ns in re.findall(r'N+', lstr):
                    count_dict[ns] += 1
    return count_dict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    args = parser.parse_args()

    count_dict = parse_fasta(args.fasta)
    print(count_dict)
    total_N = 0
    for Ns, num in count_dict.items():
        total_N += len(Ns)*num
    print('total_Ns:', total_N)

if __name__ == '__main__':
    main()
