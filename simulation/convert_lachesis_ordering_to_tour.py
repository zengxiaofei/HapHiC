#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-09 15:53

import argparse
import os

def parse_fasta(fasta):
    
    fa_list = list()

    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                fa_list.append(line.split()[0][1:])

    return fa_list


def parse_ordering_files(ordering_files, fa_list):
    
    sorted_ordering_files = sorted(ordering_files, key=lambda x: int(os.path.basename(x).split('.')[0].replace('group', '')))
    
    output_ordering = list()

    for ordering_file in sorted_ordering_files:
        with open(ordering_file) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                cols = line.split()
                if cols[1] in fa_list:
                    if cols[2] == '0':
                        output_ordering.append(cols[1]+'+')
                    else:
                        output_ordering.append(cols[1]+'-')

    # for ID in fa_list:
    #     if ID+'+' not in output_ordering and ID+'-' not in output_ordering:
    #         output_ordering.append(ID+'+')

    return output_ordering


def output_tour(output_ordering, prefix):

    with open('{}.tour'.format(prefix), 'w') as f:
        f.write('>INIT\n')
        f.write('{}\n'.format(' '.join(output_ordering)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input contigs from a specific chromosome in FASTA format')
    parser.add_argument('prefix', help='prefix for output tour file')
    parser.add_argument('ordering_files', nargs='+', help='ordering files, the program will sort them automatically')
    args = parser.parse_args()

    fa_list = parse_fasta(args.fasta)

    output_ordering = parse_ordering_files(args.ordering_files, fa_list)

    if len(output_ordering) > 1:
        output_tour(output_ordering, args.prefix)

if __name__ == '__main__':
    main()
