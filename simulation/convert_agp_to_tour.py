#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-09 16:39

import argparse


def parse_agp(agp):

    output_ordering = list()
    with open(agp) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            if cols[4] == 'W':
                output_ordering.append(cols[5]+cols[8])
    
    return output_ordering

def output_tour(output_ordering, prefix):

    with open('{}.tour'.format(prefix), 'w') as f:
        f.write('>INIT\n')
        f.write('{}\n'.format(' '.join(output_ordering)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='agp file output by YaHs or SALSA2')
    parser.add_argument('prefix', help='prefix for output tour file')
    args = parser.parse_args()
    
    # YaHs and SALSA2 will output all contigs in the final agp file,
    # so there is no need to parse fasta file
    output_ordering = parse_agp(args.agp)

    output_tour(output_ordering, args.prefix)


if __name__ == '__main__':
    main()

