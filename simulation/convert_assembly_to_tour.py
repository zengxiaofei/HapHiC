#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-09 17:02

import argparse


def parse_assembly(assembly):

    ctg_dict = dict()
    output_ordering = list()

    with open(assembly) as f:
        for line in f:
            cols = line.split()
            if line.startswith('>'):
                ctg_dict[cols[1]] = cols[0][1:]
            else:
                for ID in cols:
                    if ID.startswith('-'):
                        output_ordering.append(ctg_dict[ID[1:]]+'-')
                    else:
                        output_ordering.append(ctg_dict[ID]+'+')
                
    return output_ordering


def output_tour(output_ordering, prefix):
    
    with open('{}.tour'.format(prefix), 'w') as f:
        f.write('>INIT\n')
        f.write('{}\n'.format(' '.join(output_ordering)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('assembly', help='"genome.0.assembly" file output by 3D-DNA')
    parser.add_argument('prefix', help='prefix for output tour file')
    args = parser.parse_args()

    # although 3D-DNA has a parameter --input input_size, the contigs in all the tests for ordering 
    # and orientation in this research won't exeed this size limit , so there is no need to parse fasta file
    output_ordering = parse_assembly(args.assembly)
    output_tour(output_ordering, args.prefix)

if __name__ == '__main__':
    main()
