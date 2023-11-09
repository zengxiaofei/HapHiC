#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-10-08 16:30

import argparse


def parse_fasta(fasta):

    fa_dict = dict()

    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                fa_dict[ID] = 0
            else:
                fa_dict[ID] += len(line.strip())

    return fa_dict


def parse_clusters(clusters, fa_dict):

    group_list = list()
    
    with open(clusters) as f:
        
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            # one group should have at least two contigs inside
            if len(cols) < 2:
                continue
            group_list.append(cols)    
    
    for n, ctg_list in enumerate(group_list, 1):
        with open('group{}.txt'.format(n), 'w') as f:
            for ctg in ctg_list:
                # lachesis does not split contigs
                length = fa_dict[ctg]
                f.write('{}\tNA\t{}\n'.format(ctg, length))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('clusters', help='the output file clusters.by_name.txt of LACHESIS')
    parser.add_argument('fasta', help='the draft assembly in FASTA file (to get contig lengths)')

    args = parser.parse_args()

    fa_dict = parse_fasta(args.fasta)

    parse_clusters(args.clusters, fa_dict)


if __name__ == '__main__':
    main()
