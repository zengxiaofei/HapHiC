#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-09-10 15:26


import argparse
import gzip
import collections


def parse_liftover(liftover):
    
    id_dict = collections.defaultdict(list)
    with open(liftover) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            if cols[4] == 'W':
                id_dict[cols[5]].append((cols[0], int(cols[6])))

    return id_dict


def read_gfa(gfa, id_dict):

    def get_converted_id(old_ID):
        assert old_ID in id_dict
        if len(id_dict[old_ID]) > 1:
            id_dict[old_ID].sort(key=lambda x: x[1])
            return '_'.join([ctg for ctg, start in id_dict[old_ID]])
        else:
            assert len(id_dict[old_ID]) == 1
            return id_dict[old_ID][0][0]

    if gfa.endswith('.gz'):
        fopen = gzip.open
    else:
        fopen = open

    with fopen(gfa, 'rt') as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if cols[0] == 'S' or cols[0] == 'A':
                new_id = get_converted_id(cols[1])
                print('{}\t{}\t{}'.format(cols[0], new_id, '\t'.join(cols[2:])))
            else:
                assert cols[0] == 'L'
                new_id_1 = get_converted_id(cols[1])
                new_id_2 = get_converted_id(cols[3])
                print('L\t{}\t{}\t{}\t{}\t{}'.format(new_id_1, cols[2], new_id_2, cols[4], '\t'.join(cols[5:])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gfa', help='the GFA file output by assemblers, gzipped files are acceptable')
    parser.add_argument('liftover', help='the *.liftover.agp file output by the `juicer pre` of YaHS')
    args = parser.parse_args()

    id_dict = parse_liftover(args.liftover)
    read_gfa(args.gfa, id_dict)

if __name__ == '__main__':
    main()

