#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-05-17 16:01

import argparse
import collections

def get_seq_len(genome):

    len_dict = dict()
    with open(genome) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                len_dict[ID] = 0
            else:
                len_dict[ID] += len(line.strip())

    return len_dict


def parse_bed(bed, len_dict):

    def reverse_symbol(symbol):
        assert symbol in {'+', '-', '.'}
        if symbol == '+':
            return '-'
        elif symbol == '-':
            return '+'
        else:
            return '.'

    output_dict = collections.OrderedDict()
    with open(bed) as f:
        for line in f:
            cols = line.split()
            seq = cols[0]
            start, end = int(cols[1]), int(cols[2])
            new_start = len_dict[seq] - end
            new_end = len_dict[seq] - start
            if seq in output_dict:
                output_dict[seq].append([new_start, '{}\t{}\t{}\t{}\t{}\t{}'.format(seq, new_start, new_end, cols[3], cols[4], reverse_symbol(cols[5]))])
            else:
                output_dict[seq] = [[new_start, '{}\t{}\t{}\t{}\t{}\t{}'.format(seq, new_start, new_end, cols[3], cols[4], reverse_symbol(cols[5]))]]
    
    return output_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', help='a bed file (0-based)')
    parser.add_argument('genome', help='the genome file in FASTA format')
    args = parser.parse_args()

    len_dict = get_seq_len(args.genome)
    output_dict = parse_bed(args.bed, len_dict)

    for seq, info in output_dict.items():
        for line in sorted(info, key=lambda x:x[0]):
            print(line[1])


if __name__ == '__main__':
    main()

