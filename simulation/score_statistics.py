#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-12-06 17:15

import argparse
import sys


def parse_tour(tour, chrom, N50, program, each_iteration):

    with open(tour) as f:
        start_score = ''
        ngen = 0
        for line in f:
            if line.startswith('>GA') and not line.startswith('>GA2-0'):
                score = line.strip().split('-')[-1]
                if not start_score:
                    start_score = score
                if each_iteration:
                    print('{}\t{}\t{}\t{}\t{}'.format(program, chrom, N50, ngen, score), file=sys.stderr)
                ngen += 500
    
    print('{}\t{}\t{}\t{}'.format(program, chrom, N50, score))
    if program == 'HapHiC':
        print('{}_presort\t{}\t{}\t{}'.format(program, chrom, N50, start_score))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tour', help='input tour file')
    parser.add_argument('chrom', help='chromosome')
    parser.add_argument('N50', help='contig N50')
    parser.add_argument('program', help='program, HapHiC or ALLHiC')
    parser.add_argument('--each_iteration', default=False, action='store_true',
                        help='output the score S of each iteration (STDERR), otherwise the first '
                        '(only if the program is HapHiC) and the last iteration (STDOUT), default: %(default)s')
    args = parser.parse_args()

    parse_tour(args.tour, args.chrom, args.N50, args.program, args.each_iteration)

if __name__ == '__main__':
    main()
