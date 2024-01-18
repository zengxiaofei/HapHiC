#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-04-21 16:45

import argparse
import os
import re


def parse_bam(bam, mapq, NM, remove_dup, remove_singletons, single_end_mapq_filtering, threads):
    if not remove_dup:
        cmd = 'samtools view -h {} -@ {}'.format(bam, threads)
    else:
        # skip read pair when DUP (1024 -> '0b10000000000')
        cmd = 'samtools view -h {} -F 1024 -@ {}'.format(bam, threads)
    
    with os.popen(cmd) as f:
        line = f.readline()
        while line:
            if line.startswith('@'):
                print(line, end='')
            else:
                line2 = f.readline()
                if not line2:
                    break
                cols1 = line.split()
                cols2 = line2.split()
                if cols1[0] != cols2[0]:
                    if remove_singletons:
                        line = line2
                        continue
                    else:
                        raise Exception('BAM may be coord-sorted or has singletons. Sort it by read name or try --remove_singletons')
                # filter NM
                if NM:
                    match = re.match(r'.+NM:i:(\d+)', line)
                    match2 = re.match(r'.+NM:i:(\d+)', line2)
                    if int(match.groups()[0]) >= NM or int(match2.groups()[0]) >= NM:
                        line = f.readline()
                        continue
                # remove read pair if MAPQ < cutoff
                if single_end_mapq_filtering:
                    if int(cols1[4]) >= mapq or int(cols2[4]) >= mapq:
                        print(line+line2, end='')
                else:
                    if int(cols1[4]) >= mapq and int(cols2[4]) >= mapq:
                        print(line+line2, end='')
            line = f.readline()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'bam', help='input BAM file, should not be sorted or name-sorted, do NOT sort it by coordinate')
    parser.add_argument(
            'mapq', type=int, help='MAPQ cutoff, read pairs with both MAPQ >= this value will be kept')
    parser.add_argument(
            '--single_end_mapq_filtering', default=False, action='store_true', 
            help='when this parameter is added, either end of a read pair having a MAPQ >= `mapq` will be kept')
    parser.add_argument(
            '--NM', type=int, default=None,
            help='edit distance cutoff, read pairs with single-end NM >= this value will be removed, default:%(default)s')
    parser.add_argument(
            '--remove_dup', default=False, action='store_true',
            help='remove PCR duplicates, default: %(default)s. Note that the duplicates should have been marked in the BAM file (flag 1024)')
    parser.add_argument(
            '--remove_singletons', default=False, action='store_true',
            help='remove singletons in bam file')
    parser.add_argument(
            '--threads', type=int, default=8,
            help='threads for samtools view')
    args = parser.parse_args()

    parse_bam(args.bam, args.mapq, args.NM, args.remove_dup, args.remove_singletons, args.single_end_mapq_filtering, args.threads)


if __name__ == '__main__':
    main()

