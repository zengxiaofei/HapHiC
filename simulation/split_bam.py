#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-11 15:04

import argparse
import os


def parse_bam(bam, chrs, threads):

    out_sam_dict = {chr_: open(chr_+'.sam', 'w') for chr_ in chrs}

    with os.popen('samtools view -h {} -@ {}'.format(bam, threads)) as f:
        for line in f:
            if line.startswith('@'):
                if line.startswith('@SQ'):
                    source = line.split()[1].split('_')[0][3:]
                    out_sam_dict[source].write(line)
                else:
                    for sam in out_sam_dict.values():
                        sam.write(line)
            else:
                cols = line.split()
                source = cols[2].split('_')[0]
                if cols[6] == '=' or (cols[6].split('_')[0] == source):
                    out_sam_dict[source].write(line)

    for sam in out_sam_dict.values():
        sam.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='input bam file')
    parser.add_argument('chrs', nargs='+', help='chromosomes to be split')
    parser.add_argument('--threads', type=int, default=8, help='threads for samtools, default: %(default)s')
    args = parser.parse_args()

    parse_bam(args.bam, args.chrs, args.threads)


if __name__ == '__main__':
    main()

