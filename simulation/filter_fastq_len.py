#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-03-25 18:54

import argparse
import gzip


def parse_fastq_files(fastq_files, len_cutoff):
    
    with gzip.open('output.fq.gz', 'wt') as fout:
        for fq in fastq_files:
            if fq.endswith('fastq') or fq.endswith('fq'):
                fopen = open
            else:
                assert fq.endswith('fastq.gz') or fq.endswith('fq.gz')
                fopen = gzip.open
            with fopen(fq, 'rt') as fin:
                for line1 in fin:
                    line2 = fin.readline()
                    line3 = fin.readline()
                    line4 = fin.readline()
                    if len(line2) >= len_cutoff and len(line2) == len(line4):
                        fout.write(line1)
                        fout.write(line2)
                        fout.write(line3)
                        fout.write(line4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_files', nargs='+', help='input FASTQ file(s) for length filtering, gzipped files are acceptable')
    parser.add_argument('--len_cutoff', type=int, default=50000, help='length cutoff, reads shorter than this value will be removed, default: %(default)s')
    args = parser.parse_args()

    parse_fastq_files(args.fastq_files, args.len_cutoff)


if __name__ == '__main__':
    main()

