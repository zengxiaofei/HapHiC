#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-06-12 18:08

import argparse
import gzip

def parse_fastq(out_fq, in_fqs, length):
    
    with gzip.open(out_fq, 'wb') as fout:
        for in_fq in in_fqs:
            with gzip.open(in_fq, 'rb') as f:
                for line1 in f:
                    line2 = f.readline()
                    line3 = f.readline()
                    line4 = f.readline()

                    if len(line2.rstrip()) >= length:
                        fout.write(line1)
                        fout.write(line2)
                        fout.write(line3)
                        fout.write(line4)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('out_fq', help='output fastq file (gzipped)')
    parser.add_argument('in_fqs', nargs='+', help='input fastq file(s) (gzipped)')
    parser.add_argument('--length', type=int, default=50000, help='length threshold, default: %(default)s. Contigs shorter than this value will be removed')
    args = parser.parse_args()

    parse_fastq(args.out_fq, args.in_fqs, args.length)

if __name__ == '__main__':
    main()
