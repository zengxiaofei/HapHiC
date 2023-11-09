#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-10-12 11:16


import argparse
import os


def parse_bam(bam, threads):
    pass



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'bam', help='input name sorted bam file')
    parser.add_argument(
            '--threads', type=int, default=4,
            help='threads for samtools view')
    args = parser.parse_args()

    parse_bam(args.bam, args.threads)


if __name__ == '__main__':
    main()

