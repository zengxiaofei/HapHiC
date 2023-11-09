#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-04-23 14:33

import argparse
import os

def parse_list_file(listf):
    filtered_read_set = set()
    with open(listf) as f:
        for line in f:
            if not line.strip():
                continue
            filtered_read_set.add(line.strip())
    return filtered_read_set


def filter_bam(bam, filtered_read_set, threads):
    with os.popen('samtools view -h {} -@ {}'.format(bam, threads)) as f:
        for line in f:
            if line.startswith('@'):
                print(line, end='')
            elif line.split()[0] not in filtered_read_set:
                print(line, end='')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='input bam of sim3C read mapping for filtering')
    parser.add_argument('listf', help='input list file of filtered-out reads (output of generate_sim3C_filter_list.py)')
    parser.add_argument('--threads', type=int, default=8, help='threads for running samtools view, default: %(default)s')
    args = parser.parse_args()

    filtered_read_set = parse_list_file(args.listf)
    filter_bam(args.bam, filtered_read_set, args.threads)


if __name__ == '__main__':
    main()

