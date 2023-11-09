#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-25 15:11

import argparse


def parse_paf(pafs):
    
    for paf in pafs:
        with open(paf) as f:
            n = 0
            for line in f:
                if not line.strip():
                    continue
                n += 1
                cols = line.split()
                query, subject = cols[0], cols[5]
                q_len, s_len = int(cols[1]), int(cols[6])
                if cols[4] == '+':
                    q_start = int(cols[2]) + 1
                    q_end = int(cols[3])
                else:
                    assert cols[4] == '-'
                    q_start = int(cols[3])
                    q_end = int(cols[2]) + 1

                s_start = int(cols[7])
                s_end = int(cols[8])

                print('{}_{}\t{}\t{}\talignment_{}\t{}\t{}'.format(query, subject, q_len, s_len, n, q_start, s_start))
                print('{}_{}\t{}\t{}\talignment_{}\t{}\t{}'.format(query, subject, q_len, s_len, n, q_end, s_end))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pafs', nargs='+', help='the alignment files in PAF format')
    args = parser.parse_args()

    parse_paf(args.pafs)

if __name__ == '__main__':
    main()
