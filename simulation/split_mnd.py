#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-11 16:12

import argparse


def parse_mnd(mnd, chrs):

    out_mnd_dict = {chr_: open(chr_+'.mnd', 'w') for chr_ in chrs}

    with open(mnd) as f:
        for line in f:
            cols = line.split()
            if cols[1].split('_')[0] == cols[5].split('_')[0]:
                source = cols[1].split('_')[0]
                out_mnd_dict[source].write(line)

    for mnd in out_mnd_dict.values():
        mnd.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mnd', help='input mnd file')
    parser.add_argument('chrs', nargs='+', help='chromosomes to be split')
    args = parser.parse_args()

    parse_mnd(args.mnd, args.chrs)


if __name__ == '__main__':
    main()

