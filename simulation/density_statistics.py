#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-09-07 15:55


import argparse


def parse_density(density_file):
    
    with open(density_file) as f:
        for line in f:
            cols = line.strip().split('\t')
            tag, ctg, density = cols
            nhaps = len(ctg.split('_')[1])
            print('{}\t{}\tnhap{}'.format(tag, density, nhaps))
            # print('{}\t{}\ttotal'.format(tag, density))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('density', help='link_density.txt')
    args = parser.parse_args()

    parse_density(args.density)

if __name__ == '__main__':
    main()
