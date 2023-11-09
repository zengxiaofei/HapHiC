#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-09 14:54


import argparse
import random


def random_sampling(npairs, proportion, seed):
    
    sample_n = int(npairs * proportion)
    random.seed(seed)
    return set(random.sample(range(npairs), sample_n))


def parse_mnd(mnd, random_set):

    with open(mnd) as f:
        for n, line in enumerate(f):
            if n in random_set:
                print(line, end='')


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('mnd', help='input "merged_nodups.txt" file')
    parser.add_argument('nparis', type=int, help='number of read pairs in input mnd file')
    parser.add_argument('proportion', type=float, help='proportion of read pairs for sampling')
    parser.add_argument('--seed', type=int, default=12345, help='random seed, default: %(default)s')
    args = parser.parse_args()

    random_set = random_sampling(args.nparis, args.proportion, args.seed)

    parse_mnd(args.mnd, random_set)


if __name__ == '__main__':
    main()

