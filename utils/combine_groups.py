#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-04-30 20:59


import sys
import collections

def parse_list(list_file):
    group_dict = collections.defaultdict(list)
    with open(list_file) as f:
        for line in f:
            group_file = line.strip()
            with open(group_file) as fin:
                for l in fin:
                    if l.startswith('#') or not l.strip():
                        continue
                    group_dict[group_file.split('.')[0]].append(l.split()[0])
    print('#Group\tnContigs\tContigs')
    for group, ctg_list in group_dict.items():
        print('{}\t{}\t{}'.format(group, len(ctg_list), ' '.join(ctg_list)))


def main():
    parse_list(sys.argv[1])

if __name__ == '__main__':
    main()
