#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-09-01 11:10

import argparse
import collections

def parse_groups(groups):

    group_ctg_dict = collections.defaultdict(list)

    for n, group in enumerate(groups, 1):
        group_name = 'group{}'.format(n)
        with open(group) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                cols = line.split()
                ctg = cols[0]
                group_ctg_dict[group_name].append(ctg)
    
    return group_ctg_dict


def output_clusters(group_ctg_dict):
    
    with open('all.clusters.txt', 'w') as f:
        f.write('#Group\tnContigs\tContigs\n')
        for group_name, ctgs in group_ctg_dict.items():
            f.write('{}\t{}\t{}\n'.format(group_name, len(ctgs), ' '.join(ctgs)))

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('groups', nargs='+', help='input group files')
    args = parser.parse_args()

    group_ctg_dict = parse_groups(args.groups)
    output_clusters(group_ctg_dict)

if __name__ == '__main__':
    main()

