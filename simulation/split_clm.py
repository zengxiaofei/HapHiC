#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-10-07 01:32

import argparse
import os

def parse_group(groups):

    ctg_group_dict = dict()
    
    for group_file in groups:
        with open(group_file) as f:
            group = os.path.basename(group_file).rsplit('.')[0]
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                cols = line.split()
                ctg_group_dict[cols[0]] = group
    
    return ctg_group_dict


def parse_clm(clm, ctg_group_dict):

    os.mkdir('split_clms')
    fp_dict = {group: open('split_clms/{}.clm'.format(group), 'w') for group in set(ctg_group_dict.values())}
    
    with open(clm) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            ctg1, ctg2 = cols[0][:-1], cols[1][:-1]
            if ctg1 in ctg_group_dict and ctg2 in ctg_group_dict and ctg_group_dict[ctg1] == ctg_group_dict[ctg2]:
                fp_dict[ctg_group_dict[ctg1]].write(line)
    
    for fp in fp_dict.values():
        fp.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('clm', help='input clm file')
    parser.add_argument('groups', nargs='+', help='group files')

    args = parser.parse_args()

    ctg_group_dict = parse_group(args.groups)

    parse_clm(args.clm, ctg_group_dict)

if __name__ == '__main__':
    main()

