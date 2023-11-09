#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-04-27 16:53

import argparse
import os
import collections

parser = argparse.ArgumentParser()
parser.add_argument('wrk_dir', help='path to wrk_dir')
parser.add_argument('--for_manual', action='store_true', default=False, help='for manual result, default: %(default)s')

args = parser.parse_args()

if args.for_manual:
    rescue_dir = '05.rescue_manual'
else:
    rescue_dir = '03.rescue'

cluster_dict = collections.defaultdict(list)

for root, dirs, files in os.walk(args.wrk_dir):
    basename = os.path.basename(root)
    if basename == rescue_dir:
        for f in files:
            if f.startswith('group'):
                fpath = os.path.join(root, f)
                with open(fpath) as fin:
                    for line in fin:
                        if not line.strip() or line.startswith('#'):
                            continue
                        cols = line.split()
                        group_name = os.path.splitext(f)[0]
                        group_name = '{}_{}'.format(root.split('/')[-2], group_name)
                        cluster_dict[group_name].append(cols[0])


with open('user-prepared.clusters.txt', 'w') as fout:
    fout.write('#Group\tnContigs\tContigs\n')
    for group_name, ctg_list in cluster_dict.items():
        fout.write('{}\t{}\t{}\n'.format(group_name, len(ctg_list), ' '.join(ctg_list)))
