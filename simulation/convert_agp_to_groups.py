#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-10-08 17:05


import argparse
import collections

def parse_agp(agp):

    group_dict = collections.OrderedDict()
    len_dict = dict()
    segment_dict = dict()

    with open(agp) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if cols[4] != 'W':
                continue
            group = cols[0]
            ctg = cols[5]
            if ctg not in segment_dict:
                seg = ctg+'_seg1'
                segment_dict[ctg] = [seg]
            else:
                seg_num = len(segment_dict[ctg])
                seg = '{}_seg{}'.format(ctg, seg_num+1)
                segment_dict[ctg].append(seg)
            
            len_dict[seg] = int(cols[7]) - int(cols[6]) + 1
            
            if group in group_dict:
                group_dict[group].append(seg)
            else:
                group_dict[group] = [seg]
    n = 0
    for group, seg_list in group_dict.items():
        # one group should have at least two contigs inside
        if len(seg_list) < 2:
            continue
        n += 1
        with open('group{}.txt'.format(n), 'w') as fout:
            for seg in seg_list:
                length = len_dict[seg]
                fout.write('{}\tNA\t{}\n'.format(seg, length))

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='input agp file (output of SALSA or YaHS)')
    args = parser.parse_args()

    parse_agp(args.agp)

if __name__ == '__main__':
    main()
