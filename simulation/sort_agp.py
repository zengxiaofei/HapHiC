#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-05-23 16:00


import argparse
import collections

def parse_agp(agp, chrs):
    
    group_ctg_dict = collections.defaultdict(list)
    group_lines = collections.defaultdict(list)
    with open(agp) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            group_lines[cols[0]].append(line)
            if cols[4] != 'W':
                continue
            source = cols[5].rsplit('_', 2)[0]
            if source in chrs:
                group_ctg_dict[cols[0]].append((cols[5], int(cols[7])))
    
    return group_ctg_dict, group_lines


def get_group_chr(group_ctg_dict):

    group_chr_dict = collections.defaultdict(dict)
    for group, ctg_list in group_ctg_dict.items():
        sorted_ctg_list = sorted(ctg_list, key=lambda x: x[1], reverse=True)
        if len(sorted_ctg_list) >= 10:
            top_ten_ctg_list = sorted_ctg_list[:10]
            for ctg, length in top_ten_ctg_list:
                chrom = ctg.rsplit('_', 2)[0]
                if chrom in group_chr_dict[group]:
                    group_chr_dict[group][chrom] += length
                else:
                    group_chr_dict[group][chrom] = length

    return group_chr_dict


def parse_group(group_chr_dict):
    
    chr_to_groups = collections.defaultdict(list)

    for group, chrom_dict in group_chr_dict.items():
        # print(group, chrom_dict)
        chrom = sorted(chrom_dict.items(), key=lambda x: x[1], reverse=True)[0][0]
        chr_to_groups[chrom].append(group)

    return chr_to_groups


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='agp file output by HapHiC or ALLHiC')
    parser.add_argument('chr_ordering', nargs='+', help='ordering of output chrs')
    args = parser.parse_args()

    group_ctg_dict, group_lines = parse_agp(args.agp, args.chr_ordering)
    group_chr_dict = get_group_chr(group_ctg_dict)
    chr_to_groups = parse_group(group_chr_dict)
    output_group  = set()
    
    # ordered groups
    for chrom in args.chr_ordering:
        for group in chr_to_groups[chrom]:
            output_group.add(group)
            # print(chrom, group)
            for line in group_lines[group]:
                print(line, end='')
    # other groups
    for group, lines in group_lines.items():
        if group not in output_group:
            for line in lines:
                print(line, end='')


if __name__ == '__main__':
    main()
