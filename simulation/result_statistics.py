#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-04-03 20:27

import sys
import argparse
import collections


def parse_fasta(fasta):
    
    fa_len_dict = dict()
    ignore = False
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                if 'collapsed' in line or 'chimeric' in line:
                    ignore = True
                    continue
                else:
                    ignore = False
                    fa_len_dict[ID] = 0
            elif not ignore:
                fa_len_dict[ID] += len(line.strip())

    return fa_len_dict


def parse_groups(groups):

    n_groups = len(groups)
    anchored_len_dict = collections.defaultdict(int)
    inter_homo_err = 0
    inter_nonhomo_err = 0
    
    chr4_8_excluded_anchored_len = 0
    chr4_8_excluded_inter_homo_err = 0
    chr4_8_excluded_inter_nonhomo_err = 0

    largest_group_dict = collections.defaultdict(int)

    for group in groups:
        
        source_chr_len_dict = collections.defaultdict(int)
        chr4_8_excluded_group_len = 0

        nctgs = 0
        with open(group) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                nctgs += 1

        # one group should have at least two contigs inside
        if nctgs < 2:
            print('group file {} is skipped because of {} contig inside'.format(
                group, nctgs), file=sys.stderr)
            continue

        with open(group) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                # collapsed and chimeric contigs are skipped
                if 'collapsed' in line or 'chimeric' in line:
                    continue
                cols = line.split()
                ctg, length = cols[0], int(cols[2])
                source_chr = '_'.join(ctg.split('_')[:2])
                source_chr_len_dict[source_chr] += length
                anchored_len_dict[source_chr] += length
                if source_chr.split('_')[0] not in {'Chr4', 'Chr8'}:
                    chr4_8_excluded_group_len += length
        
        source_chr_len_list = list(source_chr_len_dict.items())
        if source_chr_len_list:
            source_chr_len_list.sort(key=lambda x: x[1])
            dominant_chr, dominant_chr_len = source_chr_len_list[-1]
        else:
            continue

        for source_chr, length in source_chr_len_list:
            if length > largest_group_dict[source_chr]:
                largest_group_dict[source_chr] = length
            if source_chr != dominant_chr:
                if source_chr.split('_')[0] == dominant_chr.split('_')[0]:
                    inter_homo_err += length
                else:
                    inter_nonhomo_err += length
        
        if dominant_chr.split('_')[0] in {'Chr4', 'Chr8'}:
            continue
        
        chr4_8_excluded_anchored_len += chr4_8_excluded_group_len
        
        for source_chr, length in source_chr_len_list:
            if source_chr.split('_')[0] in {'Chr4', 'Chr8'} or source_chr == dominant_chr:
                continue
            if source_chr.split('_')[0] == dominant_chr.split('_')[0]:
                chr4_8_excluded_inter_homo_err += length
            else:
                chr4_8_excluded_inter_nonhomo_err += length

    return n_groups, anchored_len_dict, inter_homo_err, inter_nonhomo_err, chr4_8_excluded_anchored_len, chr4_8_excluded_inter_homo_err, chr4_8_excluded_inter_nonhomo_err, largest_group_dict


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file of contigs')
    parser.add_argument('groups', nargs='+', help='input group files')
    args = parser.parse_args()

    fa_len_dict = parse_fasta(args.fasta)

    total_ctg_len = sum(fa_len_dict.values())
    
    n_groups, anchored_len_dict, inter_homo_err, inter_nonhomo_err, chr4_8_excluded_anchored_len, chr4_8_excluded_inter_homo_err, chr4_8_excluded_inter_nonhomo_err, largest_group_dict = parse_groups(args.groups)

    anchored_len = sum(anchored_len_dict.values())

    contiguity = sum([d_chr_len/anchored_len_dict[d_chr] for d_chr, d_chr_len in largest_group_dict.items()]) / len(anchored_len_dict)
    
    chr4_8_excluded_contiguity = sum([d_chr_len/anchored_len_dict[d_chr] for d_chr, d_chr_len in largest_group_dict.items() if d_chr.split('_')[0] not in {'Chr4', 'Chr8'}]) / len([source_chr for source_chr in anchored_len_dict if source_chr.split('_')[0] not in {'Chr4', 'Chr8'}])

    print('Contiguity\t{}'.format(contiguity))
    print('Inter_homo_error_rate\t{}%'.format(inter_homo_err/anchored_len*100))
    print('Inter_nonhomo_error_rate\t{}%'.format(inter_nonhomo_err/anchored_len*100))
    print('Chr4_8_excluded_contiguity\t{}'.format(chr4_8_excluded_contiguity))
    print('Chr4_8_excluded_inter_homo_error_rate\t{}%'.format(chr4_8_excluded_inter_homo_err/chr4_8_excluded_anchored_len*100))
    print('Chr4_8_excluded_inter_nonhomo_error_rate\t{}%'.format(chr4_8_excluded_inter_nonhomo_err/chr4_8_excluded_anchored_len*100))
    print('Ngroups\t{}'.format(n_groups))
    print('Anchoring rate\t{}%'.format(anchored_len/total_ctg_len*100))


if __name__ == '__main__':
    main()
