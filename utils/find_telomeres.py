#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-09-05 11:44

import argparse
import os
import sys
import re


comp_table = str.maketrans('ATCGN', 'TAGCN')


def parse_genome(genome, contigs):
    
    output = False
    fa_dict = dict()
    with open(genome) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                if contigs is None or ID in contigs:
                    fa_dict[ID] = []
                    output = True
                else:
                    output = False
            elif output:
                fa_dict[ID].append(line.strip().upper())

    for ID in fa_dict:
        fa_dict[ID] = ''.join(fa_dict[ID])

    return fa_dict


def revcom(seq):

    return seq[::-1].translate(comp_table)


def get_hist(seq, repeat):

    hist = [0] * 100
    matches = re.finditer(repeat, seq)
    for match in matches:
        hist[match.start() * 100 // len(seq)] += 1

    return ''.join(['{}\t{}\n'.format(i + 1, h) for i, h in enumerate(hist)]), max(hist)


def find_repeat(fa_dict, repeat, contigs, show_uplot):

    if contigs is None:
        ID_list = fa_dict.keys()
    else:
        ID_list = contigs

    reverse_repeat = revcom(repeat)
    forward_two_repeats = repeat * 2
    reverse_two_repeats = revcom(repeat) * 2
    reverse_two_repeats_reverse = reverse_two_repeats[::-1]
    
    print('Seq_ID\tSeq_len\tNumber_of_{0}/{1}\tNumber_of_{0}/{1}_per_Mb\tLeftmost_{0}_pos\tRightmost_{1}_pos\tLeftmost_relative_pos\tRightmost_relative_pos'.format(
        forward_two_repeats, reverse_two_repeats))

    for ID in ID_list:
        seq_len = len(fa_dict[ID])
        nrepeats = fa_dict[ID].count(forward_two_repeats) + fa_dict[ID].count(reverse_two_repeats)

        if forward_two_repeats in fa_dict[ID]:
            start_pos = fa_dict[ID].index(forward_two_repeats) + 1
            relative_start_pos = '{:.4f}'.format(start_pos / seq_len)
        else:
            start_pos, relative_start_pos = 'NA', 'NA'

        if reverse_two_repeats in fa_dict[ID]:
            end_pos = seq_len - (fa_dict[ID][::-1].index(reverse_two_repeats_reverse) + 1)
            relative_end_pos = '{:.4f}'.format(end_pos / seq_len)
        else:
            end_pos, relative_end_pos = 'NA', 'NA'

        print('{}\t{}\t{}\t{:.4f}\t{}\t{}\t{}\t{}'.format(
            ID, seq_len, nrepeats, nrepeats/seq_len*1000000, start_pos, end_pos, relative_start_pos, relative_end_pos))

        if not show_uplot:
            continue

        forward_data, forward_max = get_hist(fa_dict[ID], repeat)
        reverse_data, reverse_max = get_hist(fa_dict[ID], reverse_repeat)
        max_y = max(forward_max, reverse_max)

        os.system('echo -e "{}" | uplot line -w 50 -h 5 -t "{} {} (forward)" --xlabel "Contig position (100 bins)" --ylabel "Repeat count" --xlim 0,100 --ylim 0,{}'.format(
            forward_data, ID, repeat, max_y))
        print('', file=sys.stderr)
        os.system('echo -e "{}" | uplot line -w 50 -h 5 -t "{} {} (reverse)" --xlabel "Contig position (100 bins)" --ylabel "Repeat count" --xlim 0,100 --ylim 0,{}'.format(
            reverse_data, ID, reverse_repeat, max_y))
        print('', file=sys.stderr)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome', help='input FASTA file of genome')
    parser.add_argument('--repeat', default='CCCTAAA', help='telomere repeat, default: %(default)s')
    parser.add_argument('--show_uplot', default=False, action='store_true', help='show lineplot using YouPlot')
    parser.add_argument('--contigs', nargs='+', help='contig IDs, default: %(default)s. If this parameter is not specified, the program will find telomeres in all contigs.')

    args = parser.parse_args()
    repeat = args.repeat.upper()

    fa_dict = parse_genome(args.genome, args.contigs)
    find_repeat(fa_dict, repeat, args.contigs, args.show_uplot)

if __name__ == '__main__':
    main()

