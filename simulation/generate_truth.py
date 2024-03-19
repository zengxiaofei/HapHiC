#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-03-14 14:08


import argparse
import collections
import re




def is_chr(ID, chr_pattern):

    print(chr_pattern, ID, re.match(chr_pattern, ID))

    if re.match(chr_pattern, ID):
        return True
    else:
        return False


def parse_fasta_for_sim(fasta, chr_pattern):

    truth_dict = collections.defaultdict(list)
    ctg_len_dict = collections.defaultdict(int)
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                # for contigs broken by HapHiC
                ids = ID.rsplit('_', 5)
                source = ids[0] if is_chr(ids[0], chr_pattern) else 'unanchored'
                order = int(ids[1])
                orient = ids[4]
                truth_dict[source].append((ID, orient, order, 0, 0))
            else:
                ctg_len_dict[ID] += len(line.strip())

    return truth_dict, ctg_len_dict


def parse_fasta(fasta, chr_pattern):
    
    truth_dict = collections.defaultdict(list)
    ctg_len_dict = collections.defaultdict(int) 
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                # for contigs broken by HapHiC
                if 'break' in ID:
                    if ID.count('break') == 1:
                        ids = ID.rsplit('_', 3)
                        source = ids[0] if is_chr(ids[0], chr_pattern) else 'unanchored'
                        order = int(ids[1][3:])
                        orient = ids[2]
                        breakn1 = int(ids[3][5:]) if orient == '+' else -int(ids[3][5:])
                        truth_dict[source].append((ID, orient, order, breakn1, 0))
                    else:
                        # for all experiments, the parameter --corrent_nrounds was set to 2
                        assert ID.count('break') == 2
                        ids = ID.rsplit('_', 4)
                        source = ids[0] if is_chr(ids[0], chr_pattern) else 'unanchored'
                        order = int(ids[1][3:])
                        orient = ids[2]
                        breakn1 = int(ids[3][5:]) if orient == '+' else -int(ids[3][5:])
                        breakn2 = int(ids[4][5:]) if orient == '+' else -int(ids[4][5:])
                        truth_dict[source].append((ID, orient, order, breakn1, breakn2))
                # for contigs broken by ALLHiC
                elif ID[-1].isnumeric():
                    ids = ID.rsplit('_', 4)
                    source = ids[0] if is_chr(ids[0], chr_pattern) else 'unanchored'
                    order = int(ids[1][3:])
                    orient = ids[2]
                    start = int(ids[3]) if orient == '+' else -int(ids[3])
                    truth_dict[source].append((ID, orient, order, start, 0))
                # for uncorrected contigs
                else:
                    ids = ID.rsplit('_', 2)
                    source = ids[0] if is_chr(ids[0], chr_pattern) else 'unanchored'
                    order = int(ids[1][3:])
                    orient = ids[2]
                    truth_dict[source].append((ID, orient, order, 0, 0))
            else:
                ctg_len_dict[ID] += len(line.strip())

    return truth_dict, ctg_len_dict


def generate_truth_file(truth_dict):

    with open('truth.txt', 'w') as f:
        for source, ctg_info_list in truth_dict.items():
            f.write('>{}\n'.format(source))
            ctg_info_list.sort(key=lambda x: x[2:])
            for ctg_info in ctg_info_list:
                f.write('{} {}\n'.format(ctg_info[0], ctg_info[1]))


def generate_agp_file(truth_dict, ctg_len_dict, gap_len):
    
    with open('truth.agp', 'w') as f:
        for source, ctg_info_list in truth_dict.items():
            accumulated_len = 0
            # ctg_info_list has been sorted already
            for n, ctg_info in enumerate(ctg_info_list):
                if source != 'unanchored' and n != 0:
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        source, accumulated_len + 1, accumulated_len + gap_len,
                        2 * n, 'U', gap_len, 'contig', 'yes', 'map'))
                    accumulated_len += gap_len
                ctg = ctg_info[0]
                ctg_len = ctg_len_dict[ctg]
                if source == 'unanchored':
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        ctg, accumulated_len + 1, accumulated_len + ctg_len, 
                        1, 'W', ctg, 1, ctg_len, ctg_info[1]))
                    accumulated_len = 0
                else:
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        source, accumulated_len + 1, accumulated_len + ctg_len, 
                        2 * (n + 1) - 1, 'W', ctg, 1, ctg_len, ctg_info[1]))
                    accumulated_len += ctg_len


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='contigs in fasta file')
    parser.add_argument('chr_pattern', help='pattern for chromosome IDs')
    parser.add_argument('--gap_len', type=int, default=100, help='gap length for agp file, default: %(default)s')
    parser.add_argument('--for_simulated_assembly', default=False, action='store_true', help='contig names are parsed in a another way (simulated assembly)')
    args = parser.parse_args()
    
    if args.for_simulated_assembly:
        truth_dict, ctg_len_dict = parse_fasta_for_sim(args.fasta, args.chr_pattern)
    else:
        truth_dict, ctg_len_dict = parse_fasta(args.fasta, args.chr_pattern)
    generate_truth_file(truth_dict)
    generate_agp_file(truth_dict, ctg_len_dict, args.gap_len)

if __name__ == '__main__':
    main()

