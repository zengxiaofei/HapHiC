#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-03-23 16:48

import argparse
import collections
import scipy.stats


REVCOM_DICT = {'a': 't', 'A': 'T', 
               't': 'a', 'T': 'A',
               'c': 'g', 'C': 'G',
               'g': 'c', 'G': 'C',
               'n': 'n', 'N': 'N'}


def read_fasta(fasta):
    
    len_dict = collections.OrderedDict()
    seq_dict = collections.defaultdict(list)
    
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                len_dict[ID] = 0
            else:
                lstr = line.strip()
                len_dict[ID] += len(lstr)
                seq_dict[ID].append(lstr)

    for ID in seq_dict:
        seq_dict[ID] = ''.join(seq_dict[ID])

    return len_dict, seq_dict


def get_ctg_len_random_list(mean, CV, len_dict, seed, min_len):

    total_len = sum(len_dict.values())
    # * 5 to get some redundant numbers
    n = int(total_len // mean * 5)
    
    ctg_len_list = [int(v) for v in scipy.stats.norm.rvs(loc=mean, scale=CV*mean, size=n, random_state=seed) if int(v) >= min_len]
    ctg_ori_list = scipy.stats.bernoulli.rvs(0.5, size=len(ctg_len_list), random_state=seed).tolist()

    return ctg_len_list, ctg_ori_list


def split_chr(len_dict, ctg_len_list, min_last_len, min_len):

    output_ctg_dict = collections.defaultdict(list)
    
    for ID, chr_len in len_dict.items():
        if chr_len <= ctg_len_list[0]:
            output_ctg_dict[ID].append(chr_len)
        else:
            while ctg_len_list and chr_len > ctg_len_list[0]:
                output_ctg_dict[ID].append(ctg_len_list[0])
                chr_len -= ctg_len_list[0]
                ctg_len_list.pop(0)
            if chr_len and chr_len >= min_last_len and chr_len >= min_len:
                output_ctg_dict[ID].append(chr_len)
                ctg_len_list.pop(0)
            elif chr_len:
                output_ctg_dict[ID][-1] += chr_len

    return output_ctg_dict


def revcom(seq):
    
    seq_list = list()
    for base in seq[::-1]:
        if base in REVCOM_DICT:
            seq_list.append(REVCOM_DICT[base])
        else:
            seq_list.append('N')

    return ''.join(seq_list)


def output_fasta(output_ctg_dict, fasta, len_dict, seq_dict, ctg_ori_list, mean, CV):

    with open('{}.split_{}_{}.fa'.format(fasta.rsplit('.', 1)[0], mean, CV), 'w') as fout:
        m = 0
        for ID, split_len_list in output_ctg_dict.items():
            p = 0
            for n, split_len in enumerate(split_len_list, 1):
                split_seq = seq_dict[ID][p:split_len+p]
                ori = '+'
                if ctg_ori_list[m]:
                    split_seq = revcom(split_seq)
                    ori = '-'
                m += 1
                # format: ChromID_fragNum_start_end_orientation_fragLen
                fout.write('>{}_{}_{}_{}_{}_{}\n{}\n'.format(ID, n, p+1, split_len+p, ori, split_len, split_seq))
                p += split_len


def main():
    
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input chromosome-level genome file')
    parser.add_argument('mean', type=int, help='mean length')
    parser.add_argument('CV', type=float, help='coefficient of variation (standard deviation / mean length)')
    parser.add_argument('--min_len',type=int, default=5000, help='minimum output contig length, however still some contigs longer than 0.5*mean*(1-CV) in chromosome ends may be output, default: %(default)s bp')
    parser.add_argument('--seed', type=int, default=12345, help='seed for random number generation, default: %(default)s')
    args = parser.parse_args()

    len_dict, seq_dict = read_fasta(args.fasta)
    
    ctg_len_list, ctg_ori_list = get_ctg_len_random_list(args.mean, args.CV, len_dict, args.seed, args.min_len)
    
    output_ctg_dict = split_chr(len_dict, ctg_len_list, args.mean*(1-args.CV)/2, args.min_len)

    output_fasta(output_ctg_dict, args.fasta, len_dict, seq_dict, ctg_ori_list, args.mean, args.CV)


if __name__ == '__main__':
    main()

