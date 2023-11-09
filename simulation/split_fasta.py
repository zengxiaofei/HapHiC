#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-05-10 09:56


import argparse
import scipy.stats

REVCOM_DICT = {'a': 't', 'A': 'T',
               't': 'a', 'T': 'A',
               'c': 'g', 'C': 'G',
               'g': 'c', 'G': 'C'}


def format_output(ID, ctg, length):
    print('>{}'.format(ID))
    for x in range(length//60+1):
        line_seq = ctg[x*60:(x+1)*60]
        print(line_seq)


def split_fasta(fasta, bin_size):

    output_contig_list = list()
    seq_dict = dict()
    
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                seq_dict[ID] = list()
            else:
                seq_dict[ID].append(line.strip())

    for ID, seq_list in seq_dict.items():
        n = 0
        seq = ''.join(seq_list)
        # make it robust
        seq = seq.replace('n', 'N')
        for ctg in seq.split('N'):
            ctg_len = len(ctg)
            if ctg_len != 0:
                n += 1
                if bin_size:
                    for m in range(ctg_len//(bin_size*1000)+1):
                        sub_seq = ctg[m*bin_size*1000:(m+1)*bin_size*1000]
                        new_ID = '{}_ctg{}_bin{}'.format(ID, n, m+1)
                        output_contig_list.append((new_ID, sub_seq, len(sub_seq)))

                else:
                    new_ID = '{}_ctg{}'.format(ID, n)
                    output_contig_list.append((new_ID, ctg, ctg_len))

    return output_contig_list


def get_orientation(size, seed):
    
    ctg_ori_list = scipy.stats.bernoulli.rvs(0.5, size=size, random_state=seed).tolist()
    
    return ctg_ori_list


def revcom(seq):
    
    seq_list = list()
    for base in seq[::-1]:
        if base in REVCOM_DICT:
            seq_list.append(REVCOM_DICT[base])
        else:
            seq_list.append('N')

    return ''.join(seq_list)


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('--bin_size', default=None, type=int, help='bin size (kbp), default: %(default)s')
    parser.add_argument('--seed', default=12345, type=int, help='random seed, default: %(default)s')
    args = parser.parse_args()

    output_contig_list = split_fasta(args.fasta, args.bin_size)
    ctg_ori_list = get_orientation(len(output_contig_list), args.seed)

    for n, (ID, seq, length) in enumerate(output_contig_list):
        
        if ctg_ori_list[n]:
            new_ID = '{}_{}'.format(ID, '-')
            new_seq = revcom(seq)
        else:
            new_ID = '{}_{}'.format(ID, '+')
            new_seq = seq
        
        format_output(new_ID, new_seq, length)


if __name__ == '__main__':
    main()

