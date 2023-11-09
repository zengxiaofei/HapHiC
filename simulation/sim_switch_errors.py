#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-05-18 15:11


import argparse
import sys
import os
import random
import numpy as np

def parse_fasta(fasta):
    
    seq_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                seq_dict[ID] = list()
            else:
                seq_dict[ID].append(line.strip().upper())
    for ID, seq_list in seq_dict.items():
        seq_dict[ID] = ''.join(seq_list)

    return seq_dict


def parse_allele(allele_info):

    allele_list = list()
    with open(allele_info) as f:
        # skip the header line
        f.readline()
        for line in f:
            cols = line.split()
            chrom = cols[1]
            ref_pos = int(cols[2])
            base_list = cols[3:]
            ploidy = len(base_list) - 1
            allele_list.append([chrom, ref_pos, base_list])

    return allele_list, ploidy


def simulate_switch_error(allele_list, ploidy, rate, seed):
    
    alleles_num = len(allele_list)
    error_alleles_num = int(alleles_num * rate)
    
    random.seed(seed)
    error_allele_list = random.sample(allele_list, error_alleles_num)

    elements = list(range(ploidy-1))
    probabilities = [1/(ploidy-1)] * (ploidy-1)
    np.random.seed(seed*2)
    for n, e in enumerate(np.random.choice(elements, error_alleles_num, p=probabilities)):
        error_allele_list[n].append(e)


def output_fasta(fasta, seq_dict, allele_list, ploidy):
  
    trans_dict = dict()
    for n in range(ploidy):
        _list = list(range(ploidy))
        _list.pop(n)
        trans_dict[n] = _list

    allele_dict, error_allele_dict = dict(), dict()
    for _list in allele_list:
        if len(_list) == 3:
            chrom, ref_pos, base_list = _list
            allele_dict[(chrom, ref_pos)] = base_list
        elif len(_list) == 4:
            chrom, ref_pos, base_list, num = _list
            error_allele_dict[(chrom, ref_pos)] = (base_list, num)

    base_name = os.path.basename(fasta).rsplit('.', 1)[0]
    fp_list = [open('{}_hap{}.fa'.format(base_name, n+1), 'w') for n in range(ploidy)]
    fnew = open('new_allele_info.txt', 'w')
    fnew.write('Number\tChrom\tRef_coord\tRef_base\t{}\n'.format('\t'.join(['hap_{}'.format(n+1) for n in range(ploidy)])))
    allele_n = 0
    for ID, seq in seq_dict.items():
        # write ID
        for n, fp in enumerate(fp_list):
            fp.write('>{}_{}\n'.format(ID.rsplit('_', 1)[0], n+1))
        # write bases
        for pos, base in enumerate(seq, 1):
            if (ID, pos) in allele_dict:
                allele_n += 1
                # as allele_info.txt
                base_list = allele_dict[(ID, pos)]
                fnew.write('{}\t{}\t{}\t{}\n'.format(allele_n, ID, pos, '\t'.join(base_list)))
                for n, fp in enumerate(fp_list):
                    new_base = base_list[n+1]
                    if new_base != '-':
                        fp.write(new_base)
            elif (ID, pos) in error_allele_dict:
                allele_n += 1
                # switch error
                base_list, num = error_allele_dict[(ID, pos)]
                for n, b in enumerate(base_list[1:]):
                    if b != base_list[0]:
                        break
                m = trans_dict[n][num]
                base_list[n+1], base_list[m+1] = base_list[m+1], base_list[n+1]
                fnew.write('{}\t{}\t{}\t{}\n'.format(allele_n, ID, pos, '\t'.join(base_list)))
                for n, fp in enumerate(fp_list):
                    new_base = base_list[n+1]
                    if new_base != '-':
                        fp.write(new_base)
            else:
                for fp in fp_list:
                    # same
                    fp.write(base)
        # add line break
        for fp in fp_list:
            fp.write('\n')
    
    # close files
    fnew.close()
    for fp in fp_list:
        fp.close()


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input template genome of FASTA format')
    parser.add_argument('allele_info', help='input allele_info.txt')
    parser.add_argument('--rate', type=float, default=0.1, help='switch error rate, default: %(default)s')
    parser.add_argument('--seed', type=int, default=12345, help='random seed, default: %(default)s')
    args = parser.parse_args()

    seq_dict = parse_fasta(args.fasta)
    allele_list, ploidy = parse_allele(args.allele_info)
    simulate_switch_error(allele_list, ploidy, args.rate, args.seed)
    
    output_fasta(args.fasta, seq_dict, allele_list, ploidy)


if __name__ == '__main__':
    main()

