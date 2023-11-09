#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-06-23 16:59


import argparse
import collections
import numpy as np
import random
import math


def parse_fasta(fasta):

    source_chr_dict = collections.defaultdict(list)
    fa_dict = collections.defaultdict(list)
    
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                source_chr = '_'.join(ID.split('_')[:2])
                source_chr_dict[source_chr].append(ID)
            else:
                fa_dict[ID].append(line.rstrip())

    for ID, seq in fa_dict.items():
        fa_dict[ID] = len(''.join(seq))

    return source_chr_dict, fa_dict


def simulate_errors(source_chr_dict, ratio, error_type, seed):
    
    # prepare
    
    total_n_ctgs = 0
    total_ctg_list = list()
    for source_chrs, ctgs in source_chr_dict.items():
        total_n_ctgs += len(ctgs)
        total_ctg_list.extend(ctgs)
    total_ctg_list.sort()
    error_n_ctgs = int(ratio * total_n_ctgs)
    
    n_chrs = len(source_chr_dict)
    ploidy = max([int(sc.split('_')[1]) for sc in source_chr_dict])

    # simulation
    new_source_chr_dict = collections.defaultdict(list)

    
    # sample contigs to simulate errors
    random.seed(seed)
    error_ctg_list = random.sample(total_ctg_list, error_n_ctgs)
    
    if error_type == 'inter_homo':
        
        # find somewhere to put these contigs in
        elements = list(range(ploidy-1))
        probabilities = [1/(ploidy-1)] * (ploidy-1)
        np.random.seed(seed*2)
        element_list = np.random.choice(elements, error_n_ctgs, p=probabilities)
    
        # generate new source_chr_dict
        for source_chr, ctgs in source_chr_dict.items():
            for ctg in ctgs:
                if ctg in error_ctg_list:
                    other_chrs = ['{}_{}'.format(source_chr.split('_')[0], n+1) for n in range(ploidy)]
                    other_chrs.remove(source_chr)
                    error_ctg_index = error_ctg_list.index(ctg)
                    new_chr_index = element_list[error_ctg_index]
                    new_chr = other_chrs[new_chr_index]
                    new_source_chr_dict[new_chr].append(ctg)
                else:
                    new_source_chr_dict[source_chr].append(ctg)
   
    elif error_type == 'inter_nonhomo':

        # find somewhere to put these contigs in
        elements = list(range(n_chrs-ploidy))
        probabilities = [1/(n_chrs-ploidy)] * (n_chrs-ploidy)
        np.random.seed(seed*2)
        element_list = np.random.choice(elements, error_n_ctgs, p=probabilities)

        # generate new source_chr_dict
        for source_chr, ctgs in source_chr_dict.items():
            for ctg in ctgs:
                if ctg in error_ctg_list:
                    other_chrs = [sc for sc in source_chr_dict if sc.split('_')[0] != source_chr.split('_')[0]]
                    other_chrs.sort()
                    error_ctg_index = error_ctg_list.index(ctg)
                    new_chr_index = element_list[error_ctg_index]
                    new_chr = other_chrs[new_chr_index]
                    new_source_chr_dict[new_chr].append(ctg)
                else:
                    new_source_chr_dict[source_chr].append(ctg)

    elif error_type == 'anchoring_rate':
        
        # just remove these contigs
        for source_chr, ctgs in source_chr_dict.items():
            for ctg in ctgs:
                if ctg not in error_ctg_list:
                    new_source_chr_dict[source_chr].append(ctg)

    return new_source_chr_dict


def simulate_contiguity(source_chr_dict, ratio, error_type, seed):
   
    # prepare
    
    contiguity = ratio
    split_n = math.ceil(1/contiguity)
    new_source_chr_dict = collections.defaultdict(list)

    # simulation
    
    for m, (source_chr, ctgs) in enumerate(source_chr_dict.items()):
        n_ctgs = len(ctgs)
        max_n_ctgs = int(contiguity*n_ctgs)
        _ctgs = sorted(ctgs)
        for n in range(1, split_n):
            random.seed(seed*n+m)
            sampled_ctgs = random.sample(_ctgs, max_n_ctgs)
            for ctg in sampled_ctgs:
                new_source_chr = '{}_{}'.format(source_chr, n)
                new_source_chr_dict[new_source_chr].append(ctg)
                _ctgs.remove(ctg)

        if _ctgs:
            for ctg in _ctgs:
                new_source_chr = '{}_{}'.format(source_chr, n+1)
                new_source_chr_dict[new_source_chr].append(ctg)

    return new_source_chr_dict


def output_clusters(new_source_chr_dict, fa_dict, error_type, ratio, output_groups):
    
    with open('{}_{}.clusters.txt'.format(error_type, ratio), 'w') as fcluster:
        fcluster.write('#Group\tnContigs\tContigs\n')
        for source_chr, ctgs in new_source_chr_dict.items():
            fcluster.write('{}\t{}\t{}\n'.format(source_chr, len(ctgs), ' '.join(ctgs)))
            if output_groups:
                with open('{}_{}.group_{}.txt'.format(error_type, ratio, source_chr), 'w') as fgroup:
                    fgroup.write('#Contig\tRECounts\tLength\n')
                    for ctg in ctgs:
                        fgroup.write('{}\tNA\t{}\n'.format(ctg, fa_dict[ctg]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'fasta', 
            help='input contigs in fasta format')
    parser.add_argument(
            'ratio', type=float,
            help='ratio of errors / contiguity')
    parser.add_argument(
            '--error_type', 
            default='inter_homo',
            choices={'inter_homo', 'inter_nonhomo', 'contiguity', 'anchoring_rate'}, 
            help='error type of the simulation, default: %(default)s')
    parser.add_argument(
            '--output_groups',
            default=False, action='store_true',
            help='output group files, default: %(default)s')
    
    parser.add_argument(
            '--seed',
            default=12345, type=int,
            help='seed for random simulation, default: %(default)s')
    

    args = parser.parse_args()

    # parse fasta file to get the ground truth
    source_chr_dict, fa_dict = parse_fasta(args.fasta)

    if args.error_type != 'contiguity':
        new_source_chr_dict = simulate_errors(source_chr_dict, args.ratio, args.error_type, args.seed)
    else:
        new_source_chr_dict = simulate_contiguity(source_chr_dict, args.ratio, args.error_type, args.seed)

    # output new clusters file
    output_clusters(new_source_chr_dict, fa_dict, args.error_type, args.ratio, args.output_groups)


if __name__ == '__main__':
    main()

