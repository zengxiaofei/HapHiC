#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-03-10 18:28

import sys
import argparse
import collections

comp_table = str.maketrans('ATCG', 'TAGC')


def get_fa_dict(genome):

    fa_dict = dict()
    
    with open(genome) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                fa_dict[ID] = list()
            else:
                fa_dict[ID].append(line.strip().upper())

    for ID in fa_dict:
        fa_dict[ID] = ''.join(fa_dict[ID])
    
    return fa_dict


def split_kmers(fa_dict, k):

    kmer_list = list()
    ID_list = list()

    for ID, seq in fa_dict.items():
        print('splitting assembly kmers for {}'.format(ID), file=sys.stderr)
        ID_list.append(ID)
        kmer_list.append(list())
        l = len(seq)
        for i in range(l - k + 1):
            kmer_forward = seq[i:(i+k)]
            if 'N' in kmer_forward:
                continue
            kmer_reverse = kmer_forward.translate(comp_table)[::-1]
            kmer = kmer_forward if kmer_forward < kmer_reverse else kmer_reverse
            kmer_list[-1].append(kmer)
    
    
    return kmer_list, ID_list


def annotate_kmers(ref_fa_dict, k, kmer_list):

    kmer_set = set()
    
    for kmers in kmer_list:
        kmer_set |= set(kmers)
    
    kmer_annotation_dict = dict()

    for ID, seq in ref_fa_dict.items():
        print('annotating assembly kmers based on reference kmers from {}'.format(ID), file=sys.stderr)
        l = len(seq)
        for i in range(l - k + 1):
            kmer_forward = seq[i:(i+k)]
            if 'N' in kmer_forward:
                continue
            kmer_reverse = kmer_forward.translate(comp_table)[::-1]
            kmer = kmer_forward if kmer_forward < kmer_reverse else kmer_reverse
            if kmer in kmer_set:
                if kmer not in kmer_annotation_dict:
                    kmer_annotation_dict[kmer] = collections.defaultdict(int)
                kmer_annotation_dict[kmer][ID] += 1

    return kmer_annotation_dict


def classify_kmers(kmer_list, ID_list, kmer_annotation_dict, bin_size):

    bin_annotation_dict = dict()
    estimated_chr_dict = dict()
    
    for n, kmers in enumerate(kmer_list):
        ID = ID_list[n]
        chr_count_dict = collections.defaultdict(int)
        print('classifying kmers for {}'.format(ID), file=sys.stderr)
        for i, kmer in enumerate(kmers):
            bin_ID = '{}_bin{}'.format(ID, i//bin_size)
            if bin_ID not in bin_annotation_dict:
                bin_annotation_dict[bin_ID] = collections.defaultdict(int)
            # kmers found in ref
            if kmer in kmer_annotation_dict:
                source_dict = kmer_annotation_dict[kmer]
                # single source, unique kmers
                if len(source_dict) == 1:
                    source = list(source_dict.keys())[0]
                    bin_annotation_dict[bin_ID][source] += 1
                    chr_count_dict[source.split('_')[0]] += 1
                # multi-sources
                else:
                    sources = list(source_dict.keys())
                    # kmers shared by homologous chromosomes
                    if all([source.split('_')[0] == sources[-1].split('_')[0] for source in sources]):
                        chr_source = sources[-1].split('_')[0]
                        bin_annotation_dict[bin_ID][chr_source+'_shared'] += 1
                        chr_count_dict[chr_source] += 1
                    # nonspecific kmers
                    else:
                        bin_annotation_dict[bin_ID]['nonspecific'] += 1
            # kmers that are not found in ref
            else:
                bin_annotation_dict[bin_ID]['unknown'] += 1
        estimated_chr = sorted(list(chr_count_dict.items()), key=lambda x: x[1])[-1][0]
        estimated_chr_dict[ID] = estimated_chr

    return bin_annotation_dict, estimated_chr_dict


def output(asm_fa_dict, bin_annotation_dict, estimated_chr_dict, k, bin_size):
    
    fp_dict = dict()
    for group in asm_fa_dict:
        fp_dict[group] = open('{}_k{}_{}.txt'.format(group, k, bin_size), 'w')
    
    for bin_ID, source_dict in bin_annotation_dict.items():
        group, bin_ = bin_ID.split('_bin')
        start = int(bin_) * bin_size + 1
        end = start + bin_size - 1
        
        chr_specific_n = 0
        stat_dict = collections.defaultdict(int)
        for source, n in source_dict.items():
            # reliable kmers
            if source.endswith('_shared'):
                if source.split('_shared')[0] == estimated_chr_dict[group]:
                    stat_dict['shared'] += n
                    chr_specific_n += n
                else:
                    stat_dict['other_chrom'] += n
            elif source.startswith('chr'):
                chr_source, hap = source.split('_')
                hap = 'hap{}'.format(hap)
                if chr_source == estimated_chr_dict[group]:
                    stat_dict[hap] += n
                    chr_specific_n += n
                else:
                    stat_dict['other_chrom'] += n
            # nonspecific and unknown kmers
            else:
                stat_dict['unreliable'] += n
        stat_list = sorted(list(stat_dict.items()), key=lambda x: x[1])
        primary_source = stat_list[-1][0]
        if primary_source in {'shared', 'hap1', 'hap2', 'hap3', 'hap4'}:
            max_n = 0
            for s, n in stat_list:
                if s.startswith('hap'):
                    if n > max_n:
                        max_n = n
                        primary_source = s
            alpha = max_n / chr_specific_n
        else:
            alpha = 1
        fp_dict[group].write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(group, start, end, primary_source, alpha, source_dict.items()))

    for group in asm_fa_dict:
        fp_dict[group].close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref', help='reference genome')
    parser.add_argument('asm', help='assembly genome')
    parser.add_argument('--kmer_size', type=int, default=201, help='kmer size, default: %(default)s')
    parser.add_argument('--bin_size', type=int, default=500000, help='bin size, default: %(default)s')
    args = parser.parse_args()

    ref_fa_dict = get_fa_dict(args.ref)
    asm_fa_dict = get_fa_dict(args.asm)
    
    k = args.kmer_size
    kmer_list, ID_list = split_kmers(asm_fa_dict, k)
    kmer_annotation_dict = annotate_kmers(ref_fa_dict, k, kmer_list)
    bin_annotation_dict, estimated_chr_dict = classify_kmers(kmer_list, ID_list, kmer_annotation_dict, args.bin_size)
    output(asm_fa_dict, bin_annotation_dict, estimated_chr_dict, k, args.bin_size)


if __name__ == '__main__':
    main()

