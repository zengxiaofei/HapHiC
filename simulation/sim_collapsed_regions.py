#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-04-25 11:06

import argparse
import random
from decimal import Decimal
import math
from itertools import combinations
import numpy as np
from portion import closed, empty
import collections
import os
import sys


def parse_fasta(fasta):
    fa_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                fa_dict[ID] = [[], 0]
            else:
                seq = line.strip()
                fa_dict[ID][0].append(seq.upper())
    for ID in fa_dict:
        fa_dict[ID][0] = ''.join(fa_dict[ID][0])
        fa_dict[ID][1] = len(fa_dict[ID][0])
    
    return fa_dict


def parse_allele_info(allele_info):

    allele_info_dict = collections.defaultdict(dict)
    
    with open(allele_info) as f:
        # skip the first line
        ploidy = len(f.readline().split()[4:])
        for line in f:
            cols = line.split()
            allele_info_dict[cols[1]][int(cols[2])] = cols[3:]

    return allele_info_dict, ploidy


def simulate_collapsed_regions(fa_dict, collapsed_len, collapsed_ctg_N, weights, ploidy, seed, strict):
    
    # (1) get collapsed regions

    # number of collapsed cadidates in each chromosome
    collapsed_candidate_list = [(ID, seq_len-collapsed_len+1) for ID, (seq, seq_len) in fa_dict.items()]
    # number of all collapsed candidates
    collapsed_candidate_N = sum([n for ID, n in collapsed_candidate_list])
    
    # get some redundant regions to prevent overlapping, not robust
    random.seed(seed)
    collapsed_region_list = random.sample(range(collapsed_candidate_N), collapsed_ctg_N*50)
    random.seed(seed*2)
    random.shuffle(collapsed_region_list)
    
    collapsed_interval_dict = collections.defaultdict(empty)
    collapsed_region_N = 0
    for index in collapsed_region_list:
        # use accumulated_index to find the corresponding chromosome
        accumulated_index = 0
        for ID, n in collapsed_candidate_list:
            accumulated_index += n
            if accumulated_index-1 >= index:
                break
        # the region / index in the corresponding chromosome
        region_start = index-(accumulated_index-n)+1
        region_end = index-(accumulated_index-n)+collapsed_len
        if collapsed_interval_dict[ID].overlaps(closed(region_start, region_end)):
            continue
        else:
            collapsed_region_N += 1
            collapsed_interval_dict[ID] |= closed(region_start, region_end)
        if collapsed_region_N == collapsed_ctg_N:
            break
    if strict:
        assert collapsed_region_N == collapsed_ctg_N

    # (2) get the collapse type of each region based on weights
    
    # two-haplotype collapse C(4, 2) = 6 (1~6)
    # three-haplotype collapse  C(4, 3) = 4 (7~10)
    # four-haplotype collapse C(4, 4) = 1 (11)
    two_hap_ncomb = math.comb(ploidy, 2) if weights[0] else 0
    three_hap_ncomb = math.comb(ploidy, 3) if weights[1] else 0
    four_hap_ncomb = math.comb(ploidy, 4) if weights[2] else 0
    
    elements = list(range(1, two_hap_ncomb+three_hap_ncomb+four_hap_ncomb+1))
    probabilities = list()
    if two_hap_ncomb:
        probabilities.extend([weights[0]/two_hap_ncomb]*two_hap_ncomb)
    if three_hap_ncomb:
        probabilities.extend([weights[1]/three_hap_ncomb]*three_hap_ncomb)
    if four_hap_ncomb:
        probabilities.extend([weights[2]/four_hap_ncomb]*four_hap_ncomb)
    print('elements: {}\nprobabilities: {}'.format(str(elements), str(probabilities)), file=sys.stderr)
    
    np.random.seed(seed*3)
    collapse_type_list = np.random.choice(elements, collapsed_ctg_N, p=probabilities).tolist()

    # (3) get the hap orientation in each collapsed region
    elements = ['+', '-']
    probabilities = [0.5, 0.5]
    np.random.seed(seed*4)
    orientation_list = np.random.choice(elements, collapsed_ctg_N*ploidy, p=probabilities).tolist()

    return collapsed_interval_dict, collapse_type_list, orientation_list, (two_hap_ncomb, three_hap_ncomb, four_hap_ncomb)


def build_output_fasta(collapsed_len, collapsed_ratio, fa_dict, allele_info_dict, collapsed_interval_dict, collapse_type_list, orientation_list, ploidy, ncombs, seed):

    # get complementary base
    base_plus = 'ATCGN'
    base_minus = 'TAGCN'
    com_tab = str.maketrans(base_plus, base_minus)
    
    def get_frag_seq(seq, lower, upper, orient):
        
        if orient == '+':
            # lower & upper is 1-based, but the index of seq is 0-based,
            # so the left boundary should - 1
            return seq[lower-1: upper]
        else:
            # same region but reverse complement
            return seq[lower-1: upper].translate(com_tab)[::-1]


    def output_seq(frag_seq, lower, upper, orient, ID, nhap, f=None):
        
        output_bases = list()
        for m, base in enumerate(frag_seq):
            if orient == '+':
                # nbase is the actual index of the base in the father chromsome
                # nbase, lower, upper are 1-based, m is 0-based
                nbase = lower+m
            else:
                nbase = upper-m
            # when bases are not the same in 4 haps
            if nbase in allele_info_dict[ID]:
                # haps in collapsed_haps are 1-based
                allele_base = allele_info_dict[ID][nbase][nhap]
                # not a deletion
                if allele_base != '-':
                    if orient == '+':
                        output_bases.append(allele_base)
                    else:
                        output_bases.append(allele_base.translate(com_tab)[::-1])
            # bases are same in 4 haps
            else:
                output_bases.append(base)
        
        output_str = ''.join(output_bases)
        # for fcol
        if f:
            f.write('{}\n'.format(output_str))
        # for ftmp
        else:
            return output_str

    # collapse_type_dict = {
    #         1: [1, 2], 2: [1, 3], 3: [1, 4], 4: [2, 3], 5: [2, 4], 6: [3, 4],
    #         7: [1, 2, 3], 8: [1, 2, 4], 9: [1, 3, 4], 10: [2, 3, 4], 
    #         11: [1, 2, 3, 4]
    #         }

    all_haps = list(range(1, ploidy+1))
    
    collapse_types = list()
    if ncombs[0]:
        collapse_types.extend(list(combinations(all_haps, 2)))
    if ncombs[1]:
        collapse_types.extend(list(combinations(all_haps, 3)))
    if ncombs[2]:
        collapse_types.extend(list(combinations(all_haps, 4)))

    collapse_type_dict = {n: t for n, t in enumerate(collapse_types, 1)}
    
    print('collapse_type_dict: {}'.format(str(collapse_type_dict)), file=sys.stderr)

    # generate two genomes:
    # one is collapsed used for downstream anlaysis, 
    # another is chromosome-level and not collapsed used for Hi-C reads simulation
    collapsed_out = 'haplotypes_collapsed_{}_{}.fa'.format(collapsed_len, collapsed_ratio)
    template_out = 'haplotypes_template_{}_{}.fa'.format(collapsed_len, collapsed_ratio)
    
    with open(collapsed_out, 'w') as fcol, open(template_out, 'w') as ftmp:
        
        # index in all collalpsed intervals
        index = 0
        for nchr, (ID, intervals) in enumerate(collapsed_interval_dict.items()):
            
            seq, seq_len = fa_dict[ID]
            father_id = ID.split('_')[0]            
            
            # split and sort intervals
            total_interval = closed(1, seq_len)
            non_collapsed_intervals = total_interval - intervals
            all_intervals = [i for i in intervals] + [i for i in non_collapsed_intervals]
            all_intervals.sort(key=lambda x: x.lower)
            
            # randomly assign orientation for each non-collapsed intervals
            elements = ['+', '-']
            probabilities = [0.5, 0.5]
            np.random.seed(seed*5+nchr)
            non_collapsed_orientation_list = np.random.choice(elements, len(non_collapsed_intervals)*ploidy, p=probabilities).tolist()
            
            # index of non-collapsed intervals in a chromosome 
            index_nc = 0

            # a list to store strings for ftmp output
            ftmp_list = [[] for p in range(ploidy)]
            
            for p in range(1, ploidy+1):
                # get chromosome ID
                ftmp_list[p-1].append('>{}_{}\n'.format(father_id, p))

            for n, i in enumerate(all_intervals):
                
                # get the true lower and upper
                lower = i.lower if str(i.left) == 'CLOSED' else i.lower+1
                upper = i.upper if str(i.right) == 'CLOSED' else i.upper-1
                
                tmp_seq = get_frag_seq(seq, lower, upper, '+')

                # collapsed
                if str(i.left) == str(i.right) == 'CLOSED':
                    # get the randomly assigned collapse type
                    collapse_type = collapse_type_list[index]
                    # translate collapse type to haps
                    collapsed_haps = collapse_type_dict[collapse_type]
                    collapsed_haps_str = ''.join([str(h) for h in collapsed_haps])
                    # the orientation of collapsed region is designated to the orientation of the first hap arbitrarily
                    orient = orientation_list[index*ploidy+collapsed_haps[0]]
                    frag_seq = get_frag_seq(seq, lower, upper, orient)
                    
                    # output collapsed region
                    
                    # write sequence ID of collapsed region
                    # father_chromID, hap, ori_lower, ori_upper, ori_lower, ori_upper, order_num_in_chrom, orientation
                    fcol.write('>{}_{}_{}_{}_collapsed_ctg_{}_{}\n'.format(father_id, collapsed_haps_str, lower, upper, n+1, orient))

                    # write sequence of collapsed region
                    # Using the base from the first haplotype (collapsed_haps[0]), to make the sequence divergence 
                    # between collapsed region and remaining region consistent to that between haplotypes.
                    output_seq(frag_seq, lower, upper, orient, ID, collapsed_haps[0], fcol)
                    
                    # for ftmp
                    for p in collapsed_haps:
                        ftmp_list[p-1].append(output_seq(tmp_seq, lower, upper, '+', ID, collapsed_haps[0]))
                    
                    # output remaining haplotypes
                    
                    remaining_haps = set(all_haps) - set(collapsed_haps)
                    
                    if remaining_haps:
                        for p in remaining_haps:
                            # hap in orientation_list is 0-based
                            # but p is 1-based, so - 1
                            orient = orientation_list[index*ploidy+p-1]
                            frag_seq = get_frag_seq(seq, lower, upper, orient) 
                            
                            # write sequence ID of remaining haplotypes
                            fcol.write('>{}_{}_{}_{}_remaining_hap_{}_{}\n'.format(father_id, p, lower, upper, n+1, orient))

                            # write sequences of remaining haplotypes, p is 1-based
                            output_seq(frag_seq, lower, upper, orient, ID, p, fcol)

                            # for ftmp
                            ftmp_list[p-1].append(output_seq(tmp_seq, lower, upper, '+', ID, p))
                    index += 1
                
                # non-collapsed
                else:
                    for p in all_haps:
                        # p is 1-based
                        orient = non_collapsed_orientation_list[index_nc*ploidy+p-1]
                        frag_seq = get_frag_seq(seq, lower, upper, orient)
                        
                        # write sequence ID of non-collapsed region
                        # father_chromID, hap, ori_lower, ori_upper, order_num_in_chrom, orientation
                        fcol.write('>{}_{}_{}_{}_{}_{}\n'.format(father_id, p, lower, upper, n+1, orient))
                        
                        # write sequences of non-collapsed regions
                        output_seq(frag_seq, lower, upper, orient, ID, p, fcol)

                        # for ftmp
                        ftmp_list[p-1].append(output_seq(tmp_seq, lower, upper, '+', ID, p))
                    index_nc += 1
            
            # output ftmp
            for lines in ftmp_list:
                ftmp.writelines(lines)
                ftmp.write('\n')

            ftmp_list = [[] for p in range(ploidy)]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('allele_info', help='input allele_info.txt')
    parser.add_argument('--collapsed_len', default=500000, type=int, help='length of collapsed region, default: %(default)s')
    parser.add_argument('--collapsed_ratio', default=0.2, type=float, help='length ratio of collapsed regions in the final output genome, default: %(default)s')
    parser.add_argument('--weights', default='0.7,0.2,0.1', type=str, help='weights for two/three/four-haplotype collapse, separated with commas, default: %(default)s')
    parser.add_argument('--seed', default=12345, type=int, help='seed for random processes, default: %(default)s')
    parser.add_argument('--strict', default=False, action='store_true', help='strict collapsed ratio, default: %(default)s')

    args = parser.parse_args()
    
    # parse input files
    fa_dict = parse_fasta(args.fasta)
    allele_info_dict, ploidy = parse_allele_info(args.allele_info)

    # calculate collapsed_ctg_N:
    # collapsed_ratio  = collapsed_len * collapsed_ctg_N / (single_hap_len * 4 - collapsed_len * collapsed_ctg_N * (weights[0]*1 + weights[1]*2 + weights[2]*3))
    single_hap_len = sum([seq_len for seq, seq_len in fa_dict.values()])
    weights = [float(w) for w in args.weights.split(',')]
    if sum(Decimal(w) for w in args.weights.split(',')) != 1:
        weights = [w/sum(weights) for w in weights]
    collapsed_ctg_N = int(single_hap_len*ploidy/((1/args.collapsed_ratio+(weights[0]*1+weights[1]*2+weights[2]*3))*args.collapsed_len))

    # simulate collapsed regions
    collapsed_interval_dict, collapse_type_list, orientation_list, ncombs = simulate_collapsed_regions(
            fa_dict, args.collapsed_len, collapsed_ctg_N, weights, ploidy, args.seed, args.strict)

    # build output fasta file
    build_output_fasta(
            args.collapsed_len, args.collapsed_ratio, fa_dict, allele_info_dict, collapsed_interval_dict, 
            collapse_type_list, orientation_list, ploidy, ncombs, args.seed)


if __name__ == '__main__':
    main()

