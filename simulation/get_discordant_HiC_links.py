#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-03-15 09:45


import argparse
import collections
from math import ceil
from portion import closed, empty
import pysam
import numpy as np
from scipy import stats

pysam.set_verbosity(0)


def parse_agp(agp, bin_size):

    def change_orientation(cols_list):
        forward = 0
        reverse = 0
        for cols in cols_list:
            if cols[4] != 'W':
                continue
            if cols[8] not in cols[5].split('_'):
                reverse += int(cols[2]) - int(cols[1]) + 1
            else:
                forward += int(cols[2]) - int(cols[1]) + 1
        if reverse > forward:
            oriented_cols_list = []
            scaffold_end = int(cols[2])
            for i, cols in enumerate(cols_list[::-1]):
                start, end = int(cols[1]), int(cols[2])
                cols[2] = scaffold_end - start + 1
                cols[1] = scaffold_end - end + 1
                cols[3] = i + 1
                cols[8] = '+' if cols[8] == '-' else '-'
                oriented_cols_list.append(cols)
            return oriented_cols_list
        else:
            return cols_list
    
    agp_dict = collections.defaultdict(list)
    with open(agp) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            scaffold = cols[0]
            agp_dict[scaffold].append(cols)
    
    bin_index_dict = collections.defaultdict(set)
    ctg_bin_dict = collections.defaultdict(list)
    for scaffold, cols_list in agp_dict.items():
        gap_len = 0
        oriented_cols_list = change_orientation(cols_list)
        for cols in oriented_cols_list:
            scaffold = cols[0]
            if cols[4] == 'W':
                order = int(cols[3])
                # offset is used to eliminate the difference of gap length
                offset = int((order - 1) / 2 * gap_len)
                start, end = int(cols[1]) - offset, int(cols[2]) - offset
                ctg = cols[5]
                orient = cols[8]
                bins_in_scaffold = empty()
                for i in range(ceil(start / bin_size) - 1, ceil(end / bin_size)):
                    bins_in_scaffold |= closed(bin_size * i + 1, bin_size * (i + 1))
                    bin_index_dict[scaffold].add(i)
                for intersection in bins_in_scaffold & closed(start, end):
                    bin_index = ceil(intersection.upper / bin_size) - 1
                    assert ceil(intersection.lower / bin_size) - 1 == bin_index
                    if orient == '+':
                        ctg_range = closed(intersection.lower - start + 1, intersection.upper - start + 1)
                        ctg_bin_dict[ctg].append((ctg_range, scaffold, bin_index))
                    else:
                        ctg_range = closed(end - intersection.upper + 1, end - intersection.lower + 1)
                        ctg_bin_dict[ctg].append((ctg_range, scaffold, bin_index))
            elif cols[4] in ('U', 'N'):
                gap_len = int(cols[5])

    return ctg_bin_dict, bin_index_dict


def parse_bam(bam, truth_ctg_bin_dict, haphic_ctg_bin_dict, haphic_fastsort_ctg_bin_dict, yahs_ctg_bin_dict, threads):
    
    def get_scaffold_and_bin_index(ctg, pos, ctg_bin_dict):

        for ctg_range, scaffold, bin_index in ctg_bin_dict[ctg]:
            if pos in ctg_range:
                return scaffold, bin_index

    truth_link_dict = collections.defaultdict(int)
    haphic_link_dict = collections.defaultdict(int)
    haphic_fastsort_link_dict = collections.defaultdict(int)
    yahs_link_dict = collections.defaultdict(int)
    
    format_options = [b'filter=flag.read1']
    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as f:
        for aln in f:
            ref, mref, pos, mpos = aln.reference_name, aln.next_reference_name, aln.reference_start + 1, aln.next_reference_start + 1
            # truth
            scaffold, bin_index = get_scaffold_and_bin_index(ref, pos, truth_ctg_bin_dict)
            mscaffold, mbin_index = get_scaffold_and_bin_index(mref, mpos, truth_ctg_bin_dict)
            truth_link_dict[tuple(sorted([(scaffold, bin_index), (mscaffold, mbin_index)]))] += 1
            # haphic
            scaffold, bin_index = get_scaffold_and_bin_index(ref, pos, haphic_ctg_bin_dict)
            mscaffold, mbin_index = get_scaffold_and_bin_index(mref, mpos, haphic_ctg_bin_dict)
            haphic_link_dict[tuple(sorted([(scaffold, bin_index), (mscaffold, mbin_index)]))] += 1
            # truth
            scaffold, bin_index = get_scaffold_and_bin_index(ref, pos, haphic_fastsort_ctg_bin_dict)
            mscaffold, mbin_index = get_scaffold_and_bin_index(mref, mpos, haphic_fastsort_ctg_bin_dict)
            haphic_fastsort_link_dict[tuple(sorted([(scaffold, bin_index), (mscaffold, mbin_index)]))] += 1
            # truth
            scaffold, bin_index = get_scaffold_and_bin_index(ref, pos, yahs_ctg_bin_dict)
            mscaffold, mbin_index = get_scaffold_and_bin_index(mref, mpos, yahs_ctg_bin_dict)
            yahs_link_dict[tuple(sorted([(scaffold, bin_index), (mscaffold, mbin_index)]))] += 1

    return truth_link_dict, haphic_link_dict, haphic_fastsort_link_dict, yahs_link_dict


def make_stat(fstat, link_dict, tag, summary_intra_dict, summary_inter_dict):
    
    for ((scaffold1, bin_index1), (scaffold2, bin_index2)), nlinks in link_dict.items():
        # intra-scaffold Hi-C links
        if scaffold1 == scaffold2:
            if tag in summary_intra_dict[bin_index2 - bin_index1]:
                summary_intra_dict[bin_index2 - bin_index1][tag].append(nlinks)
            else:
                summary_intra_dict[bin_index2 - bin_index1][tag] = [nlinks]
            fstat.write('{}\t{}\t{}\t{}\t{}\tintra\t{}\n'.format(scaffold1, bin_index1, scaffold2, bin_index2, nlinks, tag))
        # inter-scaffold Hi-C links
        else:
            summary_inter_dict[tag].append(nlinks)
            fstat.write('{}\t{}\t{}\t{}\t{}\tinter\t{}\n'.format(scaffold1, bin_index1, scaffold2, bin_index2, nlinks, tag))


def summarize(summary_intra_dict, summary_inter_dict, haphic_bin_index_dict, haphic_fastsort_bin_index_dict, yahs_bin_index_dict):
    
    def not_in_interval(nlinks_list, interval, nlinks_not_in_interval, nbins_not_in_interval):
        
        for nlinks in nlinks_list:
            if nlinks not in interval:
                nlinks_not_in_interval += nlinks
                nbins_not_in_interval += 1
        
        return nlinks_not_in_interval, nbins_not_in_interval
    
    def get_diff(nlinks_list, truth_count_dict):
        
        diff_nlinks = 0
        diff_nbins = 0
        count_dict = collections.defaultdict(int)
        
        for nlinks in nlinks_list:
            count_dict[nlinks] += 1

        all_nlinks = set(truth_count_dict.keys()) | set(count_dict.keys())
        
        for nlinks in all_nlinks:
            if count_dict[nlinks] != truth_count_dict[nlinks]:
                d = count_dict[nlinks] - truth_count_dict[nlinks] if count_dict[nlinks] > truth_count_dict[nlinks] else 0
                diff_nbins += d
                diff_nlinks += d * nlinks
        
        return diff_nlinks, diff_nbins
    
    def get_nbins(bin_index_dict):
        
        intra_nbins = 0
        total_scaffold_nbins = 0
        for scaffold, bin_set in bin_index_dict.items():
            scaffold_nbins = len(bin_set)
            total_scaffold_nbins += scaffold_nbins
            intra_nbins += scaffold_nbins + int(scaffold_nbins * (scaffold_nbins - 1) / 2)
        inter_nbins = total_scaffold_nbins + int(total_scaffold_nbins * (total_scaffold_nbins - 1) / 2) - intra_nbins
        
        return intra_nbins, inter_nbins
    
    haphic_nlinks_not_in_interval, haphic_nbins_not_in_interval = 0, 0
    haphic_fastsort_nlinks_not_in_interval, haphic_fastsort_nbins_not_in_interval = 0, 0
    yahs_nlinks_not_in_interval, yahs_nbins_not_in_interval = 0, 0

    haphic_intra_nlinks = 0
    haphic_fastsort_intra_nlinks = 0
    yahs_intra_nlinks = 0
    # intra-scaffold Hi-C links
    for bin_dist in summary_intra_dict:
        truth_nlinks_list = summary_intra_dict[bin_dist]['Truth'] if 'Truth' in summary_intra_dict[bin_dist] else []
        haphic_nlinks_list = summary_intra_dict[bin_dist]['HapHiC'] if 'HapHiC' in summary_intra_dict[bin_dist] else []
        haphic_fastsort_nlinks_list = summary_intra_dict[bin_dist]['HapHiC_fastsort'] if 'HapHiC_fastsort' in summary_intra_dict[bin_dist] else []
        yahs_nlinks_list = summary_intra_dict[bin_dist]['YaHS'] if 'YaHS' in summary_intra_dict[bin_dist] else []
        
        haphic_intra_nlinks += sum(haphic_nlinks_list)
        haphic_fastsort_intra_nlinks += sum(haphic_fastsort_nlinks_list)
        yahs_intra_nlinks += sum(yahs_nlinks_list)
        # get confidence interval
        # df should >= 1
        if len(truth_nlinks_list) >= 2:
            # lower, upper = stats.t.interval(0.95, df=len(truth_nlinks_list)-1, loc=np.mean(truth_nlinks_list), scale=stats.sem(truth_nlinks_list))
            # interval = closed(lower, upper)
            interval = closed(min(truth_nlinks_list), max(truth_nlinks_list))
        elif len(truth_nlinks_list) == 1:
            interval = {truth_nlinks_list[0]}
        else:
            interval = empty()
        
        # print(bin_dist, truth_nlinks_list, interval)
        # print(haphic_nlinks_list)
        haphic_nlinks_not_in_interval, haphic_nbins_not_in_interval = not_in_interval(haphic_nlinks_list, interval, haphic_nlinks_not_in_interval, haphic_nbins_not_in_interval)
        haphic_fastsort_nlinks_not_in_interval, haphic_fastsort_nbins_not_in_interval = not_in_interval(haphic_fastsort_nlinks_list, interval, haphic_fastsort_nlinks_not_in_interval, haphic_fastsort_nbins_not_in_interval)
        yahs_nlinks_not_in_interval, yahs_nbins_not_in_interval = not_in_interval(yahs_nlinks_list, interval, yahs_nlinks_not_in_interval, yahs_nbins_not_in_interval)

    # inter-scaffold Hi-C links
    truth_inter_nlinks_list = summary_inter_dict['Truth']
    haphic_inter_nlinks_list = summary_inter_dict['HapHiC']
    haphic_fastsort_inter_nlinks_list = summary_inter_dict['HapHiC_fastsort']
    yahs_inter_nlinks_list = summary_inter_dict['YaHS']
    
    truth_count_dict = collections.defaultdict(int)
    for nlinks in truth_inter_nlinks_list:
        truth_count_dict[nlinks] += 1
    
    haphic_diff_nlinks, haphic_diff_nbins = get_diff(haphic_inter_nlinks_list, truth_count_dict)
    haphic_fastsort_diff_nlinks, haphic_fastsort_diff_nbins = get_diff(haphic_fastsort_inter_nlinks_list, truth_count_dict)
    yahs_diff_nlinks, yahs_diff_nbins = get_diff(yahs_inter_nlinks_list, truth_count_dict)
    
    haphic_inter_nlinks = sum(haphic_inter_nlinks_list)
    haphic_fastsort_inter_nlinks = sum(haphic_fastsort_inter_nlinks_list)
    yahs_inter_nlinks = sum(yahs_inter_nlinks_list)
    
    haphic_intra_nbins, haphic_inter_nbins = get_nbins(haphic_bin_index_dict)
    haphic_fastsort_intra_nbins, haphic_fastsort_inter_nbins = get_nbins(haphic_fastsort_bin_index_dict)
    yahs_intra_nbins, yahs_inter_nbins = get_nbins(yahs_bin_index_dict)

    print('\n###### HapHiC ######')
    print('\tTotal number of intra-scaffold Hi-C links: {}'.format(haphic_intra_nlinks))
    print('\tTotal number of intra-scaffold bins: {}'.format(haphic_intra_nbins))
    print('\tNumber of differential intra-scaffold Hi-C links: {} ({} %)'.format(
        haphic_nlinks_not_in_interval, haphic_nlinks_not_in_interval/haphic_intra_nlinks*100))
    print('\tNumber of differential intra-scaffold bins: {} ({} %)'.format(
        haphic_nbins_not_in_interval, haphic_nbins_not_in_interval/haphic_intra_nbins*100))
    print('------------------------------------------------------------') 
    print('\tTotal number of inter-scaffold Hi-C links: {}'.format(haphic_inter_nlinks))
    print('\tTotal number of inter-scaffold bins: {}'.format(haphic_inter_nbins))
    print('\tNumber of differential inter-scaffold Hi-C links: {} ({} %)'.format(
        haphic_diff_nlinks, haphic_diff_nlinks/haphic_inter_nlinks*100))
    print('\tNumber of differential inter-scaffold bins: {} ({} %)'.format(
        haphic_diff_nbins, haphic_diff_nbins/haphic_inter_nbins*100))
    
    print('\n###### HapHiC fast sorting ######')
    print('\tTotal number of intra-scaffold Hi-C links: {}'.format(haphic_fastsort_intra_nlinks))
    print('\tTotal number of intra-scaffold bins: {}'.format(haphic_fastsort_intra_nbins))
    print('\tNumber of differential intra-scaffold Hi-C links: {} ({} %)'.format(
        haphic_fastsort_nlinks_not_in_interval, haphic_fastsort_nlinks_not_in_interval/haphic_fastsort_intra_nlinks*100))
    print('\tNumber of differential intra-scaffold bins: {} ({} %)'.format(
        haphic_fastsort_nbins_not_in_interval, haphic_fastsort_nbins_not_in_interval/haphic_fastsort_intra_nbins*100))
    print('------------------------------------------------------------') 
    print('\tTotal number of inter-scaffold Hi-C links: {}'.format(haphic_fastsort_inter_nlinks))
    print('\tTotal number of inter-scaffold bins: {}'.format(haphic_fastsort_inter_nbins))
    print('\tNumber of differential inter-scaffold Hi-C links: {} ({} %)'.format(
        haphic_fastsort_diff_nlinks, haphic_fastsort_diff_nlinks/haphic_fastsort_inter_nlinks*100))
    print('\tNumber of differential inter-scaffold bins: {} ({} %)'.format(
        haphic_fastsort_diff_nbins, haphic_fastsort_diff_nbins/haphic_fastsort_inter_nbins*100))
    print('\n###### YaHS ######')
    print('\tTotal number of intra-scaffold Hi-C links: {}'.format(yahs_intra_nlinks))
    print('\tTotal number of intra-scaffold bins: {}'.format(yahs_intra_nbins))
    print('\tNumber of differential intra-scaffold Hi-C links: {} ({} %)'.format(
        yahs_nlinks_not_in_interval, yahs_nlinks_not_in_interval/yahs_intra_nlinks*100))
    print('\tNumber of differential intra-scaffold bins: {} ({} %)'.format(
        yahs_nbins_not_in_interval, yahs_nbins_not_in_interval/yahs_intra_nbins*100))
    print('------------------------------------------------------------') 
    print('\tTotal number of inter-scaffold Hi-C links: {}'.format(yahs_inter_nlinks))
    print('\tTotal number of inter-scaffold bins: {}'.format(yahs_inter_nbins))
    print('\tNumber of differential inter-scaffold Hi-C links: {} ({} %)'.format(
        yahs_diff_nlinks, yahs_diff_nlinks/yahs_inter_nlinks*100))
    print('\tNumber of differential inter-scaffold bins: {} ({} %)'.format(
        yahs_diff_nbins, yahs_diff_nbins/yahs_inter_nbins*100))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('agp_truth', help='agp file for truth')
    parser.add_argument('agp_haphic', help='agp file for haphic scaffolds')
    parser.add_argument('agp_haphic_fastsort', help='agp file for haphic (fast_sorting) scaffolds')
    parser.add_argument('agp_yahs', help='agp file for yahs scaffolds')
    parser.add_argument('bam', help='bam file for mapped Hi-C reads')
    parser.add_argument('--bin_size', type=int, default=500000, help='bin size for Hi-C links statistics, default: %(default)s')
    parser.add_argument('--threads', type=int, default=8, help='number of threads for bam reading, default: %(default)s')
    args = parser.parse_args()

    truth_ctg_bin_dict, _ = parse_agp(args.agp_truth, args.bin_size)
    haphic_ctg_bin_dict, haphic_bin_index_dict = parse_agp(args.agp_haphic, args.bin_size)
    haphic_fastsort_ctg_bin_dict, haphic_fastsort_bin_index_dict = parse_agp(args.agp_haphic_fastsort, args.bin_size)
    yahs_ctg_bin_dict, yahs_bin_index_dict = parse_agp(args.agp_yahs, args.bin_size)
    
    truth_link_dict, haphic_link_dict, haphic_fastsort_link_dict, yahs_link_dict = parse_bam(
            args.bam, truth_ctg_bin_dict, haphic_ctg_bin_dict, haphic_fastsort_ctg_bin_dict, yahs_ctg_bin_dict, args.threads)
    
    summary_intra_dict = collections.defaultdict(dict)
    summary_inter_dict = collections.defaultdict(list)
    with open('stat.txt', 'w') as fstat:
        make_stat(fstat, truth_link_dict, 'Truth', summary_intra_dict, summary_inter_dict)
        make_stat(fstat, haphic_link_dict, 'HapHiC', summary_intra_dict, summary_inter_dict)
        make_stat(fstat, haphic_fastsort_link_dict, 'HapHiC_fastsort', summary_intra_dict, summary_inter_dict)
        make_stat(fstat, yahs_link_dict, 'YaHS', summary_intra_dict, summary_inter_dict)

    summarize(summary_intra_dict, summary_inter_dict, haphic_bin_index_dict, haphic_fastsort_bin_index_dict, yahs_bin_index_dict)


if __name__ == '__main__':
    main()

