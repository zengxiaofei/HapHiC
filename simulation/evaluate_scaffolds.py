#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-03-12 11:20


import argparse
import collections
import random

def parse_truth(truth):

    truth_dict = dict()
    ctg_info_dict = dict()
    with open(truth) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                truth_dict[ID] = []
                index = 0
            else:
                index += 1
                ctg, orient = line.split()
                truth_dict[ID].append((ctg, orient))
                ctg_info_dict[ctg] = [ID, index, orient]

    return truth_dict, ctg_info_dict


def parse_agp(agp, ctg_info_dict):
    
    scaffold_dict = collections.defaultdict(dict)
    scaffold_nctgs_dict = collections.defaultdict(int)
    with open(agp) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            if cols[4] != 'W':
                continue
            # statistics for ctgs
            ctg = cols[5]
            ctg_len = int(cols[2]) - int(cols[1]) + 1
            assert len(ctg_info_dict[ctg]) == 3
            ctg_info_dict[ctg].append(ctg_len)
            
            # statistics for scaffolds
            scaffold = cols[0]
            scaffold_nctgs_dict[scaffold] += 1
            source_scaffold = ctg_info_dict[ctg][0]
            orient = cols[8]
            if source_scaffold in scaffold_dict[scaffold]:
                scaffold_dict[scaffold][source_scaffold].append((ctg, orient))
            else:
                scaffold_dict[scaffold][source_scaffold] = [(ctg, orient)]

    return scaffold_dict, scaffold_nctgs_dict


def evaluate(scaffold_dict, scaffold_nctgs_dict, truth_dict, ctg_info_dict):
    
    def find_lis(ctg_order_list, forward=True):
        
        order_list = []
        order_ctg_dict = dict()
        order_len_dict = dict()
        
        for ctg, order in ctg_order_list:
            if forward and order < 0:
                continue
            if not forward and order > 0:
                continue
            order_list.append(order)
            assert order not in order_ctg_dict
            order_ctg_dict[order] = ctg
            order_len_dict[order] = ctg_info_dict[ctg][3]
        
        if not order_list:
            return 0, []

        dp = [0] * len(order_list)
        sequence = [None] * len(order_list)
        max_sum_idx = 0

        for i in range(len(order_list)):
            dp[i] = order_len_dict[order_list[i]]
            for j in range(i):
                if order_list[i] > order_list[j] and dp[i] < dp[j] + order_len_dict[order_list[i]]:
                    dp[i] = dp[j] + order_len_dict[order_list[i]]
                    sequence[i] = j
            if dp[i] >= dp[max_sum_idx]:
                max_sum_idx = i
        max_sum = dp[max_sum_idx]
        max_sequence = []
        while max_sum_idx is not None:
            max_sequence.append(order_ctg_dict[order_list[max_sum_idx]])
            max_sum_idx = sequence[max_sum_idx]
        max_sequence.reverse()

        return max_sum, max_sequence

    
    def ctg_stat(ctgs):
        
        length = 0
        for ctg in ctgs:
            length += ctg_info_dict[ctg][3]
        
        return len(ctgs), length
    
    # get lengths for source scaffolds
    white_list = []
    source_scaffold_len_dict = dict()
    for ID, ctg_list in truth_dict.items():
        source_scaffold_len_dict[ID] = 0
        for ctg, _ in ctg_list:
            source_scaffold_len_dict[ID] += ctg_info_dict[ctg][-1]
        if ID != 'unanchored' and len(ctg_list) == 1:
            white_list.append(ctg_list[0][0])
    
    unanchored_list = []
    newly_anchored_list = []
    scaffold_stat_dict = dict()
    for scaffold, source_scaffold_dict in scaffold_dict.items():
        # not anchored (nctgs < 2)
        if scaffold_nctgs_dict[scaffold] < 2 and list(source_scaffold_dict.values())[0][0][0] not in white_list:
            assert scaffold_nctgs_dict[scaffold] == 1
            unanchored_list.append(list(source_scaffold_dict.values())[0][0][0])
            continue
        # anchored (nctgs >= 2)
        scaffold_stat_dict[scaffold] = []
        for source_scaffold, ctg_list in source_scaffold_dict.items():
            # newly anchored
            if source_scaffold == 'unanchored':
                for ctg, _ in ctg_list:
                    newly_anchored_list.append(ctg)
                continue
            scaffold_stat_dict[scaffold].append([source_scaffold])
            # length ratio
            len_sum = 0
            for ctg, orient in ctg_list:
                len_sum += ctg_info_dict[ctg][-1]
            source_scaffold = ctg_info_dict[ctg][0]
            source_scaffold_len = source_scaffold_len_dict[source_scaffold]
            len_ratio = len_sum / source_scaffold_len
            scaffold_stat_dict[scaffold][-1].extend([len_ratio, len_sum])
        
        # random shuffled
        random.seed(12345)
        random.shuffle(scaffold_stat_dict[scaffold])
        # sort by len_ratio and len_sum
        scaffold_stat_dict[scaffold].sort(key=lambda x: x[1:], reverse=True)
       

    # get dominant source for each scaffold
    dominant_scaffold_dict = dict()
    dominant_source_dict = dict()
    for scaffold, scaffold_info in scaffold_stat_dict.items():
        if not scaffold_info:
            continue
        source_scaffold, len_ratio = scaffold_info[0][:2]
        if source_scaffold not in dominant_scaffold_dict or len_ratio > dominant_scaffold_dict[source_scaffold][1]:
            dominant_scaffold_dict[source_scaffold] = (scaffold, len_ratio)

    for source_scaffold, (scaffold, len_ratio) in dominant_scaffold_dict.items():
        dominant_source_dict[scaffold] = source_scaffold
    
    # statistics for inversions, relocations, and translocations
    translocation_list = []
    relocation_list = []
    inversion_list = []
    inversion_and_relocation_list = []
    syntenic_list = []

    for scaffold, source_scaffold_dict in scaffold_dict.items():
        # unanchored
        if scaffold_nctgs_dict[scaffold] < 2 and list(source_scaffold_dict.values())[0][0][0] not in white_list:
            continue
        for source_scaffold, ctg_list in source_scaffold_dict.items():
            # translocations
            if source_scaffold != 'unanchored' and (scaffold not in dominant_source_dict or source_scaffold != dominant_source_dict[scaffold]):
                translocation_list.extend([ctg for ctg, _ in ctg_list])
            elif source_scaffold != 'unanchored':
                # (1) find LIS
                ctg_order_list = []
                for ctg, orient_in_scaf in ctg_list:
                    index, orient_in_truth = ctg_info_dict[ctg][1:3]
                    if orient_in_scaf == orient_in_truth:
                        ctg_order_list.append((ctg, index))
                    else:
                        ctg_order_list.append((ctg, -index))
                # print('')
                # print(scaffold, source_scaffold)
                # print('\tctg_order_list:', ctg_order_list)
                
                max_sum_f, lis_ctg_list_f = find_lis(ctg_order_list, forward=True)
                max_sum_r, lis_ctg_list_r = find_lis(ctg_order_list, forward=False)
                
                # lis_ctg_list: LIS in both orientations
                # the contigs in this list are defined as syntenic contigs
                if max_sum_f >= max_sum_r:
                    lis_ctg_list = lis_ctg_list_f
                    lis_order = 1
                else:
                    lis_ctg_list = lis_ctg_list_r
                    lis_order = -1
                # print('\tLIS:', lis_ctg_list)
                syntenic_list.extend(lis_ctg_list)
                # (2) merge contigs that are not in LIS
                last_ctg = None
                last_order = 0
                merged_ctg_order_list = []
                for ctg, order in ctg_order_list:
                    if ctg in lis_ctg_list:
                        merged_ctg_order_list.append((ctg, order))
                        last_ctg = None
                        last_order = 0
                    else:
                        if order * last_order > 0 and order == last_order + 1:
                            merged_ctg_order_list[-1].append((ctg, order))
                        else:
                            merged_ctg_order_list.append([(ctg, order)])
                        last_ctg = ctg
                        last_order = order
                
                # print('\tmerged_ctg_order_list:', merged_ctg_order_list)
                
                # (3) identify relocations and inverse the segments that have a different order with LIS
                new_ctg_order_list = []
                all_inversion_list = []
                for i, segment in enumerate(merged_ctg_order_list):
                    # contigs that are not in LIS
                    if isinstance(segment, list):
                        # segments that have the same order with LIS are defined as relocations
                        if segment[0][1] * lis_order > 0:
                            relocation_list.extend([ctg for ctg, _ in segment])
                            new_ctg_order_list.extend(segment)
                        # different order
                        else:
                            all_inversion_list.extend(segment)
                            # inverse the segment that have a different order with LIS
                            new_ctg_order_list.extend([(ctg, -order) for ctg, order in segment[::-1]])
                    else:
                        new_ctg_order_list.append(segment)
                
                # print('\tnew_ctg_order_list:', new_ctg_order_list)

                # (4) find LIS again, and compare the new one with the original one

                if lis_order == 1:
                    _, new_lis_ctg_list = find_lis(new_ctg_order_list, forward=True)
                else:
                    assert lis_order == -1
                    _, new_lis_ctg_list = find_lis(new_ctg_order_list, forward=False)
                
                # print('\tnew_lis_ctg_list:', new_lis_ctg_list)
                
                # identify inversions
                for ctg in new_lis_ctg_list:
                    if ctg not in lis_ctg_list and ctg not in relocation_list:
                        inversion_list.append(ctg)

                # re-identify relocation
                re_identified_relocation_set = set(lis_ctg_list) - set(new_lis_ctg_list)
                relocation_list.extend(list(re_identified_relocation_set))
                syntenic_list = [ctg for ctg in syntenic_list if ctg not in re_identified_relocation_set]

                # identify contigs that need both inversion and translocation
                for ctg, order in all_inversion_list:
                    if ctg not in inversion_list:
                        inversion_and_relocation_list.append(ctg)

    # print('\n###### scaffold composition ######')
    # print('len sum and ratio of each source:', scaffold_stat_dict)
    # print('dominant source for each scaffold', dominant_source_dict)
    # 
    # print('\n###### results ######')
    # print('syntenic_list:', syntenic_list)
    # print('unanchored_list:', unanchored_list)
    # print('newly_anchored_list:', newly_anchored_list)
    # print('translocation_list:', translocation_list)
    # print('relocation_list:', relocation_list)
    # print('inversion_list:', inversion_list)
    # print('inversion_and_relocation_list:', inversion_and_relocation_list)
    
    total_ctg_set = set(ctg_info_dict.keys())
    total_ctg_num, total_ctg_len = ctg_stat(total_ctg_set)
    syntenic_ctg_num, syntenic_ctg_len = ctg_stat(syntenic_list)
    unanchored_ctg_num, unanchored_ctg_len = ctg_stat(unanchored_list)
    newly_anchored_ctg_num, newly_anchored_ctg_len = ctg_stat(newly_anchored_list)
    translocation_ctg_num, translocation_ctg_len = ctg_stat(translocation_list)
    relocation_ctg_num, relocation_ctg_len = ctg_stat(relocation_list)
    inversion_ctg_num, inversion_ctg_len = ctg_stat(inversion_list)
    inversion_and_relocation_ctg_num, inversion_and_relocation_ctg_len = ctg_stat(inversion_and_relocation_list)

    print('\n###### summary ######')
    print('Number of scaffolds (at least two contigs):\n {} scaffolds'.format(len([scaffold for scaffold, nctgs in scaffold_nctgs_dict.items() if nctgs > 1 or (nctgs == 1 and list(scaffold_dict[scaffold].values())[0][0][0] in white_list)])))
    print('Total contigs:\n {} / {} bp'.format(total_ctg_num, total_ctg_len))
    print('Syntenic contigs:\n {} / {} bp / {} %'.format(syntenic_ctg_num, syntenic_ctg_len, syntenic_ctg_len/total_ctg_len*100))
    print('Unanchored contigs:\n {} / {} bp / {} %'.format(unanchored_ctg_num, unanchored_ctg_len, unanchored_ctg_len/total_ctg_len*100))
    print('Newly_anchored contigs:\n {} / {} bp / {} %'.format(newly_anchored_ctg_num, newly_anchored_ctg_len, newly_anchored_ctg_len/total_ctg_len*100))
    print('Translocation contigs:\n {} / {} bp / {} %'.format(translocation_ctg_num, translocation_ctg_len, translocation_ctg_len/total_ctg_len*100))
    print('Relocation contigs:\n {} / {} bp / {} %'.format(relocation_ctg_num, relocation_ctg_len, relocation_ctg_len/total_ctg_len*100))
    print('Inversion contigs:\n {} / {} bp / {} %'.format(inversion_ctg_num, inversion_ctg_len, inversion_ctg_len/total_ctg_len*100))
    print('Inversion and relocation contigs:\n {} / {} bp / {} %'.format(inversion_and_relocation_ctg_num, inversion_and_relocation_ctg_len, inversion_and_relocation_ctg_len/total_ctg_len*100))

    # check results
    count_dict = collections.defaultdict(list)
    
    for n, list_ in enumerate([syntenic_list, unanchored_list, newly_anchored_list, translocation_list, relocation_list, inversion_list, inversion_and_relocation_list]):
        for ctg in list_:
            count_dict[ctg].append(n)

    for ctg, index_list in count_dict.items():
        if len(index_list) > 1:
            print(ctg, index_list)
    
    assert unanchored_ctg_num + newly_anchored_ctg_num + translocation_ctg_num + relocation_ctg_num + inversion_ctg_num + inversion_and_relocation_ctg_num + syntenic_ctg_num == total_ctg_num
    assert set(unanchored_list) | set(newly_anchored_list) | set(translocation_list) | set(relocation_list) | set(inversion_list) | set(inversion_and_relocation_list) | set(syntenic_list) == total_ctg_set

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('truth', help='ordering and orientation information for ground truth')
    parser.add_argument('agp', help='agp file for scaffolds')
    args = parser.parse_args()

    truth_dict, ctg_info_dict = parse_truth(args.truth)
    scaffold_dict, scaffold_nctgs_dict = parse_agp(args.agp, ctg_info_dict)
    evaluate(scaffold_dict, scaffold_nctgs_dict, truth_dict, ctg_info_dict)

if __name__ == '__main__':
    main()

