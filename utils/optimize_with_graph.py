#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created at: 2024-10-25 10:05


import sys
import argparse
import gzip
import collections
import networkx as nx
from itertools import combinations, product
from pywfa import WavefrontAligner



def parse_assembly(assembly):
    
    ctg_list = []
    group_list = []
    ctg_dict = dict()
    unanchored_ctg_semis = set()
    frag_dict = collections.defaultdict(list)
    group = 0
    with open(assembly) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if line.startswith('>'):
                ID = cols[0][1:]
                num = cols[1]
                length = int(cols[2])
                if ':::' not in ID:
                    ctg_list.append((ID, ID, 1, length))
                else:
                    ID_split = ID.split(':::')
                    raw_ID, frag_ID = ID_split[:2]
                    frag_dict[raw_ID].append(ID)
                    if frag_ID == 'fragment_1':
                        start = 0
                    ctg_list.append([ID, raw_ID, start + 1, length + start])
                    start += length
            else:
                group_list.append([])
                for order, num in enumerate(cols):
                    num = int(num)
                    if num < 0:
                        ID, raw_ID, start, end = ctg_list[-num - 1]
                        ctg_dict[ID] = [group, order, -num, '-', raw_ID, start, end]
                        group_list[-1].append([ID, '-'])
                    else:
                        ID, raw_ID, start, end = ctg_list[num - 1]
                        ctg_dict[ID] = [group, order, num, '+', raw_ID, start, end]
                        group_list[-1].append([ID, '+'])
                group += 1
                
                # add unanchored contigs
                if order == 0:
                    unanchored_ctg_semis.add((ID, 'Head'))
                    unanchored_ctg_semis.add((ID, 'Tail'))
    
    return ctg_dict, group_list, frag_dict, unanchored_ctg_semis


def parse_collapsed_regions(collapsed_regions, ctg_dict, group_list):
    
    ctg_region_dict = collections.defaultdict(list)
    collapsed_region_list = []
    order = 0
    with open(collapsed_regions) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            region, start_ctg, end_ctg, after_ctg, before_ctg, orientation = cols[:-1]
            
            # check format
            if start_ctg != end_ctg:
                assert ctg_dict[start_ctg][0] == ctg_dict[end_ctg][0]
                assert ctg_dict[start_ctg][1] < ctg_dict[end_ctg][1]
            assert after_ctg != 'NA' or before_ctg != 'NA'
            if after_ctg != 'NA' and before_ctg != 'NA':
                assert ctg_dict[after_ctg][0] == ctg_dict[before_ctg][0]
                assert ctg_dict[after_ctg][1] + 1 == ctg_dict[before_ctg][1]
            
            # record regions in a new list and find regions with the number of collapsed haplotypes >= 3
            from_group = ctg_dict[start_ctg][0]
            to_group = ctg_dict[after_ctg][0] if after_ctg != 'NA' else ctg_dict[before_ctg][0]
            order_start = ctg_dict[start_ctg][1]
            order_end = ctg_dict[end_ctg][1]
            collapsed_ctg_list = group_list[from_group][order_start:order_end + 1]
            
            # print('\n' + line, end='')
            # print('{} {}: {}'.format(from_group, to_group, ','.join([ctg_dict[ctg][4] for ctg, orientation in collapsed_ctg_list])), file=sys.stderr)
            
            for ctg, _ in collapsed_ctg_list:
                ctg_region_dict[ctg].append((region, order))
            collapsed_region_list.append([region, collapsed_ctg_list, to_group, after_ctg, before_ctg, orientation])
            order += 1

    return ctg_region_dict, collapsed_region_list


def parse_liftover(liftover, ctg_dict, frag_dict):

    # an old ID can be mapped to multiple new IDs
    old_new_dict = collections.defaultdict(list)
    # a new ID can only be mapped to an old ID
    new_old_dict = dict()
    with open(liftover) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            if cols[4] == 'W':
                offset = int(cols[6])
                old_ID = cols[5]
                new_ID = cols[0]
                if new_ID in frag_dict:
                    for ID in frag_dict[new_ID]:
                        ctg_dict[ID][-2] += offset - 1
                        ctg_dict[ID][-1] += offset - 1
                else:
                    ctg_dict[new_ID][-2] += offset - 1
                    ctg_dict[new_ID][-1] += offset - 1

                old_new_dict[old_ID].append(new_ID)
                new_old_dict[new_ID] = old_ID

    for old_ID, new_ID_list in old_new_dict.items():
        if len(new_ID_list) > 1:
            new_ID_list.sort(key=lambda x: ctg_dict[x][-2])
            print(old_ID, new_ID_list)

    return old_new_dict, new_old_dict


def parse_gfa(gfa, frag_dict, old_new_dict):

    def select_ctg(a_dict, key, orient, from_to):
        # select the first or the last contig according to the orientation
        if from_to == 'from':
            return a_dict[key][-1] if orient == '+' else a_dict[key][0]
        else:
            assert from_to == 'to'
            return a_dict[key][0] if orient == '+' else a_dict[key][-1]

    if gfa.endswith('gz'):
        fopen = gzip.open
    else:
        fopen = open
    
    frag_rd_dict = dict()
    ctg_len_dict = dict()
    ctg_ovl_dict = dict()

    G = nx.Graph()
    with fopen(gfa, 'rt') as f:
        for line in f:
            if line.startswith('S'):
                cols = line.split()
                
                if old_new_dict:
                    ctgs = old_new_dict[cols[1]]
                else:
                    ctgs = [cols[1]]
                
                rd = int(cols[4].split(':')[-1]) + 1
                ctg_len_dict[cols[1]] = int(cols[3].split(':')[-1])

                for ctg in ctgs:
                    if ctg in frag_dict:
                        for frag in frag_dict[ctg]:
                            G.add_edge((frag, 'Head'), (frag, 'Tail'))
                            frag_rd_dict[frag] = rd
                    else:
                        G.add_edge((ctg, 'Head'), (ctg, 'Tail'))
                        frag_rd_dict[ctg] = rd
            
            elif line.startswith('L'):
                cols = line.split()
                from_ctg, from_orient, to_ctg, to_orient = cols[1:5]
                
                if old_new_dict:
                    from_ctg = select_ctg(old_new_dict, from_ctg, from_orient, from_to='from')
                    to_ctg = select_ctg(old_new_dict, to_ctg, to_orient, from_to='to')
                
                if from_ctg in frag_dict:
                    from_ctg = select_ctg(frag_dict, from_ctg, from_orient, from_to='from')
                
                if to_ctg in frag_dict:
                    to_ctg = select_ctg(frag_dict, to_ctg, to_orient, from_to='to')
                
                if from_orient == '+':
                    from_semi = 'Tail'
                else:
                    from_semi = 'Head'
                if to_orient == '+':
                    to_semi = 'Head'
                else:
                    to_semi = 'Tail'
                
                ctg_ovl_dict[(from_ctg, from_orient, to_ctg, to_orient)] = (int(cols[5][:-1]), int(cols[6].split(':')[-1]))

                if not G.has_edge((from_ctg, from_semi), (to_ctg, to_semi)):
                    G.add_edge((from_ctg, from_semi), (to_ctg, to_semi))

    return G, frag_rd_dict, ctg_len_dict, ctg_ovl_dict


def determine_haplotype_depths(G, frag_rd_dict, expected_nhaps):
    
    def get_annother_end(node):
        if node[1] == 'Head':
            return node[0], 'Tail'
        else:
            return node[0], 'Head'

    def filter_nodes(neighbors):
        if len(neighbors) == 2 and neighbors[0][0] != neighbors[1][0]:
            if G.degree(neighbors[0]) == 2 and G.degree(neighbors[1]) == 2:
                annother_end_0 = get_annother_end(neighbors[0])
                annother_end_1 = get_annother_end(neighbors[1])
                if G.degree(annother_end_0) == 2 and G.degree(annother_end_1) == 2:
                    annother_end_0_neighbor = [neighbor for neighbor in G.neighbors(annother_end_0) if neighbor[0] != annother_end_0[0]]
                    annother_end_1_neighbor = [neighbor for neighbor in G.neighbors(annother_end_1) if neighbor[0] != annother_end_1[0]]
                    if annother_end_0_neighbor == annother_end_1_neighbor:
                        return {annother_end_0[0], annother_end_1[0]}

    # find collapsed and heterozygous frags

    collapsed_depths = []
    heterozygous_depths = []

    for frag, rd in frag_rd_dict.items():
        
        # find putative collapsed contigs: >-<
        if G.degree((frag, 'Head')) == 3 and G.degree((frag, 'Tail')) == 3:
            # find putative heterozygous contigs
            head_frags = filter_nodes([neighbor for neighbor in G.neighbors((frag, 'Head')) if neighbor != (frag, 'Tail')])
            tail_frags = filter_nodes([neighbor for neighbor in G.neighbors((frag, 'Tail')) if neighbor != (frag, 'Head')])
        
            # add both types of contigs
            if head_frags and tail_frags:
                collapsed_depths.append(rd)
                for heter_frag in head_frags | tail_frags:
                    heterozygous_depths.append(frag_rd_dict[heter_frag])

    # determine haplotype depths
    heter_depth_counts = collections.Counter(heterozygous_depths)
    # 1X depth is the most frequent depth of the heterozygous contigs
    heter_m = max(heter_depth_counts, key=lambda x: heter_depth_counts[x])
    
    collapsed_depth_counts = collections.Counter(collapsed_depths)
    # 2X depth is the most frequent depth of the collapsed contigs
    collapsed_m = max(collapsed_depth_counts, key=lambda x: collapsed_depth_counts[x])
    
    # calculate the thresholds for collapsed regions
    hap_thresholds = [0, (heter_m + collapsed_m) / 2]

    for i in range(expected_nhaps - 2):
        hap_thresholds.append(collapsed_m + (collapsed_m - hap_thresholds[1]) * (i + 1))
    
    print(heter_m, collapsed_m, file=sys.stderr)
    print(hap_thresholds, file=sys.stderr)
    
    return hap_thresholds


def copy_collapsed_regions(G, collapsed_region_list, ctg_dict, group_list, unanchored_ctg_semis, frag_rd_dict, hap_thresholds, ctg_region_dict):
    
    def guess_haplotype_depth(rd):
        
        for i in range(len(hap_thresholds) - 1):
            if hap_thresholds[i] <= rd < hap_thresholds[i+1]:
                return i + 1
        # the haplotype depth will not be larger than the expected number of haplotypes
        return i + 1

    semi_set = {'Head', 'Tail'}
    original_position_ctgs = set()
    moved_ctgs = set()
    copied_collapsed_regions_dict = collections.defaultdict(list)
    for region, collapsed_ctg_list, to_group, after_ctg, before_ctg, orientation in collapsed_region_list:
        
        print('\n{}: {}'.format(region, ','.join([ctg_dict[ctg][4] for ctg, _ in collapsed_ctg_list])))
        
        if len(collapsed_ctg_list) == 1:
            print('only one contig in this region, it should be a collapsed contig, skip')
            collapsed_ctgs = {collapsed_ctg_list[0][0]}
            allelic_ctgs = set()
            singleton_allelic_ctgs = set()
        else:
            # construct the subgraph of the group
            # check if all contigs in the region are in the same group
            group = ctg_dict[collapsed_ctg_list[0][0]][0]
            for ctg, _ in collapsed_ctg_list:
                assert ctg_dict[ctg][0] == group
            # add contigs in the group to the subgraph
            subnodes = set()
            for ctg, _ in group_list[group]:
                subnodes.add((ctg, 'Head'))
                subnodes.add((ctg, 'Tail'))
            # add unanchored contigs to the subgraph
            sub_G = G.subgraph(subnodes | unanchored_ctg_semis)

            # find potential heterozygous contigs
            all_ctg_set = set()
            heter_ctg_set = set()
            for ctg, orient in collapsed_ctg_list:
                all_ctg_set.add(ctg)
                if sub_G.degree((ctg, 'Head')) <= 2 and sub_G.degree((ctg, 'Tail')) <= 2 and sub_G.degree((ctg, 'Head')) + sub_G.degree((ctg, 'Tail')) > 2:
                    heter_ctg_set.add(ctg)
        
            # find corresponding heterozygous contigs using heter_G
            heter_G = nx.Graph()
            for (ctg1, semi1), (ctg2, semi2) in combinations(product(heter_ctg_set, semi_set), 2):
                if ctg1 == ctg2:
                    continue
                # check if there is a path between ctg1 and ctg2 in sub_G, and the path length is 2
                if nx.has_path(sub_G, (ctg1, semi1), (ctg2, semi2)) and nx.shortest_path_length(sub_G, (ctg1, semi1), (ctg2, semi2)) == 2:
                    shortest_paths = nx.shortest_simple_paths(sub_G, (ctg1, semi1), (ctg2, semi2))
                    # traverse the path and find the first node that is not ctg1 or ctg2
                    for shortest_path in shortest_paths:
                        if shortest_path[1][0] not in {ctg1, ctg2}:
                            heter_G.add_edge(ctg1, ctg2)
                            break
        
            allelic_ctgs = set()
            selected_allelic_ctgs = set()
            for connected_allelic_ctgs in nx.connected_components(heter_G):
                allelic_ctgs |= connected_allelic_ctgs
                print('connected_allelic_ctgs:', [(ctg, guess_haplotype_depth(frag_rd_dict[ctg]), len(ctg_region_dict[ctg]) + 1) for ctg in connected_allelic_ctgs])
                # select one of the connected allelic contigs
                selected_ctg = sorted(list(connected_allelic_ctgs - moved_ctgs))[0]
                selected_allelic_ctgs.add(selected_ctg)
                moved_ctgs.add(selected_ctg)
            
            singleton_allelic_ctgs = set()
            collapsed_ctgs = set() 
            for ctg in heter_ctg_set - allelic_ctgs:
                # find potential singleton allelic contigs in unconnected contigs
                if guess_haplotype_depth(frag_rd_dict[ctg]) == 1:
                    singleton_allelic_ctgs.add(ctg)
                # otherwise, this contig should be a collapsed contig and be rescued to collapsed_ctgs
                else:
                    collapsed_ctgs.add(ctg)

            # for all contigs that have a degree >= 3 in at least one end
            for ctg in all_ctg_set - heter_ctg_set:
                #  if the depth is too low, this contig should be considered as a singleton allelic contig
                if guess_haplotype_depth(frag_rd_dict[ctg]) == 1:
                    singleton_allelic_ctgs.add(ctg)
                # otherwise, a collapsed contig
                else:
                    collapsed_ctgs.add(ctg)

            print('singleton_ctgs:', [(ctg, guess_haplotype_depth(frag_rd_dict[ctg]), len(ctg_region_dict[ctg]) + 1) for ctg in singleton_allelic_ctgs])
            for ctg in collapsed_ctgs:
                print('collapsed_ctgs:', ctg, guess_haplotype_depth(frag_rd_dict[ctg]), len(ctg_region_dict[ctg]) + 1)

        print(len(allelic_ctgs), len(singleton_allelic_ctgs), len(collapsed_ctgs))

        # get contigs that will be copied in this region
        copied_ctg_list = []
        # collapsed_ctg_list is already correctly ordered
        for ctg, orient in collapsed_ctg_list:
            if ctg in singleton_allelic_ctgs:
                # singleton allelic contigs will not be copied, stay in the original position
                original_position_ctgs.add(ctg)
                continue
            if ctg in collapsed_ctgs:
                # all collapsed contigs will be copied without checking the depth
                copied_ctg_list.append((ctg, orient))
                original_position_ctgs.add(ctg)
            else:
                assert ctg in allelic_ctgs
                # one of the connected allelic contigs will be randomly selected for each region
                if ctg in selected_allelic_ctgs:
                    copied_ctg_list.append((ctg, orient))
        
        print('copied_ctg_list:', copied_ctg_list)
        copied_collapsed_regions_dict[to_group].append([copied_ctg_list, orientation, after_ctg, before_ctg])
    
    return copied_collapsed_regions_dict, moved_ctgs


def get_contig_seqs(contigs):
    
    ctg_seq_dict = dict()
    with open(contigs) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                ctg_seq_dict[ID] = []
            else:
                ctg_seq_dict[ID].append(line.strip())

    for ID, seq_list in ctg_seq_dict.items():
        ctg_seq_dict[ID] = ''.join(seq_list)

    return ctg_seq_dict


def list_output_contigs(copied_collapsed_regions_dict, moved_ctgs, ctg_dict, group_list, new_old_dict, frag_rd_dict, ctg_seq_dict, ctg_len_dict, ctg_ovl_dict, unanchored_ctg_semis,  G, args):
    
    def reverse(orient):
        if orient == '+':
            return '-'
        else:
            return '+'
        
    def revcomp(seq):
        return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]

    def is_terminus(ctg, ctg_pos):
        if ctg_pos == 1:
            return True
        if ':::' in ctg:
            ctg = ctg.split(':::')[0]
        old_ID = new_old_dict[ctg] if new_old_dict else ctg
        if ctg_pos == ctg_len_dict[old_ID]:
            return True

    def get_distance(path):
        # direct connection
        path_len = len(path)
        if path_len == 0:
            return 0
        elif path_len % 2 == 0 and all([path[i*2][0] == path[i*2+1][0] for i in range(path_len//2)]):
            return path_len / 2
        return -1
    
    def optimize():
        for n in range(args.nrounds):
            changes = 0
            for group, ctg_list in enumerate(group_list):
                sub_G = sub_G_list[group]
                for i in range(len(ctg_list) - 1):
                    
                    ctg_i, orient_i = ctg_list[i][:2]
                    ctg_j, orient_j = ctg_list[i+1][:2]
                    
                    ctg_start_i, ctg_end_i = ctg_dict[ctg_i][-2:]
                    ctg_start_j, ctg_end_j = ctg_dict[ctg_j][-2:]   
                    
                    i_head = (ctg_i, 'Head')
                    i_tail = (ctg_i, 'Tail')
                    j_head = (ctg_j, 'Head')
                    j_tail = (ctg_j, 'Tail')

                    print(ctg_i, orient_i, ctg_start_i, ctg_end_i, ctg_j, orient_j, ctg_start_j, ctg_end_j)

                    # + +
                    if is_terminus(ctg_i, ctg_end_i) and is_terminus(ctg_j, ctg_start_j) and nx.has_path(sub_G, i_tail, j_head) and 0 <= get_distance(nx.shortest_path(sub_G, i_tail, j_head)[1:-1]) <= args.distance:
                        ctg_list[i][1] = '+'
                        ctg_list[i+1][1] = '+'
                    # + -
                    elif is_terminus(ctg_i, ctg_end_i) and is_terminus(ctg_j, ctg_end_j) and nx.has_path(sub_G, i_tail, j_tail) and 0 <= get_distance(nx.shortest_path(sub_G, i_tail, j_tail)[1:-1]) <= args.distance:
                        ctg_list[i][1] = '+'
                        ctg_list[i+1][1] = '-'
                    # - +
                    elif is_terminus(ctg_i, ctg_start_i) and is_terminus(ctg_j, ctg_start_j) and nx.has_path(sub_G, i_head, j_head) and 0 <= get_distance(nx.shortest_path(sub_G, i_head, j_head)[1:-1]) <= args.distance:
                        ctg_list[i][1] = '-'
                        ctg_list[i+1][1] = '+'
                    # - - 
                    elif is_terminus(ctg_i, ctg_start_i) and is_terminus(ctg_j, ctg_end_j) and nx.has_path(sub_G, i_head, j_tail) and 0 <= get_distance(nx.shortest_path(sub_G, i_head, j_tail)[1:-1]) <= args.distance:
                        ctg_list[i][1] = '-'
                        ctg_list[i+1][1] = '-'
                    else:
                        # do nothing
                        print('Round {}: {} {} not connected on assembly graph or too far, do nothing'.format(n+1, ctg_i, ctg_j))
                        pass

                    if orient_i != ctg_list[i][1] or orient_j != ctg_list[i+1][1]:
                        changes += 1
                        print('Round {0}: {1} {3} {2} {4} -> {1} {5} {2} {6}'.format(n+1, ctg_i, ctg_j, orient_i, orient_j, ctg_list[i][1], ctg_list[i+1][1]))
                    else:
                        print('Round {0}: unchanged'.format(n+1))
            
            # if there are no changes, break the loop
            if changes == 0:
                break

    def rescue():
        
        def update_rescued_group_list(sub_G, i_edge, j_edge, ctg_list, i, rescued_ctg_set):
            if i == 0:
                rescued_group_list[-1].append(ctg_list[i])
            path = []
            for x, path in enumerate(nx.all_simple_paths(sub_G, i_edge, j_edge, cutoff=20)):
                if x == 1:
                    rescued_group_list[-1].append(ctg_list[i+1])
                    return
            rescued_ctgs = []
            if path and len(path) <= 10 and all([G.degree(e) == 2 and e in unanchored_ctg_semis for e in path[1:-1]]):
                for n in range(len(path[1:-1])//2):
                    rescued_ctg, semi = path[1:-1][2*n]
                    rescued_ctgs.append(rescued_ctg)
                    if semi == 'Head':
                        if rescued_ctg in rescued_ctg_set:
                            rescued_group_list[-1].append((rescued_ctg, '+', 'copy'))
                        else:
                            rescued_group_list[-1].append((rescued_ctg, '+'))
                    else:
                        if rescued_ctg in rescued_ctg_set:
                            rescued_group_list[-1].append((rescued_ctg, '-', 'copy'))
                        else:
                            rescued_group_list[-1].append((rescued_ctg, '-'))
            rescued_group_list[-1].append(ctg_list[i+1])
            rescued_ctg_set |= set(rescued_ctgs)
            return rescued_ctgs

        rescued_ctg_set = set()
        rescued_group_list = []
        for group, ctg_list in enumerate(group_list):
            rescued_group_list.append([])
            if len(ctg_list) == 1 and ctg_list[0][0] not in rescued_ctg_set:
                rescued_group_list[-1].append(ctg_list[0])
                continue
            
            # filter subgraphs 
            sub_G = sub_G_list[group].copy()
            for ctg in {node[0] for node in sub_G.nodes() if node in unanchored_ctg_semis}:
                if frag_rd_dict[ctg] == 1:
                    sub_G.remove_node((ctg, 'Head'))
                    sub_G.remove_node((ctg, 'Tail'))
            
            for i in range(len(ctg_list) - 1):
                ctg_i, orient_i = ctg_list[i][:2]
                ctg_j, orient_j = ctg_list[i+1][:2]

                ctg_start_i, ctg_end_i = ctg_dict[ctg_i][-2:]
                ctg_start_j, ctg_end_j = ctg_dict[ctg_j][-2:]   
                
                i_head = (ctg_i, 'Head')
                i_tail = (ctg_i, 'Tail')
                j_head = (ctg_j, 'Head')
                j_tail = (ctg_j, 'Tail')

                # + +
                if orient_i == '+' and orient_j == '+' and is_terminus(ctg_i, ctg_end_i) and is_terminus(ctg_j, ctg_start_j) and nx.has_path(sub_G, i_tail, j_head) and not sub_G.has_edge(i_tail, j_head):
                    rescued_ctgs = update_rescued_group_list(sub_G, i_tail, j_head, ctg_list, i, rescued_ctg_set)
                # + -
                elif orient_i == '+' and orient_j == '-' and is_terminus(ctg_i, ctg_end_i) and is_terminus(ctg_j, ctg_end_j) and nx.has_path(sub_G, i_tail, j_tail) and not sub_G.has_edge(i_tail, j_tail):
                    rescued_ctgs = update_rescued_group_list(sub_G, i_tail, j_tail, ctg_list, i, rescued_ctg_set)
                # - +
                elif orient_i == '-' and orient_j == '+' and is_terminus(ctg_i, ctg_start_i) and is_terminus(ctg_j, ctg_start_j) and nx.has_path(sub_G, i_head, j_head) and not sub_G.has_edge(i_head, j_head):
                    rescued_ctgs = update_rescued_group_list(sub_G, i_head, j_head, ctg_list, i, rescued_ctg_set)
                # - - 
                elif orient_i == '-' and orient_j == '-' and is_terminus(ctg_i, ctg_start_i) and is_terminus(ctg_j, ctg_end_j) and nx.has_path(sub_G, i_head, j_tail) and not sub_G.has_edge(i_head, j_tail):
                    rescued_ctgs = update_rescued_group_list(sub_G, i_head, j_tail, ctg_list, i, rescued_ctg_set)
                else:
                    rescued_ctgs = None
                    if i == 0:
                        rescued_group_list[-1].extend([ctg_list[i], ctg_list[i+1]])
                    else:
                        rescued_group_list[-1].append(ctg_list[i+1])
                
                print('Rescue: Group{} {} {} {} {} {}'.format(group+1, ctg_i, orient_i, ctg_j, orient_j, rescued_ctgs)) 
        
        return rescued_group_list
    
    def close_gaps():
        
        def update_gap_closed_group_list():
            if i == 0:
                gap_closed_group_list[-1].append([ctg_list[i], ctg_list[i+1]])
            else:
                gap_closed_group_list[-1][-1].append(ctg_list[i+1])

        gap_closed_group_list = []
        for ctg_list in group_list:
            gap_closed_group_list.append([])
            if len(ctg_list) == 1:
                gap_closed_group_list[-1].append([ctg_list[0]])
                continue
            for i in range(len(ctg_list) - 1):
                ctg_i, orient_i = ctg_list[i][:2]
                ctg_j, orient_j = ctg_list[i+1][:2]
                
                ctg_start_i, ctg_end_i = ctg_dict[ctg_i][-2:]
                ctg_start_j, ctg_end_j = ctg_dict[ctg_j][-2:]   
                
                i_head = (ctg_i, 'Head')
                i_tail = (ctg_i, 'Tail')
                j_head = (ctg_j, 'Head')
                j_tail = (ctg_j, 'Tail')

                if orient_i == '+' and orient_j == '+' and is_terminus(ctg_i, ctg_end_i) and is_terminus(ctg_j, ctg_start_j) and G.has_edge(i_tail, j_head):
                    update_gap_closed_group_list()
                elif orient_i == '+' and orient_j == '-' and is_terminus(ctg_i, ctg_end_i) and is_terminus(ctg_j, ctg_end_j) and G.has_edge(i_tail, j_tail):
                    update_gap_closed_group_list()
                elif orient_i == '-' and orient_j == '+' and is_terminus(ctg_i, ctg_start_i) and is_terminus(ctg_j, ctg_start_j) and G.has_edge(i_head, j_head):
                    update_gap_closed_group_list()
                elif orient_i == '-' and orient_j == '-' and is_terminus(ctg_i, ctg_start_i) and is_terminus(ctg_j, ctg_end_j) and G.has_edge(i_head, j_tail):
                    update_gap_closed_group_list()
                else:
                    if i == 0:
                        gap_closed_group_list[-1].extend([[ctg_list[i]], [ctg_list[i+1]]])
                    else:
                        gap_closed_group_list[-1].append([ctg_list[i+1]])
        
        region_n, region_len = 0, 0
        for group, ctg_list in enumerate(gap_closed_group_list):
            print('Group {} All:'.format(group), ','.join([','.join([ctg_info[0] for ctg_info in ctg_info_list]) for ctg_info_list in ctg_list]))
            print('Group {} Gap closed:'.format(group), ' | '.join([','.join([ctg_info[0] for ctg_info in ctg_info_list]) for ctg_info_list in ctg_list if len(ctg_info_list) > 1]))
            for ctg_info_list in ctg_list:
                if len(ctg_info_list) > 1:
                    region_n += 1
                    for ctg_info in ctg_info_list:
                        start, end = ctg_dict[ctg_info[0]][-2:]
                        print(ctg_info[0], start, end)
                        region_len += end - start + 1
            print('Total gap closed regions: {}\nTotal length of gap closed regions: {}'.format(region_n, region_len))

        return gap_closed_group_list
    
    def output_files(group_list):
        # output final assembly
        counter = dict()
        with open('after_copy.agp', 'wt') as f_agp, open('after_copy_ctg.fa', 'wt') as f_ctg, open('after_copy_chr.fa', 'wt') as f_chr:
            for group, ctg_list in enumerate(group_list):
                m = 0
                start = 0
                end = 0
                scaffold_seq_list = []
                for n, ctg_info in enumerate(ctg_list):
                    ctg, orient = ctg_info[:2]
                    raw_ID = ctg_dict[ctg][4]
                    ctg_start, ctg_end = ctg_dict[ctg][-2:]
                    if ctg_info[-1] == 'copy' and ctg not in counter:
                        counter[ctg] = 1
                    elif ctg_info[-1] == 'copy':
                        counter[ctg] += 1
                    old_ID = new_old_dict[raw_ID] if new_old_dict else raw_ID
                    if ctg_info[-1] == 'copy':
                        output_ID = '{}_{}_{}_copy{}'.format(old_ID, ctg_start, ctg_end, counter[ctg])
                    else:
                        output_ID = '{}_{}_{}'.format(old_ID, ctg_start, ctg_end)
                    
                    # output AGP file
                    if n > 0:
                        m += 1
                        start = end + 1
                        end = start + 99
                        f_agp.write('scaffold{}\t{}\t{}\t{}\tU\t100\tscaffold\tyes\tproximity_ligation\n'.format(group+1, start, end, m))
                    m += 1
                    start = end + 1
                    end = start + ctg_end - ctg_start
                    f_agp.write('scaffold{}\t{}\t{}\t{}\tW\t{}\t{}\t{}\t{}\n'.format(group+1, start, end, m, output_ID, 1, ctg_end-ctg_start+1, orient))

                    # output FASTA file for contigs
                    ctg_seq = ctg_seq_dict[old_ID][ctg_start-1:ctg_end]
                    f_ctg.write('>{}\n{}\n'.format(output_ID, ctg_seq))
                    
                    # output FASTA file for scaffolds
                    if n > 0:
                        scaffold_seq_list.append('N'*100)
                    if orient == '+':
                        scaffold_seq_list.append(ctg_seq)
                    else:
                        scaffold_seq_list.append(revcomp(ctg_seq))
                
                f_chr.write('>scaffold{}\n{}\n'.format(group+1, ''.join(scaffold_seq_list)))
    
    def output_files_with_gap_closing(gap_closed_group_list):

        aligned_pos_dict = dict()

        def get_aligned_pos(ctg_1, orient_1, ctg_2, orient_2):

            if (ctg_1, orient_1, ctg_2, orient_2) in aligned_pos_dict:
                return aligned_pos_dict[(ctg_1, orient_1, ctg_2, orient_2)]

            match_len_1, l1_1 = ctg_ovl_dict[(ctg_1, orient_1, ctg_2, orient_2)]
            match_len_2, l1_2 = ctg_ovl_dict[(ctg_2, reverse(orient_2), ctg_1, reverse(orient_1))]

            if new_old_dict:
                old_ctg_1 = new_old_dict[ctg_1] 
                old_ctg_2 = new_old_dict[ctg_2]
            else:
                old_ctg_1 = ctg_1
                old_ctg_2 = ctg_2 

            if orient_1 == '+':
                match_start_1 = l1_1 + 1
                match_end_1 = match_len_1 + l1_1
                ctg_len_1 = ctg_len_dict[old_ctg_1]
                assert match_end_1 == ctg_len_1
            else:
                match_start_1 = 1
                match_end_1 = match_len_1

            if orient_2 == '+':
                match_start_2 = 1
                match_end_2 = match_len_2
            else:
                match_start_2 = l1_2 + 1
                match_end_2 = match_len_2 + l1_2
                ctg_len_2 = ctg_len_dict[old_ctg_2]
                assert match_end_2 == ctg_len_2

            seq_1 = ctg_seq_dict[old_ctg_1][match_start_1-1:match_end_1]
            seq_2 = ctg_seq_dict[old_ctg_2][match_start_2-1:match_end_2]

            if orient_1 == '-':
                seq_1 = revcomp(seq_1)
            if orient_2 == '-':
                seq_2 = revcomp(seq_2)
            
            # perform sequence alignment
            aln = WavefrontAligner(seq_1)
            res = aln(seq_2, clip_cigar=True, min_aligned_bases_left=10, min_aligned_bases_right=10)
            
            if orient_1 == '+':
                aligned_pos_1 = match_start_1 + res.pattern_start
            else:
                aligned_pos_1 = match_end_1 - res.pattern_start
            if orient_2 == '+':
                aligned_pos_2 = match_start_2 + res.text_start
            else:
                aligned_pos_2 = match_end_2 - res.text_start
            
            print('Alignment for rescue: {} {} {} {} {} {} {} {} {} {} {}'.format(
                ctg_1, orient_1, ctg_2, orient_2, (res.pattern_end - res.pattern_start + 1) / match_len_1,
                res.pattern_start, res.pattern_end, res.text_start, res.text_end, aligned_pos_1, aligned_pos_2), file=sys.stderr)
            print('seq1 {} {}:'.format(ctg_1, old_ctg_1), seq_1, file=sys.stderr)
            print('seq2 {} {}:'.format(ctg_2, old_ctg_2), seq_2, file=sys.stderr)

            aligned_pos_dict[(ctg_1, orient_1, ctg_2, orient_2)] = (aligned_pos_1, aligned_pos_2)

            return aligned_pos_1, aligned_pos_2

        # output final assembly
        counter = dict()
        with open('after_copy.agp', 'wt') as f_agp, open('after_copy_ctg.fa', 'wt') as f_ctg, open('after_copy_chr.fa', 'wt') as f_chr, open('gap_closed.bed', 'wt') as fbed:
            for group, ctg_list in enumerate((ctg_list for ctg_list in gap_closed_group_list if ctg_list)):
                m = 0
                start = 0
                end = 0
                scaffold_seq_list = []
                for n, ctg_info_list in enumerate(ctg_list):
                    seq_len = 0
                    seq_list = []
                    output_ID_list = []
                    for i, ctg_info in enumerate(ctg_info_list):
                        ctg, orient = ctg_info[:2]
                        raw_ID = ctg_dict[ctg][4]
                        old_ID = new_old_dict[raw_ID] if new_old_dict else raw_ID
                        # get ctg_start and ctg_end
                        if len(ctg_info_list) == 1:
                            ctg_start, ctg_end = ctg_dict[ctg][-2:]
                        else:
                            # the first ctg
                            if i == 0:
                                if orient == '+':
                                    ctg_start = ctg_dict[ctg][-2]
                                    ctg_end = get_aligned_pos(ctg, orient, ctg_info_list[1][0], ctg_info_list[1][1])[0]
                                else:
                                    ctg_start = get_aligned_pos(ctg, orient, ctg_info_list[1][0], ctg_info_list[1][1])[0]
                                    ctg_end = ctg_dict[ctg][-1]
                            # the last ctg
                            elif i == len(ctg_info_list) - 1:
                                if orient == '+':
                                    ctg_start = get_aligned_pos(ctg_info_list[i-1][0], ctg_info_list[i-1][1], ctg, orient)[1] + 1
                                    ctg_end = ctg_dict[ctg][-1]
                                else:
                                    ctg_start = ctg_dict[ctg][-2]
                                    ctg_end = get_aligned_pos(ctg_info_list[i-1][0], ctg_info_list[i-1][1], ctg, orient)[1] - 1
                            # ctgs in the middle
                            else:
                                if orient == '+':
                                    ctg_start = get_aligned_pos(ctg_info_list[i-1][0], ctg_info_list[i-1][1], ctg, orient)[1] + 1
                                    ctg_end = get_aligned_pos(ctg, orient, ctg_info_list[i+1][0], ctg_info_list[i+1][1])[0]
                                else:
                                    ctg_start = get_aligned_pos(ctg, orient, ctg_info_list[i+1][0], ctg_info_list[i+1][1])[0]
                                    ctg_end = get_aligned_pos(ctg_info_list[i-1][0], ctg_info_list[i-1][1], ctg, orient)[1] - 1

                        seq_len += ctg_end - ctg_start + 1
                        if orient == '+':
                            seq_list.append(ctg_seq_dict[old_ID][ctg_start-1:ctg_end])
                        else:
                            seq_list.append(revcomp(ctg_seq_dict[old_ID][ctg_start-1:ctg_end]))
                        if ctg_info[-1] == 'copy' and ctg not in counter:
                            counter[ctg] = 1
                        elif ctg_info[-1] == 'copy':
                            counter[ctg] += 1
                        if len(ctg_info_list) > 1:
                            old_ID = '{}_{}'.format(old_ID, orient)
                        if ctg_info[-1] == 'copy':
                            output_ID_list.append('{}_{}_{}_copy{}'.format(old_ID, ctg_start, ctg_end, counter[ctg]))
                        else:
                            output_ID_list.append('{}_{}_{}'.format(old_ID, ctg_start, ctg_end))
                    output_ID = '__'.join(output_ID_list)

                    # output AGP file
                    if n > 0:
                        m += 1
                        start = end + 1
                        end = start + 99
                        f_agp.write('scaffold{}\t{}\t{}\t{}\tU\t100\tscaffold\tyes\tproximity_ligation\n'.format(group+1, start, end, m))
                    m += 1
                    start = end + 1
                    end = start + seq_len - 1
                    f_agp.write('scaffold{}\t{}\t{}\t{}\tW\t{}\t{}\t{}\t{}\n'.format(group+1, start, end, m, output_ID, 1, seq_len, orient))

                    gap_pos = 0
                    if len(seq_list) > 1:
                        for x in range(len(seq_list) - 1):
                            gap_pos += len(seq_list[x])
                            fbed.write('{0}\t{1}\t{2}\t{0}_{3}\t1000\t.\n'.format(output_ID, gap_pos-1, gap_pos+1, x+1))
                    
                    seq = ''.join(seq_list)
                    # output FASTA file for contigs
                    f_ctg.write('>{}\n{}\n'.format(output_ID, seq))
                    
                    # output FASTA file for scaffolds
                    if n > 0:
                        scaffold_seq_list.append('N'*100)
                    scaffold_seq_list.append(seq)
                
                f_chr.write('>scaffold{}\n{}\n'.format(group+1, ''.join(scaffold_seq_list)))

    # round1: output copied contigs
    for to_group, info_list in copied_collapsed_regions_dict.items():
        for copied_ctg_list, orientation, after_ctg, before_ctg in info_list:
            if orientation == '-':
                copied_ctg_list = [[ctg, reverse(orient), 'copy'] for ctg, orient in copied_ctg_list[::-1]]
            else:
                copied_ctg_list = [[ctg, orient, 'copy'] for ctg, orient in copied_ctg_list]
            index = group_list[to_group].index([after_ctg, ctg_dict[after_ctg][3]]) + 1 if after_ctg != 'NA' else group_list[to_group].index([before_ctg, ctg_dict[before_ctg][3]])
            group_list[to_group] = group_list[to_group][:index] + copied_ctg_list + group_list[to_group][index:]
    
    # round2: remove the contigs that are not moved
    for ctg_list in group_list:
        for ctg_info in ctg_list:
            if ctg_info[-1] != 'copy' and ctg_info[0] in moved_ctgs:
                ctg_list.remove(ctg_info)
    
    sub_G_list = []    
    # add contigs in the group to the subgraph
    for group, ctg_list in enumerate(group_list):
        subnodes = set()
        for ctg_info in ctg_list:
            ctg = ctg_info[0]
            subnodes.add((ctg, 'Head'))
            subnodes.add((ctg, 'Tail'))
        # add unanchored contigs to the subgraph
        sub_G = G.subgraph(subnodes | unanchored_ctg_semis)
        sub_G_list.append(sub_G)

    if args.optimize:
        optimize()

    if args.rescue:
        rescued_group_list = rescue()
        group_list = rescued_group_list
    
    if args.close_gaps:
        gap_closed_group_list = close_gaps()
        output_files_with_gap_closing(gap_closed_group_list)
    else:
        output_files(group_list)
    
    
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'assembly', help='.review.assembly file output by Juicebox')
    parser.add_argument(
        'collapsed_regions', help='annotated collapsed regions, an example can be found at: HapHiC/examples/collapsed_regions.txt')
    parser.add_argument(
        'gfa', help='.p_utg.noseq.gfa file output by hifiasm')
    parser.add_argument(
        'contigs', help='unitigs in FASTA format')
    parser.add_argument(
        'expected_nhaps', type=int, help='expected number of haplotypes')
    parser.add_argument(
        '--liftover', help='.liftover.agp file')
    parser.add_argument(
        '--optimize', default=False, action='store_true',
        help='optimize the orientation of adjancent contigs based on their overlapping relationships in the .gfa file, default: %(default)s')
    parser.add_argument(
        '--nrounds', default=3, type=int,
        help='the number of optimization rounds, default: %(default)s')
    parser.add_argument(
        '--distance', default=10, type=int,
        help='the threshold of distance (number of contigs) between adjacent contigs for optimization, default: %(default)s')
    parser.add_argument(
        '--rescue', default=False, action='store_true',
        help='rescue unanchored contigs between adjacent contigs, default: %(default)s')
    parser.add_argument(
        '--close_gaps', default=False, action='store_true',
        help='close gaps between adjacent contigs if they are correctly oriented and overlapped, default: %(default)s')

    args = parser.parse_args()
    
    ctg_dict, group_list, frag_dict, unanchored_ctg_semis = parse_assembly(args.assembly)
    
    ctg_region_dict, collapsed_region_list = parse_collapsed_regions(args.collapsed_regions, ctg_dict, group_list)
    
    if args.liftover:
        old_new_dict, new_old_dict = parse_liftover(args.liftover, ctg_dict, frag_dict)
    else:
        old_new_dict, new_old_dict = None, None
    
    G, frag_rd_dict, ctg_len_dict, ctg_ovl_dict = parse_gfa(args.gfa, frag_dict, old_new_dict)

    hap_thresholds = determine_haplotype_depths(G, frag_rd_dict, args.expected_nhaps)
    
    copied_collapsed_regions_dict, moved_ctgs = copy_collapsed_regions(G, collapsed_region_list, ctg_dict, group_list, unanchored_ctg_semis, frag_rd_dict, hap_thresholds, ctg_region_dict)

    ctg_seq_dict = get_contig_seqs(args.contigs)
    
    list_output_contigs(copied_collapsed_regions_dict, moved_ctgs, ctg_dict, group_list, new_old_dict, frag_rd_dict, ctg_seq_dict, ctg_len_dict, ctg_ovl_dict, unanchored_ctg_semis, G, args)


if __name__ == '__main__':
    main()
