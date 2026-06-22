#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created at: 2024-10-25 10:05


import os
import sys
import argparse
import gzip
import collections
import re

import networkx as nx
from itertools import combinations, product
from pywfa import WavefrontAligner


def select_ctg(a_dict, key, orient, from_to):
    if from_to == 'from':
        return a_dict[key][-1] if orient == '+' else a_dict[key][0]
    else:
        assert from_to == 'to'
        return a_dict[key][0] if orient == '+' else a_dict[key][-1]


def reverse(orient):
    if orient == '+':
        return '-'
    else:
        return '+'


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def get_old_id(ctg_id, new_old_dict):
    raw_id = ctg_id.split(':::')[0] if ':::' in ctg_id else ctg_id
    if new_old_dict:
        return new_old_dict[raw_id]
    return raw_id


def compute_overlap_alignment(ctg_1, orient_1, ctg_2, orient_2,
                               ctg_ovl_dict, ctg_len_dict, ctg_seq_dict, new_old_dict):
    match_len_1, l1_1 = ctg_ovl_dict[(ctg_1, orient_1, ctg_2, orient_2)]
    match_len_2, l1_2 = ctg_ovl_dict[(ctg_2, reverse(orient_2), ctg_1, reverse(orient_1))]

    old_ctg_1 = get_old_id(ctg_1, new_old_dict)
    old_ctg_2 = get_old_id(ctg_2, new_old_dict)

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

    total_M = sum(int(num) for num, op in re.findall(r'(\d+)([MIDNSHPB=X])', res.cigarstring) if op in ('M', '='))

    return aligned_pos_1, aligned_pos_2, total_M, match_len_1, match_len_2, res, seq_1, seq_2, old_ctg_1, old_ctg_2


def parse_assembly(assembly, n_chromosomes):

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

    chromosome_groups = set()
    if n_chromosomes > 0:
        group_lengths = []
        for g_idx, g_ctg_list in enumerate(group_list):
            total_bp = sum(ctg_dict[ctg][-1] - ctg_dict[ctg][-2] + 1 for ctg, *_ in g_ctg_list)
            group_lengths.append((total_bp, g_idx))
        group_lengths.sort(reverse=True)
        chromosome_groups = {g for _, g in group_lengths[:n_chromosomes]}
        for g_idx, g_ctg_list in enumerate(group_list):
            if g_idx not in chromosome_groups:
                for ctg, *_ in g_ctg_list:
                    unanchored_ctg_semis.add((ctg, 'Head'))
                    unanchored_ctg_semis.add((ctg, 'Tail'))

    return ctg_dict, group_list, frag_dict, unanchored_ctg_semis, chromosome_groups


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


def parse_gfa(gfa, frag_dict, old_new_dict, ctg_seq_dict, new_old_dict, aln_cov):

    if gfa.endswith('gz'):
        fopen = gzip.open
    else:
        fopen = open

    frag_rd_dict = dict()
    ctg_len_dict = dict()
    ctg_ovl_dict = dict()

    G = nx.Graph()
    l_lines = []
    alignment_cache = {}
    rejected_edges = set()
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

                ctg_ovl_dict[(from_ctg, from_orient, to_ctg, to_orient)] = (int(cols[5][:-1]), int(cols[6].split(':')[-1]))

                l_lines.append((from_ctg, from_orient, to_ctg, to_orient))

    for from_ctg, from_orient, to_ctg, to_orient in l_lines:
        reverse_key = (to_ctg, reverse(to_orient), from_ctg, reverse(from_orient))
        if reverse_key not in ctg_ovl_dict:
            raise KeyError(
                'Reverse overlap not found in ctg_ovl_dict for edge: {} {} -> {} {} (reverse_key: {} {} -> {} {})'.format(
                    from_ctg, from_orient, to_ctg, to_orient,
                    to_ctg, reverse(to_orient), from_ctg, reverse(from_orient)))

        aligned_pos_1, aligned_pos_2, total_M, match_len_1, match_len_2, res, seq_1, seq_2, old_ctg_1, old_ctg_2 = \
            compute_overlap_alignment(from_ctg, from_orient, to_ctg, to_orient,
                                      ctg_ovl_dict, ctg_len_dict, ctg_seq_dict, new_old_dict)
        min_match = min(match_len_1, match_len_2)
        ratio = total_M / min_match if min_match > 0 else 0
        alignment_cache[(from_ctg, from_orient, to_ctg, to_orient)] = \
            (aligned_pos_1, aligned_pos_2, match_len_1, match_len_2, res, seq_1, seq_2, old_ctg_1, old_ctg_2, total_M, ratio)
        print('GFA edge validation: {} {} {} {} total_M={} min_match={} ratio={:.4f} cigar={}'.format(
            from_ctg, from_orient, to_ctg, to_orient,
            total_M, min_match, ratio, res.cigarstring), file=sys.stderr)

    for from_ctg, from_orient, to_ctg, to_orient in l_lines:
        from_semi = 'Tail' if from_orient == '+' else 'Head'
        to_semi = 'Head' if to_orient == '+' else 'Tail'

        reverse_key = (to_ctg, reverse(to_orient), from_ctg, reverse(from_orient))
        ratio = alignment_cache[(from_ctg, from_orient, to_ctg, to_orient)][-1]
        rev_ratio = alignment_cache[reverse_key][-1]

        if ratio >= aln_cov and rev_ratio >= aln_cov:
            if not G.has_edge((from_ctg, from_semi), (to_ctg, to_semi)):
                G.add_edge((from_ctg, from_semi), (to_ctg, to_semi))
        else:
            rejected_edges.add((from_ctg, from_orient, to_ctg, to_orient))
            rejected_edges.add((to_ctg, reverse(to_orient), from_ctg, reverse(from_orient)))
            print('GFA edge validation: {} {} {} {} REJECTED (ratio {:.4f} rev_ratio {:.4f} < {:.4f})'.format(
                from_ctg, from_orient, to_ctg, to_orient,
                ratio, rev_ratio, aln_cov), file=sys.stderr)

    return G, frag_rd_dict, ctg_len_dict, alignment_cache, rejected_edges


def filter_and_output_gfa(gfa, frag_dict, old_new_dict, rejected_edges):
    if gfa.endswith('gz'):
        fopen = gzip.open
    else:
        fopen = open

    basename = os.path.basename(gfa)
    output_gfa = basename.replace('.gz', '') + '.filtered.gfa'
    if output_gfa == basename:
        output_gfa = basename + '.filtered.gfa'

    with fopen(gfa, 'rt') as fin, open(output_gfa, 'wt') as fout:
        for line in fin:
            if not line.startswith('L'):
                fout.write(line)
                continue

            cols = line.split()
            from_ctg, from_orient, to_ctg, to_orient = cols[1:5]

            if old_new_dict:
                from_ctg = select_ctg(old_new_dict, from_ctg, from_orient, from_to='from')
                to_ctg = select_ctg(old_new_dict, to_ctg, to_orient, from_to='to')

            if from_ctg in frag_dict:
                from_ctg = select_ctg(frag_dict, from_ctg, from_orient, from_to='from')

            if to_ctg in frag_dict:
                to_ctg = select_ctg(frag_dict, to_ctg, to_orient, from_to='to')

            edge_key = (from_ctg, from_orient, to_ctg, to_orient)
            if edge_key not in rejected_edges:
                fout.write(line)

    print('Filtered GFA written to: {}'.format(output_gfa), file=sys.stderr)
    return output_gfa


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
                for heter_frag in sorted(head_frags | tail_frags):
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


def has_valid_shortest_path(G, source, target, max_contigs):
    try:
        path = nx.shortest_path(G, source, target)
    except nx.NetworkXNoPath:
        return False
    inner = path[1:-1]
    if len(inner) % 2 != 0:
        return False
    if not all(inner[i][0] == inner[i+1][0] for i in range(0, len(inner), 2)):
        return False
    if len(inner) // 2 > max_contigs:
        return False
    return True


def has_valid_path(G, source, target, cutoff):
    if source not in G or target not in G:
        return False
    for path in nx.all_simple_paths(G, source, target, cutoff=cutoff):
        inner = path[1:-1]
        if len(inner) % 2 != 0:
            continue
        for i in range(0, len(inner), 2):
            if inner[i][0] != inner[i+1][0]:
                break
        else:
            return True
    return False


def copy_collapsed_regions(G, collapsed_region_list, ctg_dict, group_list, unanchored_ctg_semis, frag_rd_dict, hap_thresholds, ctg_region_dict):

    def guess_haplotype_depth(rd):

        for i in range(len(hap_thresholds) - 1):
            if hap_thresholds[i] <= rd < hap_thresholds[i+1]:
                return i + 1
        # the haplotype depth will not be larger than the expected number of haplotypes
        return i + 1

    semi_set = {'Head', 'Tail'}
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
            for (ctg1, semi1), (ctg2, semi2) in combinations(product(sorted(heter_ctg_set), semi_set), 2):
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
                candidates = sorted(connected_allelic_ctgs - moved_ctgs)
                if not candidates:
                    print('Warning: all allelic contigs in {} have been moved to other regions, skip'.format(connected_allelic_ctgs))
                    continue
                selected_ctg = candidates[0]
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
                continue
            if ctg in collapsed_ctgs:
                # all collapsed contigs will be copied without checking the depth
                copied_ctg_list.append((ctg, orient))
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


def identify_telomeric_seqs(ctg_seq_dict, repeat):

    repeat_rev_comp = revcomp(repeat)
    telomere_info = dict()
    for ctg_id, seq in ctg_seq_dict.items():
        types = []
        if seq[:1000].count(repeat) >= 5:
            types.append('telomeric_at_5')
        if seq[-1000:].count(repeat_rev_comp) >= 5:
            types.append('telomeric_at_3')
        if types:
            telomere_info[ctg_id] = types

    with open('telomere_annotation.txt', 'w') as f:
        for ctg_id, types in telomere_info.items():
            f.write('{}\t{}\n'.format(ctg_id, ','.join(types)))
    print('Telomere annotation written to telomere_annotation.txt', file=sys.stderr)
    print('Total contigs with telomeric repeats: {}'.format(len(telomere_info)), file=sys.stderr)

    return telomere_info


def list_output_contigs(copied_collapsed_regions_dict, moved_ctgs, ctg_dict, group_list, new_old_dict, frag_rd_dict, ctg_seq_dict, ctg_len_dict, unanchored_ctg_semis, chromosome_groups, G, args, telomere_info, alignment_cache):

    def is_terminus(ctg, ctg_pos):
        if ctg_pos == 1:
            return True
        old_ID = get_old_id(ctg, new_old_dict)
        if ctg_pos == ctg_len_dict[old_ID]:
            return True


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

                    for (src, src_terminus, new_orient_i), (tgt, tgt_terminus, new_orient_j) in product(
                        [(i_tail, ctg_end_i, '+'), (i_head, ctg_start_i, '-')],
                        [(j_head, ctg_start_j, '+'), (j_tail, ctg_end_j, '-')]):
                        if is_terminus(ctg_i, src_terminus) and is_terminus(ctg_j, tgt_terminus) and has_valid_shortest_path(sub_G, src, tgt, args.distance):
                            ctg_list[i][1] = new_orient_i
                            ctg_list[i+1][1] = new_orient_j
                            break
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

    def extend_telomeres():

        def has_telomere(ctg, telomere_type):
            old_ID = get_old_id(ctg, new_old_dict)
            return old_ID in telomere_info and telomere_type in telomere_info[old_ID]

        def find_telomere_candidates(G, start_node, ctg_dict, has_telomere_func, max_dist, side):
            source_ctg = start_node[0]
            candidates = {}
            stack = [(start_node, 0, {source_ctg})]
            while stack:
                node, dist, visited = stack.pop()
                if dist > max_dist:
                    continue
                for neighbor in G.neighbors(node):
                    nbr_ctg, nbr_semi = neighbor
                    if nbr_ctg in visited:
                        continue
                    assert nbr_ctg in ctg_dict, 'Contig {} not found in ctg_dict'.format(nbr_ctg)
                    nbr_ctg_len = ctg_dict[nbr_ctg][-1] - ctg_dict[nbr_ctg][-2] + 1
                    new_dist = dist + nbr_ctg_len
                    if new_dist > max_dist:
                        continue
                    is_telomere = False
                    if side == 'left':
                        if nbr_semi == 'Tail' and has_telomere_func(nbr_ctg, 'telomeric_at_5'):
                            candidates[nbr_ctg] = '+'
                            is_telomere = True
                        elif nbr_semi == 'Head' and has_telomere_func(nbr_ctg, 'telomeric_at_3'):
                            candidates[nbr_ctg] = '-'
                            is_telomere = True
                    else:
                        if nbr_semi == 'Head' and has_telomere_func(nbr_ctg, 'telomeric_at_3'):
                            candidates[nbr_ctg] = '+'
                            is_telomere = True
                        elif nbr_semi == 'Tail' and has_telomere_func(nbr_ctg, 'telomeric_at_5'):
                            candidates[nbr_ctg] = '-'
                            is_telomere = True
                    if not is_telomere:
                        other_end = (nbr_ctg, 'Head') if nbr_semi == 'Tail' else (nbr_ctg, 'Tail')
                        stack.append((other_end, new_dist, visited | {nbr_ctg}))
            return candidates

        print('Extend telomere: targeting top {} groups by length: {}'.format(
            len(chromosome_groups), sorted(g + 1 for g in chromosome_groups)), file=sys.stderr)

        max_dist = args.extend_telomere_distance
        extended_ctgs = set()

        for group, ctg_list in enumerate(group_list):
            if group not in chromosome_groups:
                continue
            first_ctg, first_orient = ctg_list[0][:2]
            expected_first = 'telomeric_at_5' if first_orient == '+' else 'telomeric_at_3'
            if not has_telomere(first_ctg, expected_first):
                left_node = (first_ctg, 'Head') if first_orient == '+' else (first_ctg, 'Tail')
                candidates = find_telomere_candidates(G, left_node, ctg_dict, has_telomere, max_dist, 'left')
                if len(candidates) == 1:
                    telomere_ctg, orient = next(iter(candidates.items()))
                    ctg_list.insert(0, [telomere_ctg, orient, 'copy'])
                    extended_ctgs.add(telomere_ctg)
                    print('Extend telomere: group {} prepend {} {} (distance <= {} bp)'.format(
                        group + 1, telomere_ctg, orient, max_dist), file=sys.stderr)
                elif len(candidates) > 1:
                    print('Extend telomere: group {} left side has {} candidates ({}) within {} bp, skip'.format(
                        group + 1, len(candidates), ','.join(sorted(candidates.keys())), max_dist), file=sys.stderr)

            last_ctg, last_orient = ctg_list[-1][:2]
            expected_last = 'telomeric_at_3' if last_orient == '+' else 'telomeric_at_5'
            if not has_telomere(last_ctg, expected_last):
                right_node = (last_ctg, 'Tail') if last_orient == '+' else (last_ctg, 'Head')
                candidates = find_telomere_candidates(G, right_node, ctg_dict, has_telomere, max_dist, 'right')
                if len(candidates) == 1:
                    telomere_ctg, orient = next(iter(candidates.items()))
                    ctg_list.append([telomere_ctg, orient, 'copy'])
                    extended_ctgs.add(telomere_ctg)
                    print('Extend telomere: group {} append {} {} (distance <= {} bp)'.format(
                        group + 1, telomere_ctg, orient, max_dist), file=sys.stderr)
                elif len(candidates) > 1:
                    print('Extend telomere: group {} right side has {} candidates ({}) within {} bp, skip'.format(
                        group + 1, len(candidates), ','.join(sorted(candidates.keys())), max_dist), file=sys.stderr)

        # Post-process extended telomere contigs:
        # 1) Single scan through chromosome groups to collect copy status per contig.
        # 2) If a contig has no non-copy version in any chromosome group,
        #    remove the 'copy' tag from its first chromosome occurrence.
        # 3) Single scan through unanchored groups to remove all extended contigs.
        if extended_ctgs:
            ctg_has_noncopy = {ctg: False for ctg in extended_ctgs}
            ctg_first_copy = {ctg: None for ctg in extended_ctgs}

            for g in sorted(chromosome_groups):
                for idx, ctg_info in enumerate(group_list[g]):
                    ctg = ctg_info[0]
                    if ctg not in extended_ctgs:
                        continue
                    is_copy = len(ctg_info) > 2 and ctg_info[-1] == 'copy'
                    if not is_copy:
                        ctg_has_noncopy[ctg] = True
                    elif ctg_first_copy[ctg] is None:
                        ctg_first_copy[ctg] = (g, idx)

            for ctg in extended_ctgs:
                if not ctg_has_noncopy[ctg] and ctg_first_copy[ctg] is not None:
                    g, idx = ctg_first_copy[ctg]
                    group_list[g][idx].pop()
                    print('Extend telomere: removed copy tag from {} in group {} (first chromosome occurrence)'.format(
                        ctg, g + 1), file=sys.stderr)

            for g in range(len(group_list)):
                if g not in chromosome_groups:
                    before = len(group_list[g])
                    group_list[g] = [ctg_info for ctg_info in group_list[g] if ctg_info[0] not in extended_ctgs]
                    removed = before - len(group_list[g])
                    if removed > 0:
                        print('Extend telomere: removed {} contig(s) from unanchored group {}'.format(
                            removed, g + 1), file=sys.stderr)

        return group_list

    def rescue():

        def update_rescued_group_list(sub_G, i_edge, j_edge, ctg_list, i, rescued_ctg_set):
            if i == 0:
                rescued_group_list[-1].append(ctg_list[i])
            valid_path = None
            valid_count = 0
            for path in nx.all_simple_paths(sub_G, i_edge, j_edge, cutoff=2 * args.distance + 1):
                inner = path[1:-1]
                if len(inner) % 2 != 0:
                    continue
                for n in range(len(inner) // 2):
                    if inner[2*n][0] != inner[2*n+1][0]:
                        break
                else:
                    valid_count += 1
                    if valid_path is None:
                        valid_path = path
                    if valid_count > 1:
                        rescued_group_list[-1].append(ctg_list[i+1])
                        print('More than one valid path: {} -> {}'.format(i_edge, j_edge))
                        return
            if valid_path is None:
                rescued_group_list[-1].append(ctg_list[i+1])
                print('No valid path: {} -> {}'.format(i_edge, j_edge))
                return
            path = valid_path
            rescued_ctgs = []
            if path and all([e in unanchored_ctg_semis for e in path[1:-1]]):
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

        ordered_groups = [g for g in range(len(group_list)) if g in chromosome_groups] + \
                         [g for g in range(len(group_list)) if g not in chromosome_groups]

        for group in ordered_groups:
            ctg_list = group_list[group]
            rescued_group_list.append([])

            if group not in chromosome_groups:
                for ctg_info in ctg_list:
                    ctg = ctg_info[0]
                    if ctg not in rescued_ctg_set:
                        rescued_group_list[-1].append(ctg_info)
                continue

            if len(ctg_list) == 1:
                rescued_group_list[-1].append(ctg_list[0])
                continue

            # filter subgraphs
            sub_G = sub_G_list[group].copy()

            ctg_nodes = set()
            for ctg, orient in [c[:2] for c in ctg_list]:
                ctg_nodes.add((ctg, 'Head'))
                ctg_nodes.add((ctg, 'Tail'))

            reachable_nodes = set()
            for component in nx.connected_components(sub_G):
                if component & ctg_nodes:
                    reachable_nodes |= component

            connected_pairs = []
            for i in range(len(ctg_list) - 1):
                ctg_i, orient_i = ctg_list[i][:2]
                ctg_j, orient_j = ctg_list[i+1][:2]
                ctg_start_i, ctg_end_i = ctg_dict[ctg_i][-2:]
                ctg_start_j, ctg_end_j = ctg_dict[ctg_j][-2:]
                i_head, i_tail = (ctg_i, 'Head'), (ctg_i, 'Tail')
                j_head, j_tail = (ctg_j, 'Head'), (ctg_j, 'Tail')

                src = i_tail if orient_i == '+' else i_head
                src_terminus = ctg_end_i if orient_i == '+' else ctg_start_i
                tgt = j_head if orient_j == '+' else j_tail
                tgt_terminus = ctg_start_j if orient_j == '+' else ctg_end_j

                if is_terminus(ctg_i, src_terminus) and is_terminus(ctg_j, tgt_terminus) and has_valid_path(sub_G, src, tgt, 2 * args.distance + 1) and not sub_G.has_edge(src, tgt):
                    connected_pairs.append((src, tgt))

            for ctg in sorted({node[0] for node in sub_G.nodes() if node in unanchored_ctg_semis}):
                if frag_rd_dict[ctg] == 1:
                    ctg_head = (ctg, 'Head')
                    ctg_tail = (ctg, 'Tail')

                    if ctg_head not in reachable_nodes:
                        sub_G.remove_node(ctg_head)
                        sub_G.remove_node(ctg_tail)
                        print('Removed {} due to frag_rd == 1 (isolated)'.format(ctg))
                        continue

                    if not connected_pairs:
                        sub_G.remove_node(ctg_head)
                        sub_G.remove_node(ctg_tail)
                        print('Removed {} due to frag_rd == 1 (no connected pairs)'.format(ctg))
                        continue

                    temp_G = sub_G.copy()
                    temp_G.remove_node(ctg_head)
                    temp_G.remove_node(ctg_tail)
                    would_break = any(not has_valid_path(temp_G, src, tgt, 2 * args.distance + 1) for src, tgt in connected_pairs)

                    if not would_break:
                        sub_G.remove_node(ctg_head)
                        sub_G.remove_node(ctg_tail)
                        print('Removed {} due to frag_rd == 1'.format(ctg))
                    else:
                        print('Kept {} despite frag_rd == 1 to preserve connectivity'.format(ctg))

            for i in range(len(ctg_list) - 1):
                ctg_i, orient_i = ctg_list[i][:2]
                ctg_j, orient_j = ctg_list[i+1][:2]

                ctg_start_i, ctg_end_i = ctg_dict[ctg_i][-2:]
                ctg_start_j, ctg_end_j = ctg_dict[ctg_j][-2:]

                i_head = (ctg_i, 'Head')
                i_tail = (ctg_i, 'Tail')
                j_head = (ctg_j, 'Head')
                j_tail = (ctg_j, 'Tail')

                src = i_tail if orient_i == '+' else i_head
                src_terminus = ctg_end_i if orient_i == '+' else ctg_start_i
                tgt = j_head if orient_j == '+' else j_tail
                tgt_terminus = ctg_start_j if orient_j == '+' else ctg_end_j

                if is_terminus(ctg_i, src_terminus) and is_terminus(ctg_j, tgt_terminus) and has_valid_path(sub_G, src, tgt, 2 * args.distance + 1) and not sub_G.has_edge(src, tgt):
                    rescued_ctgs = update_rescued_group_list(sub_G, src, tgt, ctg_list, i, rescued_ctg_set)
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

                src = i_tail if orient_i == '+' else i_head
                src_terminus = ctg_end_i if orient_i == '+' else ctg_start_i
                tgt = j_head if orient_j == '+' else j_tail
                tgt_terminus = ctg_start_j if orient_j == '+' else ctg_end_j

                if is_terminus(ctg_i, src_terminus) and is_terminus(ctg_j, tgt_terminus) and G.has_edge(src, tgt):
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

            key = (ctg_1, orient_1, ctg_2, orient_2)
            aligned_pos_1, aligned_pos_2, match_len_1, match_len_2, res, seq_1, seq_2, old_ctg_1, old_ctg_2, total_M, _ = \
                alignment_cache[key]

            print('Alignment for rescue: {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
                ctg_1, orient_1, ctg_2, orient_2, (res.pattern_end - res.pattern_start + 1) / match_len_1,
                res.pattern_start, res.pattern_end, res.text_start, res.text_end, aligned_pos_1, aligned_pos_2,
                total_M, min(match_len_1, match_len_2), res.cigarstring), file=sys.stderr)
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
                        fwd_seq = ctg_seq_dict[old_ID][ctg_start-1:ctg_end]
                        if orient == '+':
                            seq_list.append(fwd_seq)
                        else:
                            seq_list.append(revcomp(fwd_seq))
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
                    agp_orient = '+' if len(ctg_info_list) > 1 else orient
                    f_agp.write('scaffold{}\t{}\t{}\t{}\tW\t{}\t{}\t{}\t{}\n'.format(group+1, start, end, m, output_ID, 1, seq_len, agp_orient))

                    gap_pos = 0
                    if len(seq_list) > 1:
                        for i in range(len(seq_list) - 1):
                            gap_pos += len(seq_list[i])
                            fbed.write('{0}\t{1}\t{2}\t{0}_{3}\t1000\t.\n'.format(output_ID, gap_pos-1, gap_pos+1, i+1))

                    oriented_seq = ''.join(seq_list)
                    is_gap_closed = len(ctg_info_list) > 1
                    ctg_fasta_seq = oriented_seq if is_gap_closed else fwd_seq
                    # output FASTA file for contigs
                    f_ctg.write('>{}\n{}\n'.format(output_ID, ctg_fasta_seq))

                    # output FASTA file for scaffolds
                    if n > 0:
                        scaffold_seq_list.append('N'*100)
                    scaffold_seq_list.append(oriented_seq)

                f_chr.write('>scaffold{}\n{}\n'.format(group+1, ''.join(scaffold_seq_list)))

    # round1: insert copied/moved contigs into target groups
    for to_group, info_list in copied_collapsed_regions_dict.items():
        for copied_ctg_list, orientation, after_ctg, before_ctg in info_list:
            if orientation == '-':
                new_list = []
                for ctg, orient in copied_ctg_list[::-1]:
                    if ctg in moved_ctgs:
                        new_list.append([ctg, reverse(orient), 'move'])
                    else:
                        new_list.append([ctg, reverse(orient), 'copy'])
                copied_ctg_list = new_list
            else:
                new_list = []
                for ctg, orient in copied_ctg_list:
                    if ctg in moved_ctgs:
                        new_list.append([ctg, orient, 'move'])
                    else:
                        new_list.append([ctg, orient, 'copy'])
                copied_ctg_list = new_list
            index = group_list[to_group].index([after_ctg, ctg_dict[after_ctg][3]]) + 1 if after_ctg != 'NA' else group_list[to_group].index([before_ctg, ctg_dict[before_ctg][3]])
            group_list[to_group] = group_list[to_group][:index] + copied_ctg_list + group_list[to_group][index:]

    # round2: remove the contigs that are not kept (copy/move)
    for i, ctg_list in enumerate(group_list):
        group_list[i] = [ctg_info for ctg_info in ctg_list
                         if ctg_info[-1] in ('copy', 'move') or ctg_info[0] not in moved_ctgs]

    # round3: strip 'move' tag so moved contigs become normal entries
    for i, ctg_list in enumerate(group_list):
        for j, ctg_info in enumerate(ctg_list):
            if ctg_info[-1] == 'move':
                ctg_list[j] = ctg_info[:2]

    sub_G_list = []
    # add contigs in the group to the subgraph
    for _, ctg_list in enumerate(group_list):
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

    if args.extend_telomeres:
        group_list = extend_telomeres()

    if args.rescue:
        rescued_group_list = rescue()
        group_list = rescued_group_list

    group_list = [g for g in group_list if g]

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
        'n_chromosomes', type=int,
        help='expected number of chromosomes; groups are sorted by total length, top N are treated as chromosomes, the rest as unanchored contigs')
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
        help='the maximum number of intermediate contigs allowed between adjacent contigs for both optimization and rescue, default: %(default)s')
    parser.add_argument(
        '--rescue', default=False, action='store_true',
        help='rescue unanchored contigs between adjacent contigs, default: %(default)s')
    parser.add_argument(
        '--close-gaps', default=False, action='store_true',
        help='close gaps between adjacent contigs if they are correctly oriented and overlapped, default: %(default)s')
    parser.add_argument(
        '--extend-telomeres', default=False, action='store_true',
        help='extend telomeric contigs for the top N longest groups, default: %(default)s')
    parser.add_argument(
        '--extend-telomere-distance', type=int, default=1_000_000,
        help='maximum distance in bp to search for telomere contigs from group ends, default: %(default)s')
    parser.add_argument(
        '--repeat', default='CCCTAAA',
        help='telomeric repeat motif for identifying telomeric sequences at contig ends, default: %(default)s')
    parser.add_argument(
        '--aln-cov', type=float, default=0.5,
        help='alignment coverage threshold for filtering GFA edges, default: %(default)s')

    args = parser.parse_args()

    assert args.expected_nhaps > 1, 'expected_nhaps must be > 1'
    assert args.n_chromosomes > 1, 'n_chromosomes must be > 1'
    assert args.nrounds > 0, 'nrounds must be > 0'
    assert args.distance > 0, 'distance must be > 0'
    assert args.extend_telomere_distance > 0, 'extend-telomere-distance must be > 0'
    assert args.aln_cov > 0, 'aln-cov must be > 0'

    ctg_dict, group_list, frag_dict, unanchored_ctg_semis, chromosome_groups = parse_assembly(args.assembly, args.n_chromosomes)

    ctg_region_dict, collapsed_region_list = parse_collapsed_regions(args.collapsed_regions, ctg_dict, group_list)

    if args.liftover:
        old_new_dict, new_old_dict = parse_liftover(args.liftover, ctg_dict, frag_dict)
    else:
        old_new_dict, new_old_dict = None, None

    ctg_seq_dict = get_contig_seqs(args.contigs)

    telomere_info = identify_telomeric_seqs(ctg_seq_dict, args.repeat)

    G, frag_rd_dict, ctg_len_dict, alignment_cache, rejected_edges = parse_gfa(args.gfa, frag_dict, old_new_dict, ctg_seq_dict, new_old_dict, args.aln_cov)

    filter_and_output_gfa(args.gfa, frag_dict, old_new_dict, rejected_edges)

    hap_thresholds = determine_haplotype_depths(G, frag_rd_dict, args.expected_nhaps)

    copied_collapsed_regions_dict, moved_ctgs = copy_collapsed_regions(G, collapsed_region_list, ctg_dict, group_list, unanchored_ctg_semis, frag_rd_dict, hap_thresholds, ctg_region_dict)

    list_output_contigs(copied_collapsed_regions_dict, moved_ctgs, ctg_dict, group_list, new_old_dict, frag_rd_dict, ctg_seq_dict, ctg_len_dict, unanchored_ctg_semis, chromosome_groups, G, args, telomere_info, alignment_cache)


if __name__ == '__main__':
    main()
