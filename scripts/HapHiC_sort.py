#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Email: zengxf@sustech.edu.cn

import os
import sys
import subprocess
import logging
import time
import gc

from HapHiC_build import parse_tours

import pickle
import argparse
from collections import defaultdict
from multiprocessing import Pool
from itertools import combinations, product

from numpy import float32, zeros, hstack
from networkx import Graph, connected_components, shortest_path
from networkx import tree as nxtree
from scipy.sparse import coo_matrix

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_fasta(fasta):

    logger.info('Parsing fasta file...')

    fa_dict = dict()

    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                fa_dict[ctg] = 0
            else:
                fa_dict[ctg] += len(line.strip())
    return fa_dict


def dict_to_matrix(dict_, shape, add_self_loops=False):

    # This function is different from the one in HapHiC_cluster, the keys of
    # input dict here, are the indexes of the final matrix

    # prepare three lists for dict_to_matrix conversion
    row, col, data = list(), list(), list()

    for ctg_pair, links in dict_.items():

        row.append(ctg_pair[0])
        col.append(ctg_pair[1])
        # diagonal symmetry
        row.append(ctg_pair[1])
        col.append(ctg_pair[0])

        data.append(links)
        # diagonal symmetry
        data.append(links)

    # add self loops (for Markov Clustering)
    if add_self_loops:
        for n in range(shape):
            row.append(n)
            col.append(n)
            data.append(1)

    # convert dict to matrix using coo_matrix (sparse matrix)
    return coo_matrix((data, (row, col)), shape=(shape, shape), dtype=float32).toarray()


def parse_group(group_file, clm_dir, quick_view):

    # find the corresponding clm file first
    # the prefix of the group file (often group name)
    prefix = os.path.splitext(os.path.basename(group_file))[0]
    clm = '{}/{}.clm'.format(clm_dir, prefix)

    if not quick_view and not os.path.exists(clm):
        raise IOError('Clm file check failed: CANNOT find corresponding clm file '
                'in {} for group file {}'.format(clm_dir, group_file))

    ctgs = list()

    with open(group_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            ctgs.append((cols[0], int(cols[2])))

    # sort by ctg length
    ctgs.sort(key=lambda x: x[1], reverse=True)

    return ctgs, clm, prefix


def get_sub_HT_dict(ctgs, HT_link_dict):

    suffix = ['_H', '_T']
    suffix_comb = list(product(suffix, repeat=2))

    index = 0
    sub_HT_dict = defaultdict(int)
    HT_index_dict = dict()

    for ctg_1, ctg_2 in combinations(ctgs, 2):
        ctg_1, ctg_2 = sorted([ctg_1, ctg_2])
        for suf_1, suf_2 in suffix_comb:
            HT_1, HT_2 = ctg_1 + suf_1, ctg_2 + suf_2
            if (HT_1, HT_2) in HT_link_dict:
                links = HT_link_dict[(HT_1, HT_2)]
            else:
                links = 0
            for HT in (HT_1, HT_2):
                # assign index
                if HT not in HT_index_dict:
                    HT_index_dict[HT] = index
                    index += 1
            # construct sub_HT_dict
            if links:
                sub_HT_dict[(HT_index_dict[HT_1], HT_index_dict[HT_2])] = links

    return sub_HT_dict, HT_index_dict


def get_len(HT, fa_dict):

    # a contig
    if isinstance(HT, str):
        # half length of a contig
        return fa_dict[HT.rsplit('_', 1)[0]]/2
    # a scaffold
    else:
        assert isinstance(HT, tuple)
        return sum([get_len(ht, fa_dict) for ht in HT])


def get_density_graph(sub_HT_matrix, shape, index_HT_dict, fa_dict, flank_HT_dict, density_cal_method):

    def get_HT_len(HT):
        if HT in flank_HT_dict:
            return flank_HT_dict[HT][1]
        else:
            return get_len(HT, fa_dict)

    # construct a length matrix for the calculation of link density
    HT_len_sum_dict = dict()
    indexes = index_HT_dict.keys()

    for i_1, i_2 in combinations(indexes, 2):

        HT_1 = index_HT_dict[i_1]
        HT_2 = index_HT_dict[i_2]
        HT_1_len = get_HT_len(HT_1)
        HT_2_len = get_HT_len(HT_2)

        if density_cal_method == 'sum':
            len_sum = HT_1_len + HT_2_len
        elif density_cal_method == 'multiplication':
            len_sum = (HT_1_len * HT_2_len)
        elif density_cal_method == 'geometric_mean':
            len_sum = (HT_1_len * HT_2_len)**0.5

        HT_len_sum_dict[(i_1, i_2)] = len_sum

    # add self loops to prevent division by zero error
    HT_len_sum_matrix = dict_to_matrix(HT_len_sum_dict, shape, add_self_loops=True)

    # get density graph
    density_graph = sub_HT_matrix / HT_len_sum_matrix

    return density_graph


def get_unfiltered_confidence_graph(shape, index_pairs, sub_HT_dict, density_graph):

    # generate a empty matrix
    confidence_graph = zeros((shape, shape))

    # calculate confidence

    # Confidence is the ratio of the link density of non-sister edge pair A-B,
    # and the link density of the second-largest non-sister edge incident on either A or B.
    # Only if density is not zero, confidence needs to be calculated.
    for i_1, i_2 in sub_HT_dict:

        density = density_graph[i_1, i_2]

        # calculate confidence
        # all densities of non-sister edge incident on either A or B, but density of A-B is only counted once
        all_densities = hstack((density_graph[i_1,:i_2], density_graph[i_1,i_2+1:], density_graph[:,i_2]))
        # get the index of maximum value
        argmax = all_densities.argmax()
        # remove the maximum value (based on the index)
        all_densities = hstack((all_densities[:argmax], all_densities[argmax+1:]))
        # the maximum value in the rest is the second largest density
        second_largest_density = all_densities.max()
        # get the confidence
        if density == 0:
            confidence = 0
        elif second_largest_density == 0:
            confidence = 2
        else:
            confidence = density / second_largest_density

        # diagonal symmetry
        confidence_graph[i_1, i_2] = confidence
        confidence_graph[i_2, i_1] = confidence

    # assign sister edges a weight of 2*MAXS, where MAXS is the maximum density of
    # all the non-sister edges in the unfiltered confidence graph

    # The confidence graph here does not includes the weight between sister edges,
    # so the max value of the graph is the MAXS.
    maxs = confidence_graph.max()

    # assign weights for sister edges
    for index_H, index_T in index_pairs:
        # diagonal symmetry
        sister_weight = 2 * maxs if maxs > 1 else 2
        confidence_graph[index_H, index_T] = sister_weight
        confidence_graph[index_T, index_H] = sister_weight

    return confidence_graph, maxs


def filter_confidence_graph(confidence_graph, sub_HT_dict, confidence_cutoff):

    # Only confidence greater than 1 is reliable otherwise unreliable.
    # I also tried to use the confidence <= 1 to link at least a pair of
    # non-sisiter edges in each round, but this could result in forked trees.
    for i_1, i_2 in sub_HT_dict:
        if confidence_graph[i_1, i_2] <= confidence_cutoff:
            confidence_graph[i_1, i_2] = 0
            confidence_graph[i_2, i_1] = 0


def get_HT_ends(HT):

    # make contig and scaffold a same format
    if isinstance(HT, str):
        return (HT,), HT, HT
    else:
        assert isinstance(HT, tuple)
        return HT, HT[0], HT[-1]


def split_new_scaffold(path, fa_dict, index_HT_dict, known_adjacency):

    # parse path
    path_len = len(path)
    # path length should be a even number
    assert path_len % 2 == 0

    sorted_path = list()

    for n in range(path_len//2):
        HT_1 = index_HT_dict[path[2*n]]
        HT_2 = index_HT_dict[path[2*n+1]]

        HT_1, HT_1_left,  HT_1_right = get_HT_ends(HT_1)
        HT_2, HT_2_left, HT_2_right = get_HT_ends(HT_2)

        # to figure out how they are linked
        # reverse 1
        if tuple(sorted([HT_1_left, HT_2_left])) in known_adjacency:
            sorted_path.extend(HT_1[::-1])
            sorted_path.extend(HT_2)
        # reverse 2
        elif tuple(sorted([HT_1_right, HT_2_right])) in known_adjacency:
            sorted_path.extend(HT_1)
            sorted_path.extend(HT_2[::-1])
        # reverse both 1 & 2
        elif tuple(sorted([HT_1_left, HT_2_right])) in known_adjacency:
            sorted_path.extend(HT_1[::-1])
            sorted_path.extend(HT_2[::-1])
        # do nothing
        else:
            assert tuple(sorted([HT_1_right, HT_2_left])) in known_adjacency
            sorted_path.extend(HT_1)
            sorted_path.extend(HT_2)

    # calculate half length of superscaffold
    half_scaffold_len = get_len(tuple([index_HT_dict[i] for i in path]), fa_dict) / 2

    accumulated_len = 0
    len_difference_list = list()

    for i, HT in enumerate(sorted_path):
        assert isinstance(HT, str)
        accumulated_len += get_len(HT, fa_dict)
        len_difference_list.append((i, abs(accumulated_len - half_scaffold_len)))

    len_difference_list.sort(key=lambda x: x[1])

    split_index = len_difference_list[0][0] + 1

    split_path_1 = tuple(sorted_path[:split_index])
    split_path_2 = tuple(sorted_path[split_index:])

    adjacent_HT = tuple(sorted([split_path_1[-1], split_path_2[0]]))

    if adjacent_HT not in known_adjacency:
        known_adjacency.add(adjacent_HT)

    return split_path_1, split_path_2


def format_HTs(HT):

    if isinstance(HT, str):
        return [HT]
    else:
        assert isinstance(HT, tuple)
        return HT


def update(path_list, old_sub_HT_matrix, index_HT_dict, HT_index_dict, flank_HT_dict, fa_dict, known_adjacency, flank):

    def write_output_path_list(HT):
        if isinstance(HT, tuple):
            output_path_list[-1].extend(HT)
        else:
            assert isinstance(HT, str)
            output_path_list[-1].append(HT)

    def get_flank_HT(split_path, order):

        rest_len = get_len(split_path, fa_dict)
        if rest_len > flank:
            for m, HT in enumerate(split_path[::order]):
                HT_len = get_len(HT, fa_dict)
                if rest_len - HT_len > flank:
                    rest_len -= HT_len
                else:
                    break
            if m == 0:
                rest_path = split_path
            elif order == 1:
                rest_path = split_path[:-m]
            else:
                rest_path = split_path[m:]
            flank_HT_dict[split_path] = (rest_path, get_len(rest_path, fa_dict))


    index_pairs = list()
    new_index_HT_dict = dict()
    output_path_list = list()

    for n, path in enumerate(path_list):

        # generate new indexes for the sister edges of each path
        i_1, i_2 = 2*n, 2*n+1
        index_pairs.append((i_2, i_1))

        # a list to interpret the paths
        output_path_list.append([])

        # path with two edges, could be a contig or a scaffolds that have been linked in prior rounds
        if len(path) == 2:

            HT_left, HT_right = index_HT_dict[path[0]], index_HT_dict[path[1]]

            new_index_HT_dict[i_1] = HT_left
            new_index_HT_dict[i_2] = HT_right

            write_output_path_list(HT_left)
            write_output_path_list(HT_right)

        # newly linked scaffolds in this round
        else:

            # each newly linked scaffold is split into two halves approximately
            split_path_left, split_path_right = split_new_scaffold(path, fa_dict, index_HT_dict, known_adjacency)

            new_index_HT_dict[i_1] = split_path_left
            new_index_HT_dict[i_2] = split_path_right

            output_path_list[-1].extend(split_path_left)
            output_path_list[-1].extend(split_path_right)

            if flank:
                get_flank_HT(split_path_left, -1)
                get_flank_HT(split_path_right, 1)

    # update sub_HT_dict
    sub_HT_dict = defaultdict(int)

    for i_1 in new_index_HT_dict:
        for i_2 in new_index_HT_dict:

            if i_1 <= i_2 or (i_1, i_2) in index_pairs:
                continue

            HT_1 = new_index_HT_dict[i_1]
            HT_2 = new_index_HT_dict[i_2]

            if HT_1 in flank_HT_dict:
                HT_1 = format_HTs(flank_HT_dict[HT_1][0])
            else:
                HT_1 = format_HTs(HT_1)
            if HT_2 in flank_HT_dict:
                HT_2 = format_HTs(flank_HT_dict[HT_2][0])
            else:
                HT_2 = format_HTs(HT_2)

            for ht_1, ht_2 in product(HT_1, HT_2):

                old_i_1 = HT_index_dict[ht_1]
                old_i_2 = HT_index_dict[ht_2]

                links = old_sub_HT_matrix[old_i_1, old_i_2]

                if links:
                    sub_HT_dict[(i_1, i_2)] += links

    return new_index_HT_dict, sub_HT_dict, index_pairs, output_path_list


def output_tour_file(output_path_list, prefix):

    tour_file = '{}.tour'.format(prefix)
    with open(tour_file, 'w') as ftour:
        ftour.write('>INIT\n')
        output_ctg_list = list()
        for path_list in output_path_list:
            for HT in path_list[::2]:
                ctg, tag = HT.rsplit('_', 1)
                if tag == 'H':
                    output_ctg_list.append(ctg+'+')
                else:
                    output_ctg_list.append(ctg+'-')
        ftour.write('{}\n'.format(' '.join(output_ctg_list)))


def remove_shortest_path(index_pairs, sub_HT_dict, density_graph):

    # update index_pairs
    shortest_i_1, shortest_i_2 = index_pairs.pop(-1)

    # update sub_HT_dict
    for i_1, i_2 in sub_HT_dict.copy():
        if i_1 in {shortest_i_1, shortest_i_2} or i_2 in {shortest_i_1, shortest_i_2}:
            sub_HT_dict.pop((i_1, i_2))

    # update density_graph
    return density_graph[:-2,:-2]


def fast_sort(args, fa_dict, group_specific_data, prefix):

    logger.info('[{}] Performing fast sorting...'.format(prefix))

    ctg_info_list, ctgs, sub_HT_dict, HT_index_dict = group_specific_data

    # file checks
    logger.info('[{}] Checking the content of input group file...'.format(prefix))
    # contig name not match
    for ctg, ctg_len in ctg_info_list:
        if ctg not in fa_dict:
            logger.error('[{}] CANNOT find contig {} in the FASTA file'.format(prefix, ctg))
            raise RuntimeError('Group file check failed: [{}] CANNOT find contig {} in the FASTA file'.format(prefix, ctg))
        elif ctg_len != fa_dict[ctg]:
            logger.error('[{}] Length of contig {} in the group file ({} bp) does NOT match the FASTA file ({} bp)'.format(
                prefix, ctg, ctg_len, fa_dict[ctg]))
            raise RuntimeError('Group file check failed: [{}] Length of contig {} in the group file ({} bp) does NOT match the FASTA file ({} bp)'.format(
                prefix, ctg, ctg_len, fa_dict[ctg]))
    # if only 1 contig, neither ordering nor orientation is required
    if len(ctg_info_list) == 1:
        logger.info('[{}] Only 1 contig in the group file, neither ordering nor orientation is required'.format(prefix))
        return [[ctg_info_list[0][0]+'_H', ctg_info_list[0][0]+'_T']], True
    # empty group file
    elif not ctg_info_list:
        logger.error('[{}] CANNOT find any contigs in the group file'.format(prefix))
        raise RuntimeError('Group file check failed: [{}] CANNOT find any contigs in the group file'.format(prefix))

    # prepare
    # get the shape of the matrix
    shape = len(ctgs)*2
    # convert HT_index_dict to index_HT_dict
    index_HT_dict = {index: HT for HT, index in HT_index_dict.items()}

    # get index pairs of sister edges
    index_pairs = list()
    # initialize output_path_list and set known_adjacency
    output_path_list = list()
    known_adjacency = set()
    for ctg in ctgs:
        ctg_H, ctg_T = ctg+'_H', ctg+'_T'
        output_path_list.append([ctg_H, ctg_T])
        index_pairs.append((HT_index_dict[ctg_H], HT_index_dict[ctg_T]))
        known_adjacency.add(tuple(sorted([ctg_H, ctg_T])))

    # a list to store all shortest paths in the forest
    path_list = ctgs.copy()
    # flank_HT_dict is used to store the flank region length of scaffolds
    # In the first round of iterations, it's empty.
    # flank_HT_dict will be updated at the end of each round.
    flank_HT_dict = dict()

    # a set to store removed contigs / scaffolds
    removed_path_list = list()
    skip_density_graph_construction = False
    # record number of rounds
    r = 0

    logger.info('[{}] Starting fast sorting iterations...'.format(prefix))

    while len(path_list) != 1:

        r += 1

        if not skip_density_graph_construction:
            # convert sub_HT_dict to sub_HT_matrix
            sub_HT_matrix  = dict_to_matrix(sub_HT_dict, shape)

            # in the first round, record the oldest version of sub_HT_matrix for ctg-ctg link calculation
            if r == 1:
                old_sub_HT_matrix = sub_HT_matrix

            # get density graph
            density_graph = get_density_graph(sub_HT_matrix, shape, index_HT_dict, fa_dict, flank_HT_dict, args.density_cal_method)

        # get unfiltered confidence graph
        confidence_graph, maxs = get_unfiltered_confidence_graph(
                shape, index_pairs, sub_HT_dict, density_graph)

        # nothing new, and number of paths is greater than 2
        if maxs <= args.confidence_cutoff and len(path_list) > 2:
            shape -= 2
            density_graph = remove_shortest_path(index_pairs, sub_HT_dict, density_graph)
            removed_path_list.append(output_path_list.pop(-1))
            path_list.pop(-1)
            logger.debug('[{}] Round {}, MAXS {}'.format(prefix, r, maxs))
            logger.debug('[{}] Path {}, Length {}, removed'.format(
                prefix, len(output_path_list)+1, get_len(tuple(removed_path_list[-1]), fa_dict)))
            skip_density_graph_construction = True
            continue
        # nothing new and only two paths, terminate iterations
        elif maxs <= args.confidence_cutoff:
            assert len(path_list) == 2
            break

        skip_density_graph_construction = False

        # filter confidence graph
        filter_confidence_graph(confidence_graph, sub_HT_dict, args.confidence_cutoff)

        # use Kruskal's algorithm to construct a maximal spanning forest
        forest = nxtree.maximum_spanning_tree(Graph(confidence_graph), algorithm='kruskal')

        # clear the result of previous round
        path_list.clear()

        # get subnodes in each tree
        for n, subnodes in enumerate(connected_components(forest)):
            tree = forest.subgraph(subnodes)
            # a list to store the head (source) and the tail (target) of the path
            ends_list = list()
            for node, degree in tree.degree():
                # the source and target should have degree 1
                if degree == 1:
                    ends_list.append(node)

            # there should be two nodes with degree 1 (will it be a circle?)
            assert len(ends_list) == 2

            # store the path
            source, target = ends_list
            path = dict(shortest_path(tree))[source][target]
            path_len = get_len(tuple([index_HT_dict[i] for i in path]), fa_dict)
            path_list.append((path, path_len))

        # sort path_list by length
        path_list = [path for path, length in sorted(path_list, key=lambda x: x[1], reverse=True)]

        # Updates for the next iteration
        # update shape
        shape = len(path_list*2)

        # update sub_HT_dict, index_HT_dict, index_pairs, flank_HT_dict
        flank = args.flanking_region * 1000
        index_HT_dict, sub_HT_dict, index_pairs, output_path_list = update(
                path_list, old_sub_HT_matrix, index_HT_dict, HT_index_dict,
                flank_HT_dict, fa_dict, known_adjacency, flank)

        # debug
        logger.debug('[{}] Round {}, MAXS {}'.format(prefix, r, maxs))
        for n, path in enumerate(output_path_list, 1):
            logger.debug('[{}] Path {}, Length {}: {}'.format(
                prefix, n, get_len(tuple(path), fa_dict), '->'.join(path)))

    output_path_list.extend(removed_path_list[::-1])

    return output_path_list, False


def run_allhic_optimization(args, group, prefix, clm, allhic):

    logger.info('[{}] Performing ALLHiC optimization...'.format(prefix))

    # prepare parameters & command
    cmd_list = [
            allhic, 'optimize', group, clm,
            '--mutapb', str(args.mutprob),
            '--ngen', str(args.ngen),
            '--npop', str(args.npop),
            '--seed', str(args.seed)
            ]

    if not args.skip_fast_sort:
        cmd_list.append('--resume')

    if args.skipGA:
        cmd_list.append('--skipGA')

    # run ALLHiC optimization
    subprocess.run(
            cmd_list, check=True,
            stdout=open('allhic_logs/{}.allhic_stdout.log'.format(prefix), 'w'),
            stderr=open('allhic_logs/{}.allhic_stderr.log'.format(prefix), 'w')
            )


def compare_fast_sort_and_allhic(prefix, fa_dict):

    def find_lis(compare_list, order_len_dict, forward=True):

        order_list = []
        if forward:
            order_list = [order for order in compare_list if order > 0]
        else:
            order_list = [order for order in compare_list if order < 0]

        if not order_list:
            return 0

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

        return max_sum

    fast_sort_tour = '{}.tour.sav'.format(prefix)
    allhic_tour = '{}.tour'.format(prefix)

    fast_sort_ctg_list, fast_sort_ori_list = [], []
    for ctg, ori in list(parse_tours([fast_sort_tour], fa_dict)[0].values())[0]:
        fast_sort_ctg_list.append(ctg)
        fast_sort_ori_list.append(ori)

    ctg_len_list = [fa_dict[ctg] for ctg in fast_sort_ctg_list]
    group_len = sum(ctg_len_list)
    # if contigs are short enough, always choose allhic
    group_ctg_len_ratio = group_len / max(ctg_len_list)
    if group_ctg_len_ratio > 50:
        logger.info('{}: choose allhic optimization (group length / longest contig = {})'.format(prefix, group_ctg_len_ratio))
        return False

    allhic_ctg_list, allhic_ori_list = [], []
    for ctg, ori in list(parse_tours([allhic_tour], fa_dict)[0].values())[0]:
        allhic_ctg_list.append(ctg)
        allhic_ori_list.append(ori)

    max_lis_len_ratio = 0
    for n in range(len(fast_sort_ctg_list) - 1):
        compare_list, order_len_dict = [], dict()
        for i, ctg in enumerate(fast_sort_ctg_list):
            j = allhic_ctg_list.index(ctg)
            if fast_sort_ori_list[i] == allhic_ori_list[j]:
                compare_list.append(j+1)
                order_len_dict[j+1] = fa_dict[ctg]
            else:
                compare_list.append(-j-1)
                order_len_dict[-j-1] = fa_dict[ctg]

        max_sum_f = find_lis(compare_list, order_len_dict, forward=True)
        max_sum_r = find_lis(compare_list, order_len_dict, forward=False)

        max_sum = max(max_sum_f, max_sum_r)

        lis_len_ratio = max_sum / group_len

        if lis_len_ratio >= 0.9:
            logger.info('{}: choose allhic optimization (LIS length / group length = {})'.format(prefix, max_sum / group_len))
            return False
        else:
            if lis_len_ratio > max_lis_len_ratio:
                max_lis_len_ratio = lis_len_ratio
            fast_sort_ctg_list = fast_sort_ctg_list[1:] + [fast_sort_ctg_list[0]]
            fast_sort_ori_list = fast_sort_ori_list[1:] + [fast_sort_ori_list[0]]

    logger.info('{}: choose fast sorting (maximum LIS length / group length = {})'.format(prefix, lis_len_ratio))
    return True


def run_haphic_sorting(args, group, fa_dict, group_specific_data, group_param, allhic):

    prefix, clm = group_param

    only_one_contig = False

    output_sav = False

    if not args.skip_fast_sort:

        # fast sort contigs before ALLHiC optimization
        output_path_list, only_one_contig = fast_sort(args, fa_dict, group_specific_data, prefix)

        # output .tour file
        output_tour_file(output_path_list, prefix)

    # call a modified version of ALLHiC for hotstart optimization (--resume)
    if not args.skip_allhic and not only_one_contig:
        run_allhic_optimization(args, group, prefix, clm, allhic)
        # comaprison is needed only if fast sort is not skipped
        if not args.skip_fast_sort:
            output_sav = compare_fast_sort_and_allhic(prefix, fa_dict)

    if output_sav:
        os.symlink('../{}.tour.sav'.format(prefix), os.path.join('final_tours', '{}.tour'.format(prefix)))
    else:
        os.symlink('../{}.tour'.format(prefix), os.path.join('final_tours', '{}.tour'.format(prefix)))


def check_exceptions(result_list):

    nerrors = 0

    for result in result_list:
        try:
            result.get()
        except Exception as e:
            nerrors += 1
            logger.error(e)

    if nerrors:
        raise Exception('{} exception(s) detected, please check the log above.'.format(nerrors))


def parse_arguments():

    parser = argparse.ArgumentParser(prog='haphic sort')

    # Parameters for parsing input files and pipeline control
    input_group = parser.add_argument_group('>>> Parameters for parsing input files and pipeline control')
    input_group.add_argument(
            'fasta', help='draft genome in FASTA format. Use `corrected_asm.fa` generated in the clustering step when `--correct_nrounds` was set')
    input_group.add_argument(
            'HT_links', help='`HT_links.pkl` generated in the clustering step')
    input_group.add_argument(
            'clm_dir', help='directory containing split clm files generated in the reassignment step (`split_clms`). '
            'The program will find the corresponding clm file for each input group file')
    input_group.add_argument(
            'groups', nargs='+', help='one or more group files generated in the reassignment step (`group*.txt` in the '
            'directory `final_groups`)')
    input_group.add_argument(
            '--quick_view', default=False, action='store_true',
            help="in quick view mode, HapHiC will skip clustering and reassignment steps, and order and orient all contigs with fast sorting, default: %(default)s. "
            "This is helpful when you encounter difficulties in generating ideal clusters or when you are unsure of the exact number of chromosomes")

    # Parameters for fast sorting
    fast_sort_group = parser.add_argument_group('>>> Parameters for fast sorting')
    fast_sort_group.add_argument(
            '--skip_fast_sort', default=False, action='store_true',
            help='skip fast sorting and only run ALLHiC optimization, default: %(default)s')
    fast_sort_group.add_argument(
            '--flanking_region', type=int, default=0,
            help='consider only the Hi-C links from the flanking regions (ends) of each contig / scaffold during fast sorting (unit: kbp). By default, this parameter is set to '
            '0 to use the whole contigs or scaffolds')
    fast_sort_group.add_argument(
            '--density_cal_method', choices={'multiplication', 'sum', 'geometric_mean'}, default='multiplication',
            help='method for Hi-C link density calculation during fast sorting, default: %(default)s')
    fast_sort_group.add_argument(
            '--confidence_cutoff', type=float, default=1,
            help='cutoff for confidence filtering, default: %(default)s')

    # Parameters for ALLHiC optimization
    allhic_group = parser.add_argument_group('>>> Parameters for ALLHiC optimization')
    allhic_group.add_argument(
            '--skip_allhic', default=False, action='store_true',
            help='skip the entire ALLHiC optimization step, default: %(default)s')
    allhic_group.add_argument(
            '--skipGA', default=False, action='store_true',
            help='skip the genetic algorithm optimization step in ALLHiC, default: %(default)s')
    allhic_group.add_argument(
            '--mutprob', type=float, default=0.2,
            help='mutation probability in the genetic algorithm, default: %(default)s')
    allhic_group.add_argument(
            '--ngen', type=int, default=5000,
            help='number of generations for convergence, default: %(default)s')
    allhic_group.add_argument(
            '--npop', type=int, default=100,
            help='population size, default: %(default)s')
    allhic_group.add_argument(
            '--seed', type=int, default=42,
            help='random seed, default: %(default)s')

    # Parameters for performance
    performance_group = parser.add_argument_group('>>> Parameters for performance')
    performance_group.add_argument(
            '--processes', type=int, default=8,
            help='processes for fast sorting and ALLHiC optimization, always less than or equal to the number '
            'of input group files, default: %(default)s. Be aware that multiprocessing will increase RAM consumption')

    # Parameters for logging
    logging_group = parser.add_argument_group('>>> Parameters for logging')
    logging_group.add_argument(
            '--verbose', default=False, action='store_true',
            help='verbose logging, default: %(default)s')

    args = parser.parse_args()

    return args


def run(args, log_file=None):

    # (for pipeline) if log_file is provided, add an additional file handler for logging
    if log_file:
        file_handler = logging.FileHandler(log_file, 'w')
        formatter=logging.Formatter(
                fmt='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    start_time = time.time()
    logger.info('Program started, HapHiC version: {} (update: {})'.format(__version__, __update_time__))
    logger.info('Python version: {}'.format(sys.version.replace('\n', '')))
    logger.info('Command: {}'.format(' '.join(sys.argv)))

    # Check parameters
    if args.skip_fast_sort and args.skip_allhic:
        logger.error('Conflict parameters: --skip_fast_sort and --skip_allhic')
        raise RuntimeError('Conflict parameters')

    # quick view mode
    if args.quick_view:
        args.skip_allhic = True

    # processes are always less than or equal to the number of input group files
    processes = min(args.processes, len(args.groups))

    if processes > 1:
        p = Pool(processes)

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # check ALLHiC
    logger.info('Checking the path of ALLHiC...')
    script_realpath = os.path.dirname(os.path.realpath(__file__))
    allhic = os.path.join(script_realpath, 'allhic')

    if os.path.exists(allhic):
        logger.info('ALLHiC has been found in {}'.format(script_realpath))
    else:
        logger.error('CANNOT find ALLHiC in {}'.format(allhic))
        raise RuntimeError('CANNOT find ALLHiC')

    # parse fasta & pickle file

    # Using default RE here is ok. Because in the sorting step,
    # we don't care about the restriction sites.
    fa_dict = parse_fasta(args.fasta)

    # load HT_link_dict from pickle file
    with open(args.HT_links, 'rb') as f:
        logger.info('Loading input pickle file...')
        HT_link_dict = pickle.load(f)

    # extract group-specific data
    # a dict used to store group specific data, it can speed up multiprocessing and improve memory usage
    group_specific_data_dict = dict()

    # a dict used to store prefix and sub clm file for each group
    group_param_dict = dict()

    logger.info('Parsing group files and clm files...')

    for group in args.groups:
        # parse group files and clm files
        ctg_info_list, clm, prefix = parse_group(group, args.clm_dir, args.quick_view)
        group_param_dict[group] = (prefix, clm)
        ctgs = [ctg for ctg, ctg_len in ctg_info_list]
        sub_HT_dict, HT_index_dict = get_sub_HT_dict(ctgs, HT_link_dict)
        # In sub_HT_dict, the links between sister edges should not be counted.
        # HT_index_dict record the oldest version of contig HT names and their indexes, and will NOT be updated in iterations.
        group_specific_data_dict[group] = (ctg_info_list, ctgs, sub_HT_dict, HT_index_dict)

    del HT_link_dict
    gc.collect()

    if not args.skip_allhic:
        # for ALLHiC logs
        os.mkdir('allhic_logs')

    os.mkdir('final_tours')

    # with multiprocessing
    if processes > 1:

        logger.info('Program will be executed in multiprocessing mode (processes={})'.format(processes))

        # a list used to store the result for each processe (group)
        result_list = list()

        for group in args.groups:
            group_specific_data = group_specific_data_dict[group]
            group_specific_data_dict[group] = None
            result_list.append(p.apply_async(
                run_haphic_sorting, args=(args, group, fa_dict, group_specific_data, group_param_dict[group], allhic)))

        p.close()
        p.join()

        # check exceptions
        check_exceptions(result_list)

    # without multiprocessing
    else:
        for group in args.groups:
            group_specific_data = group_specific_data_dict[group]
            group_specific_data_dict[group] = None
            run_haphic_sorting(args, group, fa_dict, group_specific_data, group_param_dict[group], allhic)

    finished_time = time.time()
    logger.info('Program finished in {}s'.format(finished_time-start_time))


def main():

    # get arguments
    args = parse_arguments()

    run(args, log_file='HapHiC_sort.log')


if __name__ == '__main__':
    main()

