#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Email: zengxf@sustech.edu.cn

import sys
import argparse
import os
import logging
import time
import gc

from copy import deepcopy
import pysam
import pickle
from collections import defaultdict, OrderedDict
import random
from itertools import combinations
from array import array
from portion import closed, empty
import gzip

from math import ceil
from numpy import inf, int32, float32, zeros, quantile, arange, power, allclose, median
from numpy import abs as npabs
from numpy import array as ndarray
from numpy.linalg import matrix_power
from scipy.sparse import coo_matrix, dok_matrix, csc_matrix
from scipy.stats import mode
from scipy.optimize import linear_sum_assignment
from sklearn.preprocessing import normalize
from networkx import Graph, find_cliques, connected_components, shortest_path
from decimal import Decimal

from _version import __version__, __update_time__

try:
    from sparse_dot_mkl import dot_product_mkl
    INTEL_MKL = True
except:
    INTEL_MKL = False

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

pysam.set_verbosity(0)


def parse_RE_sites(sites):

    output_sites = list()

    for site in sites:
        if 'N' in site:
            output_sites.append(site.replace('N', 'A', 1))
            output_sites.append(site.replace('N', 'T', 1))
            output_sites.append(site.replace('N', 'C', 1))
            output_sites.append(site.replace('N', 'G', 1))
        else:
            output_sites.append(site)

    if 'N' not in ''.join(output_sites):
        return output_sites
    else:
        return parse_RE_sites(output_sites)


def count_RE_sites(seq, RE):

    sites = [site.strip().upper() for site in RE.split(',') if site.strip()]
    parsed_sites = parse_RE_sites(sites)

    RE_sites = 0
    for site in parsed_sites:
        RE_sites += seq.count(site)

    return RE_sites


def parse_fasta(fasta, RE='GATC', keep_letter_case=False, logger=logger):

    """save sequences, lengths, and RE site counts of contigs into a dict"""
    logger.info('Parsing input FASTA file...')

    fa_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                fa_dict[ctg] = list()
            else:
                if keep_letter_case:
                    fa_dict[ctg].append(line.strip())
                else:
                    fa_dict[ctg].append(line.strip().upper())

    for ctg, seq_list in fa_dict.items():
        # joining list is faster than concatenating strings
        seq = ''.join(seq_list)
        # add pseudo-count of 1 to prevent division by zero (as what ALLHiC does)
        RE_sites = count_RE_sites(seq, RE) + 1
        fa_dict[ctg] = [seq, len(seq), RE_sites]

    return fa_dict


def parse_gfa(gfa_list, fa_dict, logger=logger):

    """get the read depth and phasing information from gfa file(s)"""
    logger.info('Parsing input gfa file(s)...')

    read_depth_dict = dict()

    for n, gfa in enumerate(gfa_list):
        with open(gfa) as f:
            for line in f:
                if line.startswith('S\t'):
                    cols = line.split('\t')
                    ctg = cols[1]
                    ctg_len = int(cols[3].split(':')[-1])
                    read_depth = int(cols[4].split(':')[-1])
                    # contig length check
                    if ctg in fa_dict and ctg_len != fa_dict[ctg][1]:
                        logger.error('The contig {} in gfa file {} has a different length than the one in the fasta file. '
                                     'Maybe the gfa file does not match the fasta file.'.format(ctg, gfa))
                        raise RuntimeError('The contig {} in gfa file {} has a different length than the one in the fasta file. '
                                           'Maybe the gfa file(s) does not match the fasta file.'.format(ctg, gfa))
                    read_depth_dict[ctg] = (n, read_depth)

    # all contigs in fa_dict should be in read_depth_dict
    for ctg in fa_dict:
        if ctg not in read_depth_dict:
            logger.error('Can not find contig {} in the gfa file(s). Maybe the gfa file(s) does not match the fasta file.'.format(ctg))
            raise RuntimeError('Can not find contig {} in the gfa file(s). Maybe the gfa file(s) does not match the fasta file.'.format(ctg))

    # log a warning if gfa_seq_num > fa_seq_num
    gfa_seq_num, fa_seq_num = len(read_depth_dict), len(fa_dict)
    if gfa_seq_num > fa_seq_num:
        logger.warning('The number of contigs in the gfa file(s) ({}) is greater than that in the fasta file ({}). '
                       'Maybe some contigs were removed in the fasta file?'.format(gfa_seq_num, fa_seq_num))

    return read_depth_dict


def stat_fragments(fa_dict, RE, read_depth_dict, whitelist, nchrs=0, flank=0, Nx=100, bin_size=0, logger=logger):

    """basic statistics for fragments (contigs / bins)"""

    def count_flank_RE_sites(sequence, length):

        # all RE sites of each bin
        if not flank or length <= 2 * flank:
            return count_RE_sites(sequence, RE) + 1
        # RE sites in flanking regions of each bin
        else:
            return count_RE_sites(sequence[:flank], RE) + count_RE_sites(sequence[length-flank:], RE) + 1

    # we use the word "fragment" to represent contig or bin.
    logger.info('Making some statistics of fragments (contigs / bins)')

    # kbp -> bp
    flank *= 1000

    # calculate bin_size
    total_len = sum([ctg_info[1] for ctg_info in fa_dict.values()])

    # bin_size is 0 means no fragments will be split, so bin_size should be +infinity in actual
    if not bin_size:
        logger.info('bin_size is set to {}, no fragments will be split'.format(bin_size))
        bin_size = inf
    # default
    elif bin_size < 0:
        bin_size = max(min(int(total_len / nchrs / 30), 2000000), 100000)
        logger.info('bin_size is calculated to be {} bp'.format(bin_size))
    # manually designated
    else:
        bin_size *= 1000
        logger.info('bin_size is manually designated to {} bp'.format(bin_size))

    frags = list()
    bin_set = set()
    split_ctg_set = set()
    RE_site_dict, frag_len_dict = dict(), dict()

    for ctg, (seq, ctg_len, RE_sites) in fa_dict.items():
        # contigs are split into bins
        if ctg_len > bin_size:
            split_ctg_set.add(ctg)
            nbins = ceil(ctg_len / bin_size)
            for m in range(nbins):
                bin_ = '{}_bin{}'.format(ctg, m+1)
                assert bin_ not in fa_dict
                frags.append(bin_)
                bin_set.add(bin_)
                # not the last bin
                if m+1 < nbins:
                    bin_len = bin_size
                    bin_seq = seq[m*bin_size:(m+1)*bin_size]
                    RE_site_dict[bin_] = count_flank_RE_sites(bin_seq, bin_len)
                # the last bin
                else:
                    bin_len = ctg_len-m*bin_size
                    bin_seq = seq[m*bin_size:]
                    RE_site_dict[bin_] = count_flank_RE_sites(bin_seq, bin_len)
                frag_len_dict[bin_] = bin_len
                # update read_depth_dict
                if read_depth_dict:
                    read_depth_dict[bin_] = read_depth_dict[ctg]
            if read_depth_dict:
                del read_depth_dict[ctg]

        # contigs are not split
        else:
            frags.append(ctg)
            frag_len_dict[ctg] = ctg_len
            # all RE sites of each contig
            if not flank or ctg_len <= 2 * flank:
                RE_site_dict[ctg] = RE_sites
            # RE sites in flanking regions of each contig
            else:
                RE_site_dict[ctg] = count_flank_RE_sites(seq, ctg_len)

        # the sequences are unnecessary anymore, release memory
        fa_dict[ctg][0] = None

    # In the clustering step, if Nx < 100, longer fragments are used to represent the whole genome.
    # However, if many fragments have the same length, the Python built-in function 'sorted' will
    # keep the relative order of them and make the process biased (unable to represent the whole genome).
    # Shuffling the fragments before sorting can make the result more robust
    random.seed(12345)
    random.shuffle(frags)
    sorted_frag_list = sorted([(frag, frag_len_dict[frag]) for frag in frags], key=lambda x: x[1], reverse=True)

    # get total length of the draft genome
    len_sum = 0
    Nx_frag_set = set()

    for frag, frag_len in sorted_frag_list:
        len_sum += frag_len
        if len_sum/total_len*100 < Nx or Nx == 100:
            Nx_frag_set.add(frag)

    if Nx != 100:
        # add one more fragment to make len_sum/total*100 >= Nx
        Nx_frag_set.add(sorted_frag_list[len(Nx_frag_set)][0])

    # add fragments in whitelist
    if whitelist:
        for frag, _ in sorted_frag_list:
            if frag.rsplit('_bin', 1)[0] in whitelist:
                Nx_frag_set.add(frag)

    return sorted_frag_list, bin_set, bin_size, frag_len_dict, Nx_frag_set, RE_site_dict, split_ctg_set


def is_flank(coord, length, flank):

    """determine whether the coord is inside the flanking regions of a given fragment"""
    if flank and (coord <= flank or coord > length - flank):
        return True
    elif not flank:
        return True
    else:
        return False


def dict_to_matrix(link_dict, frag_set, dense_matrix=True, add_self_loops=False):

    """convert the Hi-C link dict to a adjacency matrix"""

    def add_data(x, y, c):
        row.append(x)
        col.append(y)
        data.append(c)

    # prepare three lists for dict to matrix conversion
    row, col, data = list(), list(), list()

    index = 0
    shape = len(frag_set)
    frag_index_dict = dict()
    frags_in_dict = set()

    for frag_name_pair, links in link_dict.items():

        frag_i, frag_j = frag_name_pair

        if frag_i not in frag_set or frag_j not in frag_set:
            continue

        frags_in_dict.add(frag_i)
        frags_in_dict.add(frag_j)

        if frag_i in frag_index_dict:
            i = frag_index_dict[frag_i]
        else:
            i = index
            frag_index_dict[frag_i] = i
            index += 1

        if frag_j in frag_index_dict:
            j = frag_index_dict[frag_j]
        else:
            j = index
            frag_index_dict[frag_j] = j
            index += 1

        add_data(i, j, links)
        # diagonal symmetry
        add_data(j, i, links)

    # record the index of remaining frags (no links but in frag_set)
    assert len(frags_in_dict) == index
    for frag in frag_set - frags_in_dict:
        frag_index_dict[frag] = index
        index += 1

    # add self loops (for Markov clustering)
    if add_self_loops:
        for n in range(shape):
            add_data(n, n, 1)

    if not dense_matrix:
        # from coo to csc (sparse matrix)
        matrix = coo_matrix((data, (row, col)), shape=(shape, shape), dtype=float32).tocsc()
    else:
        # from coo to ndarray (dense matrix)
        matrix = coo_matrix((data, (row, col)), shape=(shape, shape), dtype=float32).toarray()

    return matrix, frag_index_dict


def output_clm(clm_dict):

    logger.info('Writing clm_dict to paired_links.clm...')

    ori_tuple = (('+', '+'), ('+', '-'), ('-', '+'), ('-', '-'))

    with open('paired_links.clm', 'w') as fout:
        for ctg_name_pair, list_ in clm_dict.items():
            # minLinks == 3, links = len(list_)/2
            if len(list_) < 8:
                continue
            for n in range(4):
                new_list = ['{0} {0}'.format(v) for v in sorted(list_[n::4])]
                fout.write('{}{} {}{}\t{}\t{}\n'.format(
                    ctg_name_pair[0], ori_tuple[n][0],
                    ctg_name_pair[1], ori_tuple[n][1],
                    len(new_list)*2, ' '.join(new_list)))


def update_clm_dict(clm_dict, ctg_name_pair, len_i, len_j, coord_i_0, coord_j_0):

    clm_dict[ctg_name_pair].extend((
        len_i - coord_i_0 + coord_j_0,
        len_i - coord_i_0 + len_j - coord_j_0,
        coord_i_0 + coord_j_0,
        coord_i_0 + len_j - coord_j_0))


def update_HT_link_dict(HT_link_dict, ctg_i, ctg_j, len_i, len_j, coord_i, coord_j):

    def add_suffix(ctg, ctg_len, coord):

        if coord * 2 > ctg_len:
            return ctg + '_T'
        else:
            return ctg + '_H'

    ctg_name_HT_i = add_suffix(ctg_i, len_i, coord_i)
    ctg_name_HT_j = add_suffix(ctg_j, len_j, coord_j)

    HT_link_dict[(ctg_name_HT_i, ctg_name_HT_j)] += 1


def cal_concordance_ratio(coord_list, shorter_len, nwindows):

    bin_width = shorter_len // nwindows
    npairs = len(coord_list)//2
    # y = x + b
    y_minus_x_list = [(coord_list[2*n+1] - coord_list[2*n])//bin_width for n in range(npairs)]
    # y = -x + b
    y_plus_x_list = [(coord_list[2*n+1] + coord_list[2*n])//bin_width for n in range(npairs)]

    return max(mode(y_minus_x_list, keepdims=False)[1] / npairs, mode(y_plus_x_list, keepdims=False)[1] / npairs)


def cal_concentration_adj_ratio(coord_list, bin_width=10000):

    npairs = len(coord_list)//2

    x_bin_dict = defaultdict(int)
    y_bin_dict = defaultdict(int)

    for n in range(npairs):
        x_bin_dict[coord_list[2*n]//bin_width] += 1
        y_bin_dict[coord_list[2*n+1]//bin_width] += 1

    x_bin_list = x_bin_dict.values()
    y_bin_list = y_bin_dict.values()

    x_bin_median = median([links for links in x_bin_list if links])
    y_bin_median = median([links for links in y_bin_list if links])

    concentration_ratio_x = sum([links for links in x_bin_list if links >= 10 * x_bin_median])/npairs
    concentration_ratio_y = sum([links for links in y_bin_list if links >= 10 * y_bin_median])/npairs

    return (1-concentration_ratio_x) * (1-concentration_ratio_y)


def record_coord_pairs(ctg_coord_dict, ctg_name_pair, coord_i, coord_j, max_read_pairs, fa_dict, args):

    # record coord pairs and calculate concordance ratio once there are enough coord pairs
    if not isinstance(ctg_coord_dict[ctg_name_pair], list):

        ctg_coord_dict[ctg_name_pair].extend((coord_i, coord_j))

        if len(ctg_coord_dict[ctg_name_pair]) >= max_read_pairs * 2:
            if args.remove_allelic_links:
                shorter_len = min(fa_dict[ctg_name_pair[0]][1], fa_dict[ctg_name_pair[1]][1])
                concordance_ratio = cal_concordance_ratio(ctg_coord_dict[ctg_name_pair], shorter_len, args.nwindows)
                ctg_coord_dict[ctg_name_pair] = [concordance_ratio, 1]
            if args.remove_concentrated_links:
                adj_ratio = cal_concentration_adj_ratio(ctg_coord_dict[ctg_name_pair])
                if args.remove_allelic_links:
                    ctg_coord_dict[ctg_name_pair][1] = adj_ratio
                else:
                    ctg_coord_dict[ctg_name_pair] = [0, adj_ratio]


def remove_allelic_HiC_links(fa_dict, ctg_coord_dict, full_link_dict, args, flank_link_dict=None, filtered_frags=None, ctg_pair_to_frag=None, logger=logger):

    logger.info("Removing Hi-C links between alleic contig pairs...")

    # Similar to ALLHiC, here we remove two kinds of allelic Hi-C links:

    # 1) the Hi-C links between allelic contigs which are identified by "concordance ratio"

    # Then, we use these inter-allele info to construct a matrix, and identify allele groups
    # by searching and split cliques (networkx). Finally, we find best contig-contig matches
    # across allele groups (Hungarian algorithm, maximum bipartite matching problem), and remove:

    # 2) the Hi-C links between non-max matches

    def update_link_dicts(ctg_name_pair, type_):

        if type_ == 1:
            # record inter-allele info
            inter_allele_dict[ctg_name_pair] = full_link_dict[ctg_name_pair]
            allelic_ctg_set.add(ctg_name_pair[0])
            allelic_ctg_set.add(ctg_name_pair[1])

        # update full_link_dict
        del full_link_dict[ctg_name_pair]

        if flank_link_dict:
            # for bins
            if ctg_pair_to_frag:
                for frag_name_pair in ctg_pair_to_frag[ctg_name_pair]:
                    if frag_name_pair in flank_link_dict and frag_name_pair[0] in filtered_frags and frag_name_pair[1] in filtered_frags:
                        # update flank_link_dict
                        del flank_link_dict[frag_name_pair]
            # for contigs
            elif ctg_name_pair in flank_link_dict and ctg_name_pair[0] in filtered_frags and ctg_name_pair[1] in filtered_frags:
                # update flank_link_dict
                del flank_link_dict[ctg_name_pair]

    def get_weakest_edge(graph):

        weakest_edge = (None, None, inf)

        for node1, node2, data in graph.edges(data=True):
            if node1 == node2:
                continue
            if data['weight'] < weakest_edge[-1]:
                weakest_edge = (node1, node2, data['weight'])

        assert weakest_edge[0] is not None

        return weakest_edge

    def split_cliques(graph, cliques, ploidy, cached_cliques):

        new_cliques = set()

        for clique in cliques:

            clique = tuple(clique)

            if len(clique) > ploidy:
                if clique not in cached_cliques:
                    subgraph = graph.subgraph(clique)
                    node1, node2, _ = get_weakest_edge(subgraph)
                    # unfreeze graph
                    subgraph = Graph(subgraph)
                    # remove weakest edges
                    subgraph.remove_edge(node1, node2)
                    # find sub cliques after removing weakest edges
                    sub_cliques = find_cliques(subgraph)
                    # cache sub_cliques
                    cached_cliques.add(clique)
                    # call split_cliques recursively
                    new_cliques |= split_cliques(subgraph, sub_cliques, ploidy, cached_cliques)
            else:
                new_cliques.add(tuple(clique))

        return new_cliques

    def solve_max_matching(group_pair):

        group_1, group_2 = group_pair

        # the input matrix of Hungarian algorithm should be n*n
        degree = max(len(group_1), len(group_2))
        matching_matrix = zeros((degree, degree), dtype=int)

        for i1, c1 in enumerate(group_1):
            for i2, c2 in enumerate(group_2):
                ctg_name_pair = tuple(sorted([c1, c2]))
                if ctg_name_pair in full_link_dict:
                    matching_matrix[i1, i2] = full_link_dict[ctg_name_pair]

        # Hungarian algorithm is designed to find minimum cost, here we use
        # negative weight to convert maximum matching problem to minimum cost problem
        return linear_sum_assignment(-matching_matrix)

    ploidy = args.remove_allelic_links
    min_read_pairs = args.min_read_pairs
    concordance_ratio_cutoff = args.concordance_ratio_cutoff

    # a dict storing the inter-allele info
    inter_allele_dict = dict()
    allelic_ctg_set = set()

    # 1) the Hi-C links between allelic contigs which are identified by "concordance ratio"
    for ctg_name_pair, data in ctg_coord_dict.items():

        if isinstance(data, list):
            logger.debug('{} {} links={} concordance_ratio={}'.format(
                *ctg_name_pair, full_link_dict[ctg_name_pair], data[0]))

            if data[0] > concordance_ratio_cutoff:
                update_link_dicts(ctg_name_pair, 1)

        else:
            if len(data) >= min_read_pairs * 2:
                shorter_len = min(fa_dict[ctg_name_pair[0]][1], fa_dict[ctg_name_pair[1]][1])
                concordance_ratio = cal_concordance_ratio(data, shorter_len, args.nwindows)
                logger.debug('{} {} links={} concordance_ratio={}'.format(
                    *ctg_name_pair, full_link_dict[ctg_name_pair], concordance_ratio))

                if concordance_ratio > concordance_ratio_cutoff:
                    update_link_dicts(ctg_name_pair, 1)
            else:
                logger.debug('{} {} links={} concordance_ratio=0'.format(
                    *ctg_name_pair, full_link_dict[ctg_name_pair]))

    # only when ploidy > 2, finding & splitting cliques are necessary
    if ploidy > 2:
        allele_matrix, ctg_index_dict = dict_to_matrix(inter_allele_dict, allelic_ctg_set)
        index_ctg_dict = {i: c for c, i in ctg_index_dict.items()}

        allele_graph = Graph(allele_matrix)
        allele_groups = find_cliques(allele_graph)
        cached_cliques = set()
        allele_groups = split_cliques(allele_graph, allele_groups, ploidy, cached_cliques)
        del cached_cliques
        gc.collect()
        # use a set to remove redundant allele groups
        unique_allele_groups = set()

        for group in allele_groups:
            unique_allele_groups.add(tuple(sorted([index_ctg_dict[i] for i in group])))
    # when ploidy == 2, ctg-ctg pairs are equal to allele groups
    else:
        unique_allele_groups = set(inter_allele_dict.keys())

    # 2) the Hi-C links between non-max matches

    # get ctg_allele_group_dict
    ctg_allele_group_dict = defaultdict(set)
    for group in unique_allele_groups:
        for ctg in group:
            ctg_allele_group_dict[ctg].add(group)

    # cache solutions
    solution_dict = dict()
    # store non-max ctg pair
    nonmax_ctg_pair_set = set()

    for ctg_name_pair in full_link_dict:

        # allelic contig pairs have been already removed by 'del'
        # if one of the ctgs is not in any allele group, do nothing
        ctg_1, ctg_2 = ctg_name_pair
        if ctg_1 not in ctg_allele_group_dict or ctg_2 not in ctg_allele_group_dict:
            continue

        for group_1 in ctg_allele_group_dict[ctg_1]:
            for group_2 in ctg_allele_group_dict[ctg_2]:
                group_pair = tuple(sorted([group_1, group_2]))

                if group_pair in solution_dict:
                    solution = solution_dict[group_pair]
                else:
                    solution = solve_max_matching(group_pair)
                    solution_dict[group_pair] = solution

                if ctg_1 in group_pair[0]:
                    assert ctg_2 in group_pair[1]
                    index_1 = group_pair[0].index(ctg_1)
                    index_2 = group_pair[1].index(ctg_2)
                else:
                    assert ctg_2 in group_pair[0]
                    index_1 = group_pair[0].index(ctg_2)
                    index_2 = group_pair[1].index(ctg_1)

                if solution[1][index_1] != index_2:
                    nonmax_ctg_pair_set.add(ctg_name_pair)
                    break
            # break nested loops
            else:
                continue
            break

    # remove links between non-max ctg pair
    for ctg_name_pair in nonmax_ctg_pair_set:

        logger.debug('{} {} links={} non-maximum matching'.format(
            *ctg_name_pair, full_link_dict[ctg_name_pair]))

        update_link_dicts(ctg_name_pair, 2)

    # remove isolated fragments (maybe unecessary)
    if flank_link_dict:
        remaining_frags =  set()
        for frag_1, frag_2 in flank_link_dict:
            if frag_1 in filtered_frags and frag_2 in filtered_frags:
                remaining_frags.add(frag_1)
                remaining_frags.add(frag_2)

        removed_frags = filtered_frags - remaining_frags

        logger.info('Removing isolated fragments after filtering out allelic Hi-C links...')
        logger.info('{} fragments removed, {} fragments kept'.format(len(removed_frags), len(remaining_frags)))
        for frag in removed_frags:
            logger.debug('Fragment {} is isolated and removed'.format(frag))

        return remaining_frags


def reduce_inter_hap_HiC_links(link_dict, read_depth_dict, phasing_weight, target='flank_link_dict'):

    logger.info('Reducing inter-haplotype Hi-C links in {}...'.format(target))

    removed_frag_name_pairs = set()
    for frag_name_pair, _ in link_dict.items():
        if read_depth_dict[frag_name_pair[0]][0] != read_depth_dict[frag_name_pair[1]][0]:
            link_dict[frag_name_pair] -= link_dict[frag_name_pair] * phasing_weight
            if link_dict[frag_name_pair] == 0:
                removed_frag_name_pairs.add(frag_name_pair)

    for frag_name_pair in removed_frag_name_pairs:
        del link_dict[frag_name_pair]


def output_pickle(dict_, from_, to):

    logger.info('Writing {} to {}...'.format(from_, to))

    with open(to, 'wb') as fpkl:
        pickle.dump(dict_, fpkl)


def normalize_by_nlinks(flank_link_dict, frag_link_dict):

    logger.info('Normalizing flank_link_dict by the number of links to other contigs...')

    for ctg_i, ctg_j in flank_link_dict:
        # normalized by the geometric average of links of two contigs
        flank_link_dict[(ctg_i, ctg_j)] /= (frag_link_dict[ctg_i] * frag_link_dict[ctg_j]) ** 0.5


def normalize_by_length(flank_link_dict, frag_len_dict, flank):

    logger.info('Normalizing flank_link_dict by length...')

    two_flanks = flank * 2000

    for ctg_i, ctg_j in flank_link_dict:
        len_i, len_j = frag_len_dict[ctg_i], frag_len_dict[ctg_j]
        flank_len_i = len_i if len_i <= two_flanks else two_flanks
        flank_len_j = len_j if len_j <= two_flanks else two_flanks
        # length multiplication
        flank_link_dict[(ctg_i, ctg_j)] /= (flank_len_i/1000000) * (flank_len_j/1000000)


def filter_fragments(
        Nx_frag_set, RE_site_dict, RE_site_cutoff, frag_link_dict, density_lower, density_upper,
        topN, rank_sum_upper, rank_sum_hard_cutoff, flank_link_dict, read_depth_dict, read_depth_upper, whitelist):

    logger.info('Filtering fragments...')

    frags_in_whitelist = set()
    frag_density_list = list()

    total_links = 0
    total_RE_sites = 1
    # (1) Nx filtering
    for frag in Nx_frag_set:
        RE_sites = RE_site_dict[frag]
        # (2) RE sites >= RE_site_cutoff (RE sites have already + 1)
        if RE_sites > RE_site_cutoff:
            if frag in frag_link_dict:
                frag_links = frag_link_dict[frag]
                total_links += frag_links
                total_RE_sites += RE_sites - 1
                frag_density_list.append((frag, frag_links/RE_sites))
            else:
                frag_density_list.append((frag, 0))
        if whitelist and frag.rsplit('_bin', 1)[0] in whitelist:
            frags_in_whitelist.add(frag)

    Nx_frag_num = len(Nx_frag_set)

    logger.info('[Nx filtering] {} fragments kept'.format(Nx_frag_num))
    logger.info('[RE sites filtering] {} fragments removed, {} fragments kept'.format(
        Nx_frag_num-len(frag_density_list), len(frag_density_list)))

    # (3) fragment link density rank ranges from density_lower ~ density_upper
    frag_density_list.sort(key=lambda x: x[1])

    param_density_lower = check_param('--density_lower', density_lower, {'X', 'x'})
    param_density_upper = check_param('--density_upper', density_upper, {'X', 'x'})
    remaining_nfrags = len(frag_density_list)
    average_density = total_links/total_RE_sites

    if param_density_lower[-1] in {'X', 'x'}:
        for lower, (frag, density) in enumerate(frag_density_list):
            if density >= average_density * param_density_lower[0]:
                break
        else:
            lower += 1
        logger.info('[link density filtering] Parameter --density_lower {} is set to "multiple" mode and equivalent to {} in "fraction" mode'.format(
            density_lower, lower/remaining_nfrags))
    else:
        lower = int(remaining_nfrags * float(density_lower))
        logger.info('[link density filtering] Parameter --density_lower {} is set to "fraction" mode and equivalent to {}X in "multiple" mode'.format(
            density_lower, frag_density_list[max(0, lower-1)][1]/average_density))

    if param_density_upper[-1] in {'X', 'x'}:
        for upper, (frag, density) in enumerate(frag_density_list):
            if density > average_density * param_density_upper[0]:
                break
        else:
            upper += 1
        logger.info('[link density filtering] Parameter --density_upper {} is set to "multiple" mode and equivalent to {} in "fraction" mode'.format(
            density_upper, upper/remaining_nfrags))
    else:
        upper = int(remaining_nfrags * float(density_upper))
        logger.info('[link density filtering] Parameter --density_upper {} is set to "fraction" mode and equivalent to {}X in "multiple" mode'.format(
            density_upper, frag_density_list[max(0, upper-1)][1]/average_density))

    # get link density filtered fragments
    filtered_frags = {frag for frag, density in frag_density_list[lower:upper]}
    nfiltered_frags = len(filtered_frags)
    logger.info('[link density filtering] {} fragments removed, {} fragments kept'.format(
        remaining_nfrags-nfiltered_frags, nfiltered_frags))

    for frag, density in frag_density_list[:lower] + frag_density_list[upper:]:
        logger.debug('[link density filtering] Fragment {} is removed, density={}'.format(frag, density))

    frag_density_list_unfiltered = frag_density_list
    frag_density_list = frag_density_list[lower:upper]

    # (4) read depth filtering
    # sort frag_density_list by read depth
    if read_depth_dict:
        read_depth_list = [(frag, read_depth_dict[frag][1]) for frag, density in frag_density_list_unfiltered]
        read_depth_list.sort(key=lambda x: x[1])

        param_read_depth_upper = check_param('--read_depth_upper', read_depth_upper, {'X', 'x'})

        # use IQR to filter out outliers, values > Q3 + 1.5*(Q3-Q1) are regard as outliers
        q1, m, q3 = quantile([read_depth for frag, read_depth in read_depth_list], (0.25, 0.5, 0.75))
        iqr = q3 - q1
        logger.info('[read depth filtering] Q1={}, median={}, Q3={}, IQR=Q3-Q1={}'.format(q1, m, q3, iqr))

        if param_read_depth_upper[-1]:
            read_depth_upper_limit = q3 + param_read_depth_upper[0] * iqr
            for upper, (frag, read_depth) in enumerate(read_depth_list):
                if read_depth > read_depth_upper_limit:
                    break
            else:
                upper += 1
            logger.info('[read depth filtering] Parameter --read_depth_upper {} is set to "multiple" mode and equivalent to {} in "fraction" mode'.format(
                read_depth_upper, upper/remaining_nfrags))
        else:
            upper = int(remaining_nfrags * float(read_depth_upper))
            logger.info('[read depth filtering] Parameter --read_depth_upper {} is set to "fraction" mode and equivalent to {}X in "multiple" mode'.format(
                 read_depth_upper, (read_depth_list[max(0, upper-1)][1]-q3)/iqr))

        # get read depth filtered fragments
        filtered_frags &= {frag for frag, read_depth in read_depth_list[:upper]}
        nfiltered_frags = len(filtered_frags)

        read_depth_removed_frags = {frag for frag, read_depth in read_depth_list[upper:]}
        link_density_removed_frags = {frag for frag, density in frag_density_list_unfiltered[:lower] + frag_density_list_unfiltered[upper:]}
        specific_read_depth_removed_frags = read_depth_removed_frags - link_density_removed_frags

        logger.info('[read depth filtering] {} fragments removed, {} fragments kept'.format(
            len(specific_read_depth_removed_frags), nfiltered_frags))

        for frag, read_depth in read_depth_list[upper:]:
            if frag in specific_read_depth_removed_frags:
                logger.debug('[read depth filtering] Fragment {} is removed, read depth={}'.format(frag, read_depth))

        # update frag_density_list for the following filtering step
        frag_density_list = [(frag, density) for frag, density in frag_density_list if frag in filtered_frags]

    # (5) for each fragment, get fragments with top 10 link density, then
    #     calculate the sum of link ranks between these fragments

    # convert link dict to matrix (for statistics)
    flank_link_matrix, frag_index_dict = dict_to_matrix(flank_link_dict, filtered_frags)
    index_frag_dict = {i: f for f, i in frag_index_dict.items()}

    # construct a dict storing sorted frags by number of links
    frag_rank_dict = dict()

    for frag, _ in frag_density_list:
        index = frag_index_dict[frag]
        link_list = [(i, links) for i, links in enumerate(flank_link_matrix[index, :])]
        link_list.sort(key=lambda x: x[1], reverse=True)
        frag_rank_dict[frag] = [index_frag_dict[i] for i, _ in link_list]

    # calculate the sum of link ranks between topN nearest fragments
    rank_sum_list = list()
    hard_filtered_nfrags = 0
    for frag, _ in frag_density_list:
        top_N_frags = frag_rank_dict[frag][:topN]
        rank_sum = 0
        for frag_1, frag_2 in combinations(top_N_frags, 2):
            rank_sum += min(frag_rank_dict[frag_1].index(frag_2), frag_rank_dict[frag_2].index(frag_1))
        if rank_sum_hard_cutoff and rank_sum > rank_sum_hard_cutoff:
            hard_filtered_nfrags += 1
            logger.debug('[rank sum filtering] Fragment {} is removed by hard filtering, rank sum={}'.format(frag, rank_sum))
            continue
        rank_sum_list.append((frag, rank_sum))

    rank_sum_list.sort(key=lambda x: x[1])
    remaining_nfrags = len(rank_sum_list)
    if rank_sum_hard_cutoff:
        logger.info('[rank sum filtering] {} fragments removed by hard filtering, {} fragments kept'.format(hard_filtered_nfrags, remaining_nfrags))

    param_rank_sum_upper = check_param('--rank_sum_upper', rank_sum_upper, {'X', 'x'})

    # use IQR to filter out outliers, values > Q3 + 1.5*(Q3-Q1) are regard as outliers
    q1, m, q3 = quantile([rank_sum for _, rank_sum in rank_sum_list], (0.25, 0.5, 0.75))
    iqr = q3 - q1
    logger.info('[rank sum filtering] Q1={}, median={}, Q3={}, IQR=Q3-Q1={}'.format(q1, m, q3, iqr))

    if param_rank_sum_upper[-1]:
        rank_sum_upper_limit = q3 + param_rank_sum_upper[0] * iqr
        for upper, (frag, rank_sum) in enumerate(rank_sum_list):
            if rank_sum > rank_sum_upper_limit:
                break
        else:
            upper += 1
        logger.info('[rank sum filtering] Parameter --rank_sum_upper {} is set to "multiple" mode and equivalent to {} in "fraction" mode'.format(
            rank_sum_upper, upper/remaining_nfrags))
    else:
        upper = int(remaining_nfrags * float(rank_sum_upper))
        logger.info('[rank sum filtering] Parameter --rank_sum_upper {} is set to "fraction" mode and equivalent to {}X in "multiple" mode'.format(
            rank_sum_upper, (rank_sum_list[max(0, upper-1)][1]-q3)/iqr))

    # get rank sum filtered fragments
    filtered_frags = {frag for frag, rank_sum in rank_sum_list[:upper]}
    nfiltered_frags = len(filtered_frags)
    logger.info('[rank sum filtering] {} fragments removed, {} fragments kept'.format(
        len(rank_sum_list)-nfiltered_frags, nfiltered_frags))

    for frag, rank_sum in rank_sum_list[upper:]:
        logger.debug('[rank sum filtering] Fragment {} is removed, rank sum={}'.format(frag, rank_sum))

    if frags_in_whitelist:
        n = 0
        for frag in frags_in_whitelist:
            if frag in filtered_frags:
                continue
            n += 1
            logger.debug('[rank sum filtering] Fragment {} is added since it is on the whitelist'.format(frag))
            filtered_frags.add(frag)
        logger.info('[rank sum filtering] {} fragments added, {} fragments are used to perform Markov clustering'.format(
            n, len(filtered_frags)))

    return filtered_frags


def detect_break_points(ctg_cov_dict, fa_dict, args):

    resolution = args.correct_resolution

    # use a dict to store the breakpoints of each contig
    ctg_break_point_dict = dict()

    for ctg, cov_list in ctg_cov_dict.items():

        median_cov = median(cov_list)
        # if median is 0, skip this contig
        if not median_cov:
            continue

        high_cov_regions = empty()
        filtered_high_cov_regions = empty()
        cov_cutoff = median_cov * args.median_cov_ratio
        region_cutoff = max(args.min_region_cutoff, fa_dict[ctg][1]*args.region_len_ratio)

        # find regions that covs >= cov_cutoff (high-coverage regions)
        for n, cov in enumerate(cov_list):
            if cov >= cov_cutoff:
                high_cov_regions |= closed(n*resolution, (n+1)*resolution)

        # consider only the low-coverage regions bounded by two high-coverage regions
        if len(high_cov_regions) < 2:
            continue

        # the high-coverage regions should be longer than region_cutoff
        for region in high_cov_regions:
            if region.upper-region.lower >= region_cutoff:
                filtered_high_cov_regions |= region

        # consider only the low-coverage regions bounded by two high-coverage regions
        if len(filtered_high_cov_regions) < 2:
            continue

        # find low-coverage regions (valleys) bounded by filtered high-coverage regions
        valley_intervals = closed(filtered_high_cov_regions.lower, filtered_high_cov_regions.upper) - filtered_high_cov_regions
        any_zero = False

        # get candidate breakpoints
        candidate_break_points = list()
        for interval in valley_intervals:

            # get the coverage of each bin in the valley
            start_bin, end_bin = interval.lower//resolution, interval.upper//resolution
            valley_bin_cov = cov_list[start_bin: end_bin]

            # if a valley contains one or more 0-coverage bins, there will be one (and only one) breakpoint
            if 0 in valley_bin_cov:
                any_zero = True
                # we arbitrarily regard the leftmost 0-coverage bin as the breakpoint
                candidate_break_points.append((valley_bin_cov.tolist().index(0)+start_bin, 0))
            else:
                # get the minimum-coverage (deepest) bin of each valley, as a candidate breakpoint
                min_valley_bin = valley_bin_cov.argmin()
                candidate_break_points.append((min_valley_bin+start_bin, valley_bin_cov[min_valley_bin]))

        # get final breakpoints
        if any_zero:
            # break 0-coverage bins
            ctg_break_point_dict[ctg] = [(site[0]*resolution, 0) for site in candidate_break_points if site[1] == 0]
        else:
            # break the deepest valley
            deepest_bin, deepest_cov = sorted(candidate_break_points, key=lambda x: x[1])[0]
            ctg_break_point_dict[ctg] = [(deepest_bin*resolution, deepest_cov)]

        # debug
        # print('stat: {}\t{}\t{}\t{}\t{}\t{}'.format(ctg, median_cov, cov_cutoff, valley_intervals, candidate_break_points, ctg_break_point_dict[ctg]))

    return ctg_break_point_dict


def break_and_update_ctgs(
        ctg_break_point_dict, ctg_link_pos_dict, ctg_cov_dict,
        frag_source_dict, final_break_pos_dict, final_break_frag_dict,
        fa_dict, read_depth_dict, unbroken_ctgs, args, last_round=False):

    # these dicts should be updated:
    # (1) fa_dict -> updated in all rounds
    # (2) ctg_cov_dict -> not necessary in the last round
    # (3) ctg_link_pos_dict -> not necessary in the last round
    # (4) frag_source_dict, final_break_pos_dict, final_break_frag_dict -> updated in all rounds

    def update_fa_dict(ctg, new_frag, start, end):

        new_seq = fa_dict[ctg][0][start:end]
        RE_sites = count_RE_sites(new_seq, args.RE)
        fa_dict[new_frag] = [new_seq, end-start, RE_sites]

        return end

    def pos_shift_for_multi_break_points(ctg, coord, pos_shift_list):

        if ctg not in unbroken_ctgs:
            assert ':' in ctg
            ctg = ctg.rsplit(':', 1)[0]

        for n, p in enumerate(pos_shift_list):
            new_coord = coord - p
            if new_coord >= 0:
                if n:
                    end = pos_shift_list[n-1]
                else:
                    end = fa_dict[ctg][1]
                return '{}:{}-{}'.format(ctg, p+1, end), new_coord

    def pos_shift_for_single_break_point(ctg, coord, pos_shift_list):

        if ctg in unbroken_ctgs:
            start = 1
            end = fa_dict[ctg][1]
        else:
            assert ':' in ctg
            ctg, pos_range = ctg.rsplit(':', 1)
            start, end = [int(p) for p in pos_range.split('-')]

        if coord >= pos_shift_list[0]:
            return '{}:{}-{}'.format(ctg, start+pos_shift_list[0], end), coord-pos_shift_list[0]
        else:
            return '{}:{}-{}'.format(ctg, start, pos_shift_list[0]), coord

    logger.info('Breaking contigs and updating data...')

    resolution = args.correct_resolution

    # to get unbroken contigs
    ctgs_before_breaking = set(fa_dict.keys())
    broken_ctgs = set()

    for ctg, break_points in ctg_break_point_dict.items():

        if not last_round:
            broken_ctgs.add(ctg)

            if break_points[0][1] == 0:
                # if cov of one breakpoint is zero, all will be zero
                # no need to remove coverages spanning breakpoints from ctg_cov_dict
                any_zero = True
            else:
                any_zero = False
                # if no zero-coverage breakpoints, there must be only one breakpoint
                # get the bin interval of the breakpoint
                break_point_interval = closed(break_points[0][0], break_points[0][0]+resolution)

            # a list used to calculate coordinate shift of Hi-C links after breaking contig
            pos_shift_list = [point for point, cv in break_points][::-1] + [0]

            # single breakpoint
            if len(pos_shift_list) == 2:
                pos_shift = pos_shift_for_single_break_point
            # multiple breakpoints
            else:
                pos_shift = pos_shift_for_multi_break_points

            npairs = len(ctg_link_pos_dict[ctg])//2
            # if cov of the breakpoint is not zero
            if not any_zero:
                for n in range(npairs):
                    coord_i, coord_j = ctg_link_pos_dict[ctg][2*n:2*n+2]
                    spanning_region = closed(coord_i, coord_j)

                    # if a Hi-C link spanning a breakpoint
                    if spanning_region.overlaps(break_point_interval):
                        # update ctg_cov_dict
                        spanning_start_bin, spanning_end_bin = coord_i//resolution, coord_j//resolution
                        ctg_cov_dict[ctg][spanning_start_bin:spanning_end_bin+1] -= 1
                    else:
                        # break ctg_link_pos_dict
                        new_frag_i, new_coord_i = pos_shift(ctg, coord_i, pos_shift_list)
                        new_frag_j, new_coord_j = pos_shift(ctg, coord_j, pos_shift_list)
                        # consider only intra-frag Hi-C links
                        if new_frag_i == new_frag_j:
                            ctg_link_pos_dict[new_frag_i].extend((new_coord_i, new_coord_j))
            # if cov of the breakpoint is zero
            else:
                # no need to get spanning_region
                for n in range(npairs):
                    coord_i, coord_j = ctg_link_pos_dict[ctg][2*n:2*n+2]

                    # break ctg_link_pos_dict
                    new_frag_i, new_coord_i = pos_shift(ctg, coord_i, pos_shift_list)
                    new_frag_j, new_coord_j = pos_shift(ctg, coord_j, pos_shift_list)
                    # consider only intra-frag Hi-C links
                    if new_frag_i == new_frag_j:
                        ctg_link_pos_dict[new_frag_i].extend((new_coord_i, new_coord_j))

        start = 0
        frag_source = frag_source_dict[ctg]
        father_index = final_break_frag_dict[frag_source].index(ctg)
        father_pos = final_break_pos_dict[frag_source][father_index]
        final_break_frag_dict[frag_source].pop(father_index)
        final_break_pos_dict[frag_source].pop(father_index)

        for n, (point, _) in enumerate(break_points, 1):

            # the first breakpoint
            if n == 1:
                # start: 1
                s = 1
                last_point = point
            else:
                # start: last breakpoint + 1
                s = last_point + 1
                last_point = point

            if ctg not in unbroken_ctgs:
                assert ':' in ctg
                raw_ctg, pos_range = ctg.rsplit(':', 1)
                shift = int(pos_range.split('-')[0]) - 1
                new_frag = '{}:{}-{}'.format(raw_ctg, s+shift, point+shift)
            else:
                new_frag = '{}:{}-{}'.format(ctg, s, point)

            # update frag_source_dict, final_break_frag_dict and final_break_pos_dict
            frag_source_dict[new_frag] = frag_source
            final_break_frag_dict[frag_source].insert(father_index, new_frag)
            final_break_pos_dict[frag_source].insert(father_index, father_pos+start)

            # update read_depth_dict
            if read_depth_dict:
                read_depth_dict[new_frag] = read_depth_dict[ctg]

            if not last_round:
                # break ctg_cov_dict
                ctg_cov_dict[new_frag] = ctg_cov_dict[ctg][start//resolution: point//resolution]

            # break & update fa_dict
            start = update_fa_dict(ctg, new_frag, start, point)

        # the last frag of a contig
        if ctg not in unbroken_ctgs:
            assert ':' in ctg
            raw_ctg, pos_range = ctg.rsplit(':', 1)
            shift = int(pos_range.split('-')[0]) - 1
            new_frag = '{}:{}-{}'.format(raw_ctg, shift+last_point+1, shift+fa_dict[ctg][1])
        else:
            new_frag = '{}:{}-{}'.format(ctg, last_point+1, fa_dict[ctg][1])

        # update frag_source_dict, final_break_frag_dict and final_break_pos_dict
        frag_source_dict[new_frag] = frag_source
        final_break_frag_dict[frag_source].insert(father_index, new_frag)
        final_break_pos_dict[frag_source].insert(father_index, father_pos+start)

        # update read_depth_dict
        if read_depth_dict:
            read_depth_dict[new_frag] = read_depth_dict[ctg]

        if not last_round:
            # break ctg_cov_dict
            ctg_cov_dict[new_frag] = ctg_cov_dict[ctg][start//resolution:]

        # break & update fa_dict
        ctg_len = fa_dict[ctg][1]
        update_fa_dict(ctg, new_frag, start, ctg_len)

        # remove original contig
        del fa_dict[ctg]
        if read_depth_dict:
            del read_depth_dict[ctg]
        if not last_round:
            del ctg_cov_dict[ctg]
            # del ctg_link_pos_dict[ctg]

    if not last_round:
        # if ctg has no breakpoint, remove it from ctg_cov_dict
        # these ctgs will not be checked again in the next round
        for ctg in ctgs_before_breaking - broken_ctgs:
            if ctg in ctg_cov_dict:
                del ctg_cov_dict[ctg]


def correct_assembly(ctg_cov_dict, ctg_link_pos_dict, fa_dict, read_depth_dict, args):

    logger.info('Performing assembly correction...')

    unbroken_ctgs = set(fa_dict.keys())

    frag_source_dict = dict()
    final_break_pos_dict = dict()
    final_break_frag_dict = dict()

    for nround in range(args.correct_nrounds):

        # find breakpoints in contigs
        ctg_break_point_dict = detect_break_points(ctg_cov_dict, fa_dict, args)
        logger.info('Correction round {}, breakpoints are detected in {} contig(s)'.format(
            nround+1, len(ctg_break_point_dict)))

        if nround == 0:
            nbroken_ctgs = len(ctg_break_point_dict)

        # if no breakpoint found, break the iteration
        if not ctg_break_point_dict:
            break

        # the first round, initialization
        if nround == 0:
            for ctg in ctg_break_point_dict:
                frag_source_dict[ctg] = ctg
                final_break_pos_dict[ctg] = [0]
                final_break_frag_dict[ctg] = [ctg]

        # break contigs and update info
        # if it is the last round, only fa_dict will be updated to save time
        if nround + 1 == args.correct_nrounds:
            last_round = True
        else:
            last_round = False

        break_and_update_ctgs(
                ctg_break_point_dict, ctg_link_pos_dict, ctg_cov_dict,
                frag_source_dict, final_break_pos_dict, final_break_frag_dict,
                fa_dict, read_depth_dict, unbroken_ctgs, args, last_round)

        unbroken_ctgs -= set(ctg_break_point_dict.keys())
        # print(nround)
        # for frag, source in frag_source_dict.items():
        #     print('source:', frag, source)
        # for ctg, pos_list in final_break_pos_dict.items():
        #     print('pos:', ctg, pos_list)
        # for ctg, frag_list in final_break_frag_dict.items():
        #     print('frag:', ctg, frag_list)

    # output corrected assembly in FASTA format, even if no contigs are corrected (for pipeline)
    corrected_assembly_file = 'corrected_asm.fa'
    corrected_ctgs_file = 'corrected_ctgs.txt'
    logger.info('Generating corrected assembly file...')

    # check file name (rare case, for security concern)
    if os.path.exists(corrected_assembly_file):
        bak_assembly_file = '{}.bak.{}'.format(corrected_assembly_file, time.time())
        logger.info('File {} already exists! Rename it as {}'.format(corrected_assembly_file, bak_assembly_file))
        os.rename(corrected_assembly_file, bak_assembly_file)

    if nbroken_ctgs:
        logger.info('{} contigs were broken into {} contigs. Writing corrected assembly to {}...'.format(
            nbroken_ctgs, len(fa_dict) - len(unbroken_ctgs), corrected_assembly_file))
        # output corrected assembly in FASTA format
        with open(corrected_assembly_file, 'w') as f:
            for ctg, ctg_info in fa_dict.items():
                f.write('>{}\n{}\n'.format(ctg, ctg_info[0]))
        # output a list for corrected contigs
        with open(corrected_ctgs_file, 'w') as f:
            for ctg, ctg_info in fa_dict.items():
                if ctg not in unbroken_ctgs:
                    assert ':' in ctg
                    f.write(ctg + '\n')
        # generate new gfa files for reassignment (for quick view only)
        if args.quick_view and read_depth_dict and len(args.gfa.split(',')) >= 2:
            for n, gfa in enumerate(args.gfa.split(',')):
                corrected_gfa = 'corrected_' + os.path.basename(gfa)
                with open(corrected_gfa, 'w') as f:
                    for ctg, (hap, read_depth) in read_depth_dict.items():
                        if hap == n:
                            # only keep necessary information
                            f.write('S\t{}\t*\tLN:i:{}\trd:i:{}\n'.format(ctg, fa_dict[ctg][1], read_depth))
    else:
        logger.info('No corrected contigs were found. Simply create a symbolic link of the input assembly')
        os.symlink(args.fasta, corrected_assembly_file)
        # generate an empty corrected_ctgs.txt
        with open(corrected_ctgs_file, 'w') as f:
            pass
        # generate new gfa files for reassignment (for quick view only)
        if args.quick_view and read_depth_dict and len(args.gfa.split(',')) >= 2:
            for gfa in args.gfa.split(','):
                corrected_gfa = 'corrected_' + os.path.basename(gfa)
                os.symlink(gfa, corrected_gfa)

    return nbroken_ctgs, final_break_pos_dict, final_break_frag_dict


def parse_pairs_for_correction(fa_dict, args):

    """parse pairs file for contig correction"""
    logger.info('Parsing input pairs file for contig correction...')

    resolution = args.correct_resolution

    ctg_cov_dict = dict()
    ctg_link_pos_dict = defaultdict(lambda: array('i'))

    for ctg, ctg_info in fa_dict.items():
        ctg_cov_dict[ctg] = ndarray([0] * (ctg_info[1]//resolution+1), dtype=int32)

    if args.aln_format == 'pairs':
        fopen = open
    else:
        assert args.aln_format == 'bgzipped_pairs'
        fopen = gzip.open

    with fopen(args.alignments, 'rt') as f:

        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            ref, mref = cols[1], cols[3]

            if ref != mref:
                continue

            # skip contigs that not in the fasta file
            if ref not in fa_dict:
                continue

            pos, mpos = int(cols[2]) - 1, int(cols[4]) - 1

            # get the spanning region of an alignment
            spanning_region = sorted([pos, mpos])
            spanning_start_bin = spanning_region[0]//resolution
            spanning_end_bin = spanning_region[1]//resolution

            ctg_cov_dict[ref][spanning_start_bin:spanning_end_bin+1] += 1
            ctg_link_pos_dict[ref].extend(spanning_region)

    return ctg_cov_dict, ctg_link_pos_dict


def check_sorting_order(fbam):

    hd_line = fbam.header.get('HD')

    if hd_line is not None and 'SO' in hd_line:
        if hd_line['SO'] in {'unsorted', 'queryname'}:
            logger.info('The sorting order of the BAM file is {}'.format(hd_line['SO']))
            return
        elif hd_line['SO'] == 'coordinate':
            logger.error('The sorting order of the BAM file is {}. It should be unsorted or name-sorted'.format(hd_line['SO']))
            raise RuntimeError('The sorting order of the BAM file is {}. It should be unsorted or name-sorted'.format(hd_line['SO']))

    logger.warning('The sorting order of the BAM file is unknown, but the program will continue')


def parse_bam_for_correction(fa_dict, args):

    """parse BAM file for contig correction"""
    # It will make the program read the bam file twice, but it's easier for implementation
    logger.info('Parsing input BAM file for contig correction...')

    resolution = args.correct_resolution

    ctg_cov_dict = dict()
    ctg_link_pos_dict = defaultdict(lambda: array('i'))

    for ctg, ctg_info in fa_dict.items():
        ctg_cov_dict[ctg] = ndarray([0] * (ctg_info[1]//resolution+1), dtype=int32)

    format_options = [b'filter=flag.read1 && refid == mrefid']

    with pysam.AlignmentFile(args.alignments, mode='rb', threads=args.threads, format_options=format_options) as f:
        check_sorting_order(f)
        for aln in f:

            ref = aln.reference_name

            # skip contigs that not in the fasta file
            if ref not in fa_dict:
                continue

            pos, mpos = aln.reference_start, aln.next_reference_start

            # get the spanning region of an alignment
            spanning_region = sorted([pos, mpos])
            spanning_start_bin = spanning_region[0]//resolution
            spanning_end_bin = spanning_region[1]//resolution

            ctg_cov_dict[ref][spanning_start_bin:spanning_end_bin+1] += 1
            ctg_link_pos_dict[ref].extend(spanning_region)

    return ctg_cov_dict, ctg_link_pos_dict


def pairs_generator_for_correction_ctg(pairs, aln_format, final_break_pos_dict, final_break_frag_dict):

    # convert the coordinates of Hi-C links in the contigs that have been corrected (when no contigs are split into bins)

    def convert_ctg(ctg, coord):
        n = 0
        for pos in final_break_pos_dict[ctg]:
            new_coord = coord - pos
            if new_coord >= 0:
                return final_break_frag_dict[ctg][n], new_coord
            n += 1

    if aln_format == 'pairs':
        fopen = open
    else:
        assert aln_format == 'bgzipped_pairs'
        fopen = gzip.open

    # generate a BED file for juice pre
    # The file is not completely accurate as some information is missing in the pairs files compared to the bam files
    with fopen(pairs, 'rt') as f, open('alignments.bed', 'w') as fbed:

        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            # pysam and BED format use the 0-based coordinate system
            ref, pos, mref, mpos = cols[1], int(cols[2]) - 1, cols[3], int(cols[4]) -1
            fbed.write('{0}\t{1}\t{2}\t{3}/1\t255\t.\n{4}\t{5}\t{6}\t{3}/2\t255\t.\n'.format(ref, pos, pos, cols[0], mref, mpos, mpos))

            if ref in final_break_frag_dict:
                ref, pos = convert_ctg(ref, pos)

            if mref in final_break_frag_dict:
                mref, mpos = convert_ctg(mref, mpos)

            if ref == mref:
                continue

            yield ref, mref, pos, mpos


def bam_generator_for_correction_ctg(bam, threads, format_options, final_break_pos_dict, final_break_frag_dict):

    # convert the coordinates of Hi-C links in the contigs that have been corrected (when no contigs are split into bins)

    def convert_ctg(ctg, coord):
        n = 0
        for pos in final_break_pos_dict[ctg]:
            new_coord = coord - pos
            if new_coord >= 0:
                return final_break_frag_dict[ctg][n], new_coord
            n += 1

    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as f:
        check_sorting_order(f)
        for aln in f:
            ref, mref, pos, mpos = (
                    aln.reference_name, aln.next_reference_name, aln.reference_start, aln.next_reference_start)

            if ref in final_break_frag_dict:
                ref, pos = convert_ctg(ref, pos)

            if mref in final_break_frag_dict:
                mref, mpos = convert_ctg(mref, mpos)

            if ref == mref:
                continue

            yield ref, mref, pos, mpos


def pairs_generator_for_correction(pairs, aln_format, final_break_pos_dict, final_break_frag_dict):

    # convert the coordinates of Hi-C links in the contigs that have been corrected

    def convert_ctg(ctg, coord):
        n = 0
        for pos in final_break_pos_dict[ctg]:
            new_coord = coord - pos
            if new_coord >= 0:
                return final_break_frag_dict[ctg][n], new_coord
            n += 1

    if aln_format == 'pairs':
        fopen = open
    else:
        assert aln_format == 'bgzipped_pairs'
        fopen = gzip.open

    # generate a BED file for juice pre
    # The file is not completely accurate as some information is missing in the pairs files compared to the bam files
    with fopen(pairs, 'rt') as f, open('alignments.bed', 'w') as fbed:

        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            # pysam and BED format use the 0-based coordinate system
            ref, pos, mref, mpos = cols[1], int(cols[2]) - 1, cols[3], int(cols[4]) - 1
            fbed.write('{0}\t{1}\t{2}\t{3}/1\t255\t.\n{4}\t{5}\t{6}\t{3}/2\t255\t.\n'.format(ref, pos, pos, cols[0], mref, mpos, mpos))

            if ref in final_break_frag_dict:
                ref, pos = convert_ctg(ref, pos)

            if mref in final_break_frag_dict:
                mref, mpos = convert_ctg(mref, mpos)

            yield ref, mref, pos, mpos


def bam_generator_for_correction(bam, threads, format_options, final_break_pos_dict, final_break_frag_dict):

    # convert the coordinates of Hi-C links in the contigs that have been corrected

    def convert_ctg(ctg, coord):
        n = 0
        for pos in final_break_pos_dict[ctg]:
            new_coord = coord - pos
            if new_coord >= 0:
                return final_break_frag_dict[ctg][n], new_coord
            n += 1

    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as f:
        check_sorting_order(f)
        for aln in f:
            ref, mref, pos, mpos = (
                    aln.reference_name, aln.next_reference_name, aln.reference_start, aln.next_reference_start)

            if ref in final_break_frag_dict:
                ref, pos = convert_ctg(ref, pos)

            if mref in final_break_frag_dict:
                mref, mpos = convert_ctg(mref, mpos)

            yield ref, mref, pos, mpos


def pairs_generator(pairs, aln_format):

    if aln_format == 'pairs':
        fopen = open
    else:
        assert aln_format == 'bgzipped_pairs'
        fopen = gzip.open

    # generate a BED file for juice pre
    # The file is not completely accurate as some information is missing in the pairs files compared to the bam files
    with fopen(pairs, 'rt') as f, open('alignments.bed', 'w') as fbed:

        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            # pysam and BED format use the 0-based coordinate system
            ref, pos, mref, mpos = cols[1], int(cols[2]) - 1, cols[3], int(cols[4]) - 1
            fbed.write('{0}\t{1}\t{2}\t{3}/1\t255\t.\n{4}\t{5}\t{6}\t{3}/2\t255\t.\n'.format(ref, pos, pos, cols[0], mref, mpos, mpos))

            yield ref, mref, pos, mpos


def pairs_generator_inter_ctgs(pairs, aln_format):

    if aln_format == 'pairs':
        fopen = open
    else:
        assert aln_format == 'bgzipped_pairs'
        fopen = gzip.open

    # generate a BED file for juice pre
    # The file is not completely accurate as some information is missing in the pairs files compared to the bam files
    with fopen(pairs, 'rt') as f, open('alignments.bed', 'w') as fbed:

        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            # pysam and BED format use the 0-based coordinate system
            ref, pos, mref, mpos = cols[1], int(cols[2]) - 1, cols[3], int(cols[4]) - 1
            fbed.write('{0}\t{1}\t{2}\t{3}/1\t255\t.\n{4}\t{5}\t{6}\t{3}/2\t255\t.\n'.format(ref, pos, pos, cols[0], mref, mpos, mpos))

            if ref != mref:
                yield ref, mref, pos, mpos


def bam_generator(bam, threads, format_options):

    # just a wrapper of pysam.AlignmentFile used to keep the style of Hi-C link parsing consistent

    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as f:
        check_sorting_order(f)
        for aln in f:
            yield aln.reference_name, aln.next_reference_name, aln.reference_start, aln.next_reference_start


def parse_alignments_for_ctgs(alignments, fa_dict, args, ctg_len_dict, Nx_ctg_set):

    """parse alignments (when no contigs are split into bins)"""
    logger.info('Parsing input alignments...')

    # set default parameters
    # kbp -> bp
    flank = args.flank * 1000

    full_link_dict = defaultdict(int)

    # a dict storing coords of Hi-C links between contigs
    ctg_coord_dict = defaultdict(lambda: array('i'))

    # if not reassignment:
    flank_link_dict = defaultdict(int)
    HT_link_dict = defaultdict(int)
    ctg_link_dict = defaultdict(int)
    clm_dict = defaultdict(lambda: array('i'))

    for ref, mref, pos, mpos in alignments:

        # skip contigs that not in the fasta file
        if ref not in fa_dict or mref not in fa_dict:
            continue

        # sort by contig name
        (ctg_i, coord_i), (ctg_j, coord_j) = sorted(((ref, pos+1), (mref, mpos+1)))
        ctg_name_pair = (ctg_i, ctg_j)

        # get contig index pairs
        len_i, len_j = ctg_len_dict[ctg_i], ctg_len_dict[ctg_j]

        # update flank_link_dict & ctg_link_dict
        if ctg_i in Nx_ctg_set and ctg_j in Nx_ctg_set and is_flank(coord_i, len_i, flank) and is_flank(coord_j, len_j, flank):
            flank_link_dict[ctg_name_pair] += 1
            ctg_link_dict[ctg_i] += 1
            ctg_link_dict[ctg_j] += 1

        # statistics for clm file and HT link pickle file output
        # clm file
        update_clm_dict(clm_dict, ctg_name_pair, len_i, len_j, coord_i-1, coord_j-1)

        # HT link pickle file
        update_HT_link_dict(HT_link_dict, ctg_i, ctg_j, len_i, len_j, coord_i, coord_j)

        # update full_link_dict
        full_link_dict[ctg_name_pair] += 1

        # record coord pairs and calculate concordance ratios / concentration adjustment ratio once there are enough coord pairs
        if args.remove_allelic_links or args.remove_concentrated_links:
            record_coord_pairs(ctg_coord_dict, ctg_name_pair, coord_i, coord_j, args.max_read_pairs, fa_dict, args)

    return full_link_dict, flank_link_dict, HT_link_dict, clm_dict, ctg_link_dict, ctg_coord_dict


def parse_alignments(alignments, fa_dict, args, bin_size, frag_len_dict, Nx_frag_set, split_ctg_set):

    """parse alignments (when some contigs are split into bins)"""

    def convert_frags(ctg, coord):

        if ctg in split_ctg_set:
            nbins = ceil(coord / bin_size)
            bin_name = '{}_bin{}'.format(ctg, nbins)
            bin_coord = coord-(nbins-1)*bin_size
            return bin_name, bin_coord, True
        else:
            return ctg, coord, False

    logger.info('Parsing input alignments...')

    # kbp -> bp
    flank = args.flank * 1000

    full_link_dict = defaultdict(int)

    # a dict storing coords of Hi-C links between contigs
    ctg_coord_dict = defaultdict(lambda: array('i'))
    # a dict mapping ctg pairs to frag pairs
    ctg_pair_to_frag = defaultdict(set)


    flank_link_dict = defaultdict(int)
    HT_link_dict = defaultdict(int)
    frag_link_dict = defaultdict(int)
    clm_dict = defaultdict(lambda: array('i'))

    for ref, mref, pos, mpos in alignments:

        # skip intra-bin Hi-C links
        if ref == mref and ref not in split_ctg_set:
            continue

        # skip contigs that not in the fasta file
        if ref not in fa_dict or mref not in fa_dict:
            continue

        # sort by contig name
        (ctg_i, coord_i), (ctg_j, coord_j) = sorted(((ref, pos+1), (mref, mpos+1)))
        ctg_name_pair = (ctg_i, ctg_j)

        # convert names and coords
        frag_i, frag_coord_i, i_is_bin = convert_frags(ctg_i, coord_i)
        frag_j, frag_coord_j, j_is_bin = convert_frags(ctg_j, coord_j)

        # intra-bin Hi-C links are not considered
        if frag_i == frag_j:
            continue

        # sort by bin name
        if i_is_bin or j_is_bin:
            (frag_i, frag_coord_i), (frag_j, frag_coord_j) = sorted(((frag_i, frag_coord_i), (frag_j, frag_coord_j)))

        frag_name_pair = (frag_i, frag_j)
        frag_len_i, frag_len_j = frag_len_dict[frag_i], frag_len_dict[frag_j]

        # update flank_link_dict & frag_link_dict
        if frag_i in Nx_frag_set and frag_j in Nx_frag_set and is_flank(frag_coord_i, frag_len_i, flank) and is_flank(frag_coord_j, frag_len_j, flank):
            flank_link_dict[frag_name_pair] += 1
            frag_link_dict[frag_i] += 1
            frag_link_dict[frag_j] += 1

        # udpate ctg_pair_to_frag
        if args.remove_allelic_links:
            ctg_pair_to_frag[ctg_name_pair].add(frag_name_pair)

        # intra-contig Hi-C links are NOT considered
        if ref != mref:

            # update full_link_dict
            full_link_dict[ctg_name_pair] += 1

            # statistics for clm file and HT link pickle file output
            len_i, len_j = fa_dict[ctg_i][1], fa_dict[ctg_j][1]
            # clm file
            update_clm_dict(clm_dict, ctg_name_pair, len_i, len_j, coord_i-1, coord_j-1)
            # HT link pickle file
            update_HT_link_dict(HT_link_dict, ctg_i, ctg_j, len_i, len_j, coord_i, coord_j)

            # record coord pairs and calculate concordance ratios / concentration adjustment ratio once there are enough coord pairs
            if args.remove_allelic_links or args.remove_concentrated_links:
                record_coord_pairs(ctg_coord_dict, ctg_name_pair, coord_i, coord_j, args.max_read_pairs, fa_dict, args)

    return full_link_dict, flank_link_dict, HT_link_dict, clm_dict, frag_link_dict, ctg_coord_dict, ctg_pair_to_frag


def add_edge(graph, node1, node2):

    if not graph.has_edge(node1, node2):
        graph.add_edge(node1, node2, weight=1)
    else:
        graph[node1][node2]['weight'] += 1


def parse_ul_alignments(args):

    """parse ultra-long alignments"""

    logger.info('Parsing input ultra-long alignments...')

    G_HT = Graph()

    def get_query_alignment_termini(aln):

        if aln.is_forward:
            # '5' means hard-clipping
            if aln.cigartuples[0][0] == 5:
                hard_clip = aln.cigartuples[0][1]
            else:
                hard_clip = 0
            return aln.query_alignment_start + hard_clip, aln.query_alignment_end + hard_clip
        else:
            if aln.cigartuples[-1][0] == 5:
                hard_clip = aln.cigartuples[-1][1]
            else:
                hard_clip = 0
            read_len = aln.infer_read_length()
            return read_len - aln.query_alignment_end + hard_clip, read_len - aln.query_alignment_start + hard_clip

    def parse_supplementary_aln_list(primary_aln, supplementary_aln_list, f):

        # get the best supplementary alignment based on the AS tag (score)
        if len(supplementary_aln_list) > 1:
            supplementary_aln_list.sort(key=lambda x: x.get_tag('AS'), reverse=True)
        supplementary_aln = supplementary_aln_list[0]

        # link semi-contigs using UL reads
        semi_ctg_list = [[(primary_aln, '_H'), (primary_aln, '_T')], [(supplementary_aln, '_H'), (supplementary_aln, '_T')]]
        # get the adjacent semi-contigs of different contigs
        if primary_aln.is_reverse:
            semi_ctg_list[0].reverse()
        if supplementary_aln.is_reverse:
            semi_ctg_list[1].reverse()
        semi_ctg_list.sort(key=lambda x: get_query_alignment_termini(x[0][0])[0])
        # semi_ctg_list.sort(key=lambda x: get_query_alignment_termini(x[0][0])[0], reverse=primary_aln.is_reverse)
        left_semi_ctg = semi_ctg_list[0][1][0].reference_name + semi_ctg_list[0][1][1]
        right_semi_ctg = semi_ctg_list[1][0][0].reference_name + semi_ctg_list[1][0][1]
        # add edges for semi-contigs of different contigs
        add_edge(G_HT, left_semi_ctg, right_semi_ctg)
        # to avoid breaking inside a contig, the edge weight within a contig should always be equal to or greater than the edge weight between different contigs
        add_edge(G_HT, primary_aln.reference_name+'_H',  primary_aln.reference_name+'_T')
        add_edge(G_HT, supplementary_aln.reference_name+'_H',  supplementary_aln.reference_name+'_T')

    min_ul_mapq = args.min_ul_mapq
    min_ul_alignment_length = args.min_ul_alignment_length
    max_distance_to_end = args.max_distance_to_end
    max_overlap_ratio = args.max_overlap_ratio
    max_gap_len = args.max_gap_len
    min_ul_support = args.min_ul_support

    primary_set = {0, 16}
    primary_aln = None
    supplementary_aln_list = []

    with pysam.AlignmentFile(args.ul, 'rb', format_options=[b'filter=!flag.unmap'], threads=args.threads) as f:
        for aln in f:
            # MAPQ and alignment length filtering
            if aln.mapq < min_ul_mapq or aln.reference_length < min_ul_alignment_length:
                continue
            # the alignment should be close to at least one of the end of the reference
            if aln.reference_start > max_distance_to_end and f.get_reference_length(aln.reference_name) - aln.reference_end > max_distance_to_end:
                continue
            # calculate start and end coordinates on the read based on the alignment orientation
            query_start, query_end = get_query_alignment_termini(aln)

            # # if the alignment is close to both ends of the reference, the alignment should be far away from the ends of the query
            # if aln.reference_start <= max_distance_to_end and f.get_reference_length(aln.reference_name) - aln.reference_end <= max_distance_to_end:
            #     if query_start <= max_distance_to_end or read_length - query_end <= max_distance_to_end:
            #         continue

            # for primary alignments
            if aln.flag in primary_set:
                if supplementary_aln_list:
                    parse_supplementary_aln_list(primary_aln, supplementary_aln_list, f)
                primary_aln = aln
                primary_query_start, primary_query_end = query_start, query_end
                supplementary_aln_list.clear()
            # for supplementary alignments
            elif aln.is_supplementary and primary_aln and aln.query_name == primary_aln.query_name and aln.reference_name != primary_aln.reference_name:
                # read alignment interval filtering
                primary_read_interval = closed(primary_query_start + 1, primary_query_end)
                supplementary_read_interval = closed(query_start + 1, query_end)
                overlap = primary_read_interval & supplementary_read_interval
                # if there is a overlap in read intervals between primary and supplementary alignments
                if overlap:
                    overlap_len = overlap.upper - overlap.lower + 1
                    primary_interval_len = primary_read_interval.upper - primary_read_interval.lower + 1
                    supplementary_interval_len = supplementary_read_interval.upper - supplementary_read_interval.lower + 1
                    # the ratio of overlap length to the shorter one of primary and supplementary intervals should <= max_overlap_ratio
                    if overlap_len / min(primary_interval_len, supplementary_interval_len) > max_overlap_ratio:
                        continue
                # if there is no overlap, the gap between primary and supplementary intervals should <= max_gap_len
                else:
                    gap = (primary_read_interval | supplementary_read_interval).complement()[1]
                    if gap.upper - gap.lower + 1 > max_gap_len:
                        continue
                supplementary_aln_list.append(aln)

        # the last alignment
        if supplementary_aln_list:
            parse_supplementary_aln_list(primary_aln, supplementary_aln_list, f)

        ## get the contig linking path
        # remove edges with a weight less than min_ul_support
        for node1, node2, d in deepcopy(G_HT).edges(data=True):
            if d['weight'] < min_ul_support and node1.rsplit('_', 0) != node2.rsplit('_', 0):
                G_HT.remove_edge(node1, node2)

        # remove edges linked to a node with degree > 2
        G_HT_copy = deepcopy(G_HT)
        for node1, node2 in G_HT_copy.edges():
            if (G_HT_copy.degree(node1) > 2 or G_HT_copy.degree(node2) > 2) and node1.rsplit('_', 1)[0] != node2.rsplit('_', 1)[0]:
                G_HT.remove_edge(node1, node2)

        # find connected subgraphs
        path_list = []
        connected_subgraphs = connected_components(G_HT)
        for subnodes in connected_subgraphs:
            # link at least two contigs
            if len(subnodes) < 4:
                continue
            logger.debug([(node, G_HT.degree(node)) for node in subnodes])
            assert len(subnodes) % 2 == 0
            subgraph = G_HT.subgraph(subnodes).copy()
            degree_list = [subgraph.degree(node) for node in subnodes]
            # linear
            if degree_list.count(1) == 2:
                node1, node2 = [node for node in subnodes if subgraph.degree(node) == 1]
            # circular
            else:
                assert degree_list.count(1) == 0 and degree_list.count(2) == len(degree_list)
                # break the edge of minimum weight
                weight_list = [(node1, node2, G_HT[node1][node2]['weight']) for node1, node2 in subgraph.edges()]
                weight_list.sort(key=lambda x: x[2])
                node1, node2, _ = weight_list[0]
                subgraph.remove_edge(node1, node2)
            path = shortest_path(subgraph)[node1][node2]
            path_list.append(path)
            logger.debug('{}\t{}'.format('->'.join(path), '-'.join([str(f.get_reference_length(node.rsplit('_', 1)[0])//2) for node in path])))

    return path_list


def add_HT_links_based_on_ul(path_list, HT_link_dict):

    for path in path_list:
        for i in range(len(path) - 1):
            if i % 2 != 0:
                # get linked semi-contigs between different contigs
                node1, node2 = path[i], path[i+1]

                # get corresponding contig names
                ctg1 = node1.rsplit('_', 1)[0]
                ctg2 = node2.rsplit('_', 1)[0]

                # sort by contig name
                (ctg1, node1), (ctg2, node2) = sorted(((ctg1, node1), (ctg2, node2))) 

                # update HT_link_dict
                assert (node2, node1) not in HT_link_dict
                if (node1, node2) in HT_link_dict:
                    logger.debug('update HT_link_dict: {} {}'.format(node1, node2))
                    HT_link_dict[(node1, node2)] *= 2
                else:
                    logger.debug('{} {} not in HT_link_dict'.format(node1, node2))


def add_flank_and_full_links_based_on_ul(path_list, flank_link_dict, full_link_dict, bin_set):

    ul_linked_ctgs = []
    for path in path_list:
        ul_linked_ctgs.append(set())
        for i in range(len(path) - 1):
            if i % 2 != 0:
                # get linked semi-contigs between different contigs
                node1, node2 = path[i], path[i+1]

                # get corresponding contig names
                ctg1 = node1.rsplit('_', 1)[0]
                ctg2 = node2.rsplit('_', 1)[0]
                ul_linked_ctgs[-1].add(ctg1)
                ul_linked_ctgs[-1].add(ctg2)

                # sort by contig name
                (ctg1, node1), (ctg2, node2) = sorted(((ctg1, node1), (ctg2, node2))) 

                # update full_link_dict
                assert (ctg2, ctg1) not in full_link_dict
                if (ctg1, ctg2) in full_link_dict:
                    logger.debug('update full_link_dict: {} {}'.format(ctg1, ctg2))
                    full_link_dict[(ctg1, ctg2)] *= 2
                else:
                    logger.debug('{} {} not in full_link_dict'.format(node1, node2))

    ul_linked_ctg_pairs = set()
    for ctgs in ul_linked_ctgs:
        for ctg1, ctg2 in combinations(ctgs, 2):
            ul_linked_ctg_pairs.add((ctg1, ctg2))
            ul_linked_ctg_pairs.add((ctg2, ctg1))

    for frag_i, frag_j in flank_link_dict:
        if frag_i in bin_set:
            ctg_i = frag_i.rsplit('_bin', 1)[0]
        else:
            ctg_i = frag_i
        if frag_j in bin_set:    
            ctg_j = frag_j.rsplit('_bin', 1)[0]
        else:
            ctg_j = frag_j
        # update flank_link_dict
        if (ctg_i, ctg_j) in ul_linked_ctg_pairs:
            if (frag_i, frag_j) in flank_link_dict:
                print('update flank_link_dict: {} {}'.format(frag_i, frag_j))
                flank_link_dict[(frag_i, frag_j)] *= 2
            else:
                print('{} {} not in flank_link_dict'.format(frag_i, frag_j))


def prune(matrix, pruning, dense_matrix):

    # when the matrix is sparse enough (otherwise the value assignment of dok will be very slow)
    if not dense_matrix and matrix.nnz/matrix.shape[0]**2 < 0.05:
        # create an empty sparse matrix
        pruned_matrix = dok_matrix(matrix.shape, dtype=float32)
        # set value of each element to the corresponding value of raw matrix if larger than pruning threshold
        bool_matrix = matrix >= pruning
        pruned_matrix[bool_matrix] = matrix[bool_matrix]
        # re-convert dok matrix to csc
        pruned_matrix = pruned_matrix.tocsc()
    else:
        # make a copy for data manipulation
        if not dense_matrix:
            pruned_matrix = matrix.toarray()
        else:
            pruned_matrix = matrix.copy()
        # set value of each element to zero if is smaller than pruning threshold
        pruned_matrix[pruned_matrix < pruning] = 0
        if not dense_matrix:
            pruned_matrix = csc_matrix(pruned_matrix)

    # the maximum value of each column should be kept to ensure the sum of column > 0
    ncols = matrix.shape[1]
    cols = arange(ncols)
    rows = matrix.argmax(axis=0).reshape((ncols,))
    pruned_matrix[rows, cols] = matrix[rows, cols]
    return normalize(pruned_matrix, norm='l1', axis=0)


def mkl_matrix_power(matrix, n):

    # a recursive function achieving matrix power in mkl
    if n == 2:
        return dot_product_mkl(matrix, matrix)
    else:
        return dot_product_mkl(matrix, mkl_matrix_power(matrix, n-1))


def mcl(matrix, expansion, inflation, iters, pruning, dense_matrix):

    # iterate until convergence
    for n in range(iters):
        if n != 0:
            # 2) expand
            if not dense_matrix:
                matrix = mkl_matrix_power(matrix, expansion)
            else:
                matrix = matrix_power(matrix, expansion)
        # 3) inflate
        if not dense_matrix:
            matrix = normalize(matrix.power(inflation), norm='l1', axis=0)
        else:
            matrix = normalize(power(matrix, inflation), norm='l1', axis=0)
        # 4) prune
        matrix = prune(matrix, pruning, dense_matrix)
        # 5) convergence verification
        if n > 1 and not dense_matrix:
            d = npabs(matrix-last_matrix)-1e-5*npabs(last_matrix)
            if d.max() <= 1e-8:
                logger.info('The matrix has converged after {} rounds of iterations ' \
                        '(expansion: {}, inflation: {}, maximum iterations: {}, pruning threshold: {})'.format(
                            n+1, expansion, inflation, iters, pruning))
                return matrix
        elif n > 1 and allclose(matrix, last_matrix):
            logger.info('The matrix has converged after {} rounds of iterations ' \
                    '(expansion: {}, inflation: {}, maximum iterations: {}, pruning threshold: {})'.format(
                        n+1, expansion, inflation, iters, pruning))
            return matrix
        # make a copy for convergence verification in the next iteration
        last_matrix = matrix.copy()
    logger.info('The matrix does not converge after {} rounds of iterations ' \
            '(expansion: {}, inflation: {}, maximum iterations: {}, pruning threshold: {})'.format(
                n+1, expansion, inflation, iters, pruning))

    return matrix


def interpret_result(result_matrix, dense_matrix):

    # get the shape of the matrix
    shape = result_matrix.shape[0]

    # attractors are rows in which the elements are diagonal distributed and non-zero
    attractors = result_matrix.diagonal().nonzero()[0]

    clusters = set()

    for attractor in attractors:
        # get the col number of the non-zero element in the row of attractor
        if not dense_matrix:
            cluster = tuple(result_matrix.getrow(attractor).nonzero()[1].tolist())
        else:
            cluster = tuple([n for n in range(shape) if result_matrix[attractor,n] != 0])

        clusters.add(cluster)

    # make sure:
    # (1) all nodes are in clusters;
    # (2) no repeated nodes
    nodes = set()
    for cluster in clusters:
        for n in cluster:
            if n in nodes:
                return None
            nodes.add(n)
    if len(nodes) != shape:
        return None
    return list(clusters)


def get_main_groups(result_clusters, len_ratio):

    len_sum = 0
    main_groups = len(result_clusters)
    for n in range(len(result_clusters)-1):
        len_sum += result_clusters[n][1]
        if result_clusters[n+1][1]/result_clusters[n][1] < len_ratio:
            main_groups = n + 1
            break
    return main_groups


def recommend_inflation(result_stat, nchrs, len_ratio):

    separated = [(inflation, main_groups-nchrs) for inflation, main_groups in result_stat if main_groups >= nchrs]

    if separated:
        # separated.sort(key=lambda x: (x[1], x[0]))
        separated.sort(key=lambda x: x[0])
        rcm_inflation = separated[0][0]
        logger.info('You could try inflation from {} (length ratio = {})'.format(rcm_inflation, len_ratio))
        return True
    else:
        if len_ratio > 0.5:
            logger.info('The length ratio ({}) might be too strict, trying a lower one...'.format(len_ratio))
            return False
        else:
            logger.info('It seems that some chromosomes were grouped together (length ratio = {}) '
                        'You could check whether the parameters used are correct / appropriate and '
                        'then try to tune the parameters for assembly correction, contig / Hi-C link '
                        'filtration, or Markov clustering'.format(len_ratio))
            return True


def run_mcl_clustering(link_matrix, bin_set, frag_len_dict, frag_index_dict, expansion, min_inflation,
        max_inflation, inflation_step, max_iter, pruning, fa_dict, nchrs, dense_matrix):

    logger.info('Performing Markov clustering...')

    index_frag_dict = {i: f for f, i in frag_index_dict.items()}

    start = Decimal(str(min_inflation))
    step = Decimal(str(inflation_step))
    end = Decimal(str(max_inflation)) + step

    # 1) pre-normalize the count matrix to probability matrix
    matrix = normalize(link_matrix, norm='l1', axis=0)
    # pre-expand the probability matrix
    if not dense_matrix:
        matrix = mkl_matrix_power(matrix, expansion)
    else:
        matrix = matrix_power(matrix, expansion)
    result_clusters_list = list()

    mcl_nrounds = 0

    # try different inflation values
    for inflation in arange(start, end, step):

        # run Markov clustering
        result_matrix = mcl(matrix, expansion, float(inflation), max_iter, pruning, dense_matrix)
        mcl_nrounds += 1

        # get clusters
        clusters = interpret_result(result_matrix, dense_matrix)
        if not clusters:
            logger.info('Some fragments are missing / redundant, result of inflation {} will NOT be output'.format(inflation))
            continue

        # a dict storing contigs and length for each cluster
        result_clusters = defaultdict(lambda: [[], 0])
        # a dict to determine the best cluster for each contig split into bins
        ctg_clusters = defaultdict(dict)

        for n, indexes in enumerate(clusters):
            for i in indexes:
                frag = index_frag_dict[i]
                # for bins
                if frag in bin_set:
                    frag_len = frag_len_dict[frag]
                    # TODO: there may be a risk, if input fasta IDs contain '_bin'?
                    ctg = frag.rsplit('_bin', 1)[0]
                    if n in ctg_clusters[ctg]:
                        ctg_clusters[ctg][n] += frag_len
                    else:
                        ctg_clusters[ctg][n] = frag_len
                # for contigs
                else:
                    result_clusters[n][0].append(frag)
                    result_clusters[n][1] += fa_dict[frag][1]

        # assign contigs to best clusters & stat cluster lengths
        if len(ctg_clusters):
            for ctg, stats in ctg_clusters.items():
                best_cluster = sorted(stats.keys(), key=lambda x: stats[x], reverse=True)[0]
                result_clusters[best_cluster][0].append(ctg)
                result_clusters[best_cluster][1] += fa_dict[ctg][1]

        # sort all clusters by length
        result_clusters = sorted(tuple(result_clusters.values()), key=lambda x: x[1], reverse=True)

        # make directory for results
        outdir = 'inflation_{}'.format(inflation)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        # output a cluster file for result overview and reassignment
        with open('{0}/mcl_{0}.clusters.txt'.format(outdir), 'w') as fout:
            fout.write('#Group\tnContigs\tContigs\n')
            for n, (ctgs, group_len) in enumerate(result_clusters, 1):
                # sort contigs by length in each cluster
                result_clusters[n-1][0].sort(key=lambda x: fa_dict[x][1], reverse=True)
                fout.write('group{}_{}bp\t{}\t{}\n'.format(n, group_len, len(ctgs), ' '.join(ctgs)))

        # generate group*.txt files for ALLHiC optimize
        for n, (ctgs, group_len) in enumerate(result_clusters, 1):
            with open('{}/group{}_{}bp.txt'.format(outdir, n, group_len), 'w') as fout:
                fout.write('#Contig\tRECounts\tLength\n')
                for ctg in ctgs:
                    length, RE_sites = fa_dict[ctg][1:3]
                    fout.write('{}\t{}\t{}\n'.format(ctg, RE_sites, length))

        result_clusters_list.append((inflation, result_clusters))

    # get the maximum number of clusters
    max_nclusters = max([len(result_clusters) for inflation, result_clusters in result_clusters_list])

    if max_nclusters < nchrs:
        logger.warning('The maximum number of clusters ({}) is even less than the expected number of '
                    'chromosomes ({}). You could try higher inflation.'.format(max_nclusters, nchrs))

    else:
        for len_ratio in (0.75, 0.7, 0.65, 0.6, 0.55, 0.5):
            result_stat = list()
            # result statistics for recommendation
            for inflation, result_clusters in result_clusters_list:
                main_groups = get_main_groups(result_clusters, len_ratio)
                result_stat.append((inflation, main_groups))

            # get recommended inflation
            end_loop = recommend_inflation(result_stat, nchrs, len_ratio)
            if end_loop:
                break

    return result_clusters_list, mcl_nrounds


def add_ungrouped_ctgs(fa_dict, ctg_group_dict):

    for ctg in fa_dict:
        if ctg not in ctg_group_dict:
            ctg_group_dict[ctg] = 'ungrouped'


def parse_link_dict(link_dict, ctg_group_dict):

    def add_ctg_group(ctg, group):
        if group != 'ungrouped':
            if group in ctg_group_link_dict[ctg]:
                ctg_group_link_dict[ctg][group] += links
            else:
                ctg_group_link_dict[ctg][group] = links

    ctg_group_link_dict = defaultdict(dict)

    for (ctg_i, ctg_j), links in link_dict.items():
        group_i, group_j = ctg_group_dict[ctg_i], ctg_group_dict[ctg_j]
        add_ctg_group(ctg_i, group_j)
        add_ctg_group(ctg_j, group_i)

    return ctg_group_link_dict


def cal_link_density(max_group, current_group, max_links, group_RE_sites, ctg_RE_sites):

    if max_group == current_group:
        return max_links / group_RE_sites
    else:
        return max_links / (group_RE_sites + ctg_RE_sites - 1)


def output_statistics(fa_dict, link_dict, result_clusters_list):

    def generate_axes(sorted_list):
        filtered_ctg_n_dict = OrderedDict({0: 0})
        filtered_ctg_len_dict = OrderedDict({0: 0})

        last_value = 0
        for ctg, value in sorted_list:
            if value in filtered_ctg_n_dict:
                filtered_ctg_n_dict[value] += 1
                filtered_ctg_len_dict[value] += fa_dict[ctg][1]
            else:
                filtered_ctg_n_dict[value] = filtered_ctg_n_dict[last_value] + 1
                filtered_ctg_len_dict[value] = filtered_ctg_len_dict[last_value] + fa_dict[ctg][1]
                last_value = value

        x, y1, y2 = list(), list(), list()
        for k, v in filtered_ctg_n_dict.items():
            x.append(k)
            y1.append(v/total_ctg_n*100)
            y2.append((total_ctg_len-filtered_ctg_len_dict[k])/total_ctg_len*100)

        return x, y1, y2

    def write_result(x, y1, y2, title, inflation):
        with open('inflation_{}/{}_statistics.txt'.format(inflation, title), 'w') as fout:
            fout.write('{}\tFiltered_ctg_n\tRest_ctg_len\n'.format(title))
            for n, value in enumerate(x):
                fout.write('>{}\t{}\t{}\n'.format(value, y1[n], y2[n]))

    logger.info('Making some statistics for the next HapHiC reassignment step...')

    # contig RE site statistics
    total_ctg_n = len(fa_dict)
    total_ctg_len = 0
    # extract contig RE site info from fa_dict
    RE_site_list = list()
    for ctg, ctg_info in fa_dict.items():
        total_ctg_len += ctg_info[1]
        RE_site_list.append((ctg, ctg_info[2]))
    # sort by RE site
    RE_site_list.sort(key=lambda x: x[1])

    x_RE, y1_RE, y2_RE = generate_axes(RE_site_list)

    # matplotlib lazy import, although it's not PEP8 compliant
    try:
        import matplotlib
        # use non-GUI backend, save figure files only
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        # suppress annoying matplotlib warnings
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)
        MATPLOTLIB = True
    except:
        MATPLOTLIB = False
        logger.warning('Module matplotlib is not correctly installed, HapHiC will NOT draw statistical plots')

    for inflation, result_clusters in result_clusters_list:

        # print RE site statistics result
        write_result(x_RE, y1_RE, y2_RE, 'RE_site_threshold', inflation)

        # links & link dict statistics
        ctg_group_dict = dict()
        group_RE_dict = dict()

        for n, (ctgs, _) in enumerate(result_clusters):
            group_RE_dict[n] = 1
            for ctg in ctgs:
                ctg_group_dict[ctg] = n
                group_RE_dict[n] += fa_dict[ctg][2] - 1

        add_ungrouped_ctgs(fa_dict, ctg_group_dict)

        ctg_group_link_dict = parse_link_dict(link_dict, ctg_group_dict)

        ctg_max_group_link_list = list()
        ctg_max_group_link_density_list = list()
        ctg_max_avg_group_link_density_ratio_list = list()

        for ctg in fa_dict:
            if ctg in ctg_group_link_dict:
                group_links = ctg_group_link_dict[ctg]
                # links
                sorted_group_links = sorted(group_links.items(), key=lambda x: x[1], reverse=True)
                max_group, max_links = sorted_group_links[0]
                ctg_max_group_link_list.append((ctg, max_links))
                # link density
                current_group = ctg_group_dict[ctg]
                group_RE_sites = group_RE_dict[max_group]
                ctg_RE_sites = fa_dict[ctg][2]
                max_link_density = cal_link_density(max_group, current_group, max_links, group_RE_sites, ctg_RE_sites)
                ctg_max_group_link_density_list.append((ctg, max_link_density))
                # link density ratio
                if len(group_RE_dict) > 1:
                    other_group_density_sum = sum([cal_link_density(group, current_group, links, group_RE_dict[group], ctg_RE_sites) for group, links in sorted_group_links[1:]])
                    avg_other_group_density = other_group_density_sum / (len(group_RE_dict) - 1)
                else:
                    avg_other_group_density = 0
                if avg_other_group_density:
                    link_density_ratio = max_link_density/avg_other_group_density
                else:
                    link_density_ratio = 1000000
                ctg_max_avg_group_link_density_ratio_list.append((ctg, link_density_ratio))
            else:
                # links
                ctg_max_group_link_list.append((ctg, 0))
                # link density
                ctg_max_group_link_density_list.append((ctg, 0))
                # link density ratio
                ctg_max_avg_group_link_density_ratio_list.append((ctg, 0))

        # link
        ctg_max_group_link_list.sort(key=lambda x: x[1])
        x_link, y1_link, y2_link = generate_axes(ctg_max_group_link_list)
        write_result(x_link, y1_link, y2_link, 'Link_threshold', inflation)

        # link density
        ctg_max_group_link_density_list.sort(key=lambda x: x[1])
        x_density, y1_density, y2_density = generate_axes(ctg_max_group_link_density_list)
        write_result(x_density, y1_density, y2_density, 'Link_density_threshold', inflation)

        # link density ratio
        ctg_max_avg_group_link_density_ratio_list.sort(key=lambda x: x[1])
        x_ratio, y1_ratio, y2_ratio = generate_axes(ctg_max_avg_group_link_density_ratio_list)
        write_result(x_ratio, y1_ratio, y2_ratio, 'Link_density_ratio_threshold', inflation)

        # if matplotlib has been corrected installed
        if MATPLOTLIB:

            # generate statistics figure
            fig = plt.figure(figsize=(8, 7))

            # figure 1
            ax1_RE = fig.add_subplot(221)
            ax1_RE.plot(x_RE, y1_RE, 'b')
            ax1_RE.tick_params(axis='y', colors='b')
            ax1_RE.set_xlim([0, 500])
            ax1_RE.set_ylim([0, 50])
            ax1_RE.set_ylabel('Number of contigs filtered out (%)', color='b')
            ax1_RE.set_title('RE site threshold')
            ax1_RE.set_xlabel('Number of RE sites')
            ax2_RE = ax1_RE.twinx()
            ax2_RE.plot(x_RE, y2_RE, 'r')
            ax2_RE.tick_params(axis='y', colors='r')
            ax2_RE.set_ylim([90, 100])
            ax2_RE.set_ylabel('Length of remaining contigs (%)', color='r')

            # figure 2
            ax1_link = fig.add_subplot(222)
            ax1_link.plot(x_link, y1_link, 'b')
            ax1_link.tick_params(axis='y', colors='b')
            ax1_link.set_xlim([0, 500])
            ax1_link.set_ylim([0, 50])
            ax1_link.set_ylabel('Number of contigs filtered out (%)', color='b')
            ax1_link.set_title('Hi-C link threshold')
            ax1_link.set_xlabel('Number of links to the best group')
            ax2_link = ax1_link.twinx()
            ax2_link.plot(x_link, y2_link, 'r')
            ax2_link.tick_params(axis='y', colors='r')
            ax2_link.set_ylim([90, 100])
            ax2_link.set_ylabel('Length of remaining contigs (%)', color='r')

            # figure 3
            ax1_density = fig.add_subplot(223)
            ax1_density.plot(x_density, y1_density, 'b')
            ax1_density.tick_params(axis='y', colors='b')
            ax1_density.set_xlim([0, 0.001])
            ax1_density.set_ylim([0, 50])
            ax1_density.set_ylabel('Number of contigs filtered out (%)', color='b')
            ax1_density.set_title('Link density threshold')
            ax1_density.set_xlabel('Link density to the best group')
            ax2_density = ax1_density.twinx()
            ax2_density.plot(x_density, y2_density, 'r')
            ax2_density.tick_params(axis='y', colors='r')
            ax2_density.set_ylim([90, 100])
            ax2_density.set_ylabel('Length of remaining contigs (%)', color='r')

            # figure 4
            ax1_ratio = fig.add_subplot(224)
            ax1_ratio.plot(x_ratio, y1_ratio, 'b')
            ax1_ratio.tick_params(axis='y', colors='b')
            ax1_ratio.set_xlim([0, 20])
            ax1_ratio.set_ylim([0, 50])
            ax1_ratio.set_ylabel('Number of contigs filtered out (%)', color='b')
            ax1_ratio.set_title('Link density ratio threshold')
            ax1_ratio.set_xlabel('Link density ratio (best/average)')
            ax2_ratio = ax1_ratio.twinx()
            ax2_ratio.plot(x_ratio, y2_ratio, 'r')
            ax2_ratio.tick_params(axis='y', colors='r')
            ax2_ratio.set_ylim([90, 100])
            ax2_ratio.set_ylabel('Length of remaining contigs (%)', color='r')

            fig.tight_layout(w_pad=1, h_pad=1)

            plt.savefig('inflation_{}/statistics.pdf'.format(inflation))
            # plt.show()
            plt.close()


def check_param(param, string, suffix, true_suffix=''):

    if len(string) > 1:
        if suffix and string[-1] in suffix:
            return check_param(param, string[:-1], None, string[-1])
        else:
            try:
                num = float(string)
                if not true_suffix and not 0 <= num <= 1:
                    logger.error('Parameter {} {} is illegal'.format(param, string+true_suffix))
                else:
                    return num, true_suffix
            except ValueError:
                logger.error('Parameter {} {} is illegal'.format(param, string+true_suffix))
    elif len(string) == 1:
        try:
            num = float(string)
            if not true_suffix and not 0 <= num <= 1:
                logger.error('Parameter {} {} is illegal'.format(param, string+true_suffix))
            else:
                return num, true_suffix
        except ValueError:
            logger.error('Parameter {} {} is illegal'.format(param, string+true_suffix))
    else:
        logger.error('Parameter {} is empty'.format(param))

    raise RuntimeError('Parameter check failed')


def detect_format(args):

    # check file format for Hi-C read alignments

    if args.alignments.endswith('.bam'):
        args.aln_format = 'bam'
        logger.info('The file for Hi-C read alignments is detected as being in BAM format')

    elif args.alignments.endswith('.pairs'):
        args.aln_format = 'pairs'
        logger.info('The file for Hi-C read alignments is detected as being in pairs format')

    elif args.alignments.endswith('.pairs.gz'):
        args.aln_format = 'bgzipped_pairs'
        logger.info('The file for Hi-C read alignments is detected as being in bgzipped pairs format')

    else:
        raise RuntimeError('Unknown file format for Hi-C read alignments')


def parse_arguments():

    parser = argparse.ArgumentParser(prog='haphic cluster')

    # Parameters for parsing input files and pipeline control
    input_group = parser.add_argument_group('>>> Parameters for parsing input files and pipeline control')
    input_group.add_argument(
            'fasta', help='draft genome in FASTA format')
    input_group.add_argument(
            'alignments', help='filtered Hi-C read alignments in BAM/pairs format (DO NOT sort it by coordinate)')
    input_group.add_argument(
            'nchrs', type=int, help='expected number of chromosomes')
    input_group.add_argument(
            '--aln_format', choices={'bam', 'pairs', 'bgzipped_pairs', 'auto'}, default='auto', help='file format for Hi-C read alignments, default: %(default)s')
    input_group.add_argument(
            '--RE', default='GATC',
            help='restriction enzyme site(s) (e.g., GATC for MboI & AAGCTT for HindIII), default: %(default)s. If more than one enzyme was used '
            'in the Hi-C library construction, such as with the Arima genomics kit, separate the RE sites with commas (e.g., `--RE "GATC,GANTC"` for Arima '
            'two-enzyme chemistry and `--RE "GATC,GANTC,CTNAG,TTAA"` for Arima four-enzyme chemistry)')
    input_group.add_argument(
            '--quick_view', default=False, action='store_true',
            help='in quick view mode, HapHiC will skip the clustering and reassignment steps, and order and orient all contigs with fast sorting, default: %(default)s. '
            'This is helpful when you encounter difficulties in generating ideal clusters or when you are unsure of the exact number of chromosomes')
    input_group.add_argument(
            '--gfa', default=None,
            help='(experimental) matched GFA file(s) of the phased hifiasm assembly, separated with commas (e.g., `--gfa "gfa1,gfa2"`), default: %(default)s. '
            'Here, the "phased hifiasm assemblies" include the haplotype-resolved primary contigs assembled using the trio binning or Hi-C-based algorithm '
            '(`*.hap*.p_ctg.gfa`) and the phased unitigs (`*.p_utg.gfa`). HapHiC uses the read depth information in each gfa file to filter out contigs '
            '(see `--read_depth_upper`). If multiple GFA file is provided, HapHiC assumes these GFA files are haplotype-specific, and then artificially '
            'reduces the Hi-C links between the haplotypes according to this phasing information (see `--phasing_weight`). This parameter can also work with '
            '`--quick_view` to separate contigs into different haplotype groups')
    input_group.add_argument(
            '--ul', default=None,
            help='ultra-long read alignments in BAM format, default: %(default)s')

    # Parameters for assembly correction
    correct_group = parser.add_argument_group('>>> Parameters for assembly correction')
    correct_group.add_argument(
            '--correct_nrounds', type=int, default=0,
            help='maximum number of rounds for assembly correction, default: %(default)s. If no breakpoints are detected, the iteration will be '
            'interrupted. Set the value to 0 to skip the assembly correction step')
    correct_group.add_argument(
            '--correct_resolution', type=int, default=500,
            help='resolution for breakpoint detection and assembly correction, default: %(default)s (bp)')
    correct_group.add_argument(
            '--median_cov_ratio', type=float, default=0.2,
            help='for each contig, the coverage cutoff is calculated by multiplying the median coverage by this value, default: %(default)s. '
            'A continuous region with coverages >= the cutoff is considered as a "high-coverage region". Candidate breakpoints are determined in low-coverage regions '
            'bounded by two "large" high-coverage regions')
    correct_group.add_argument(
            '--region_len_ratio', type=float, default=0.1,
            help='a "large" high-coverage region should span a region longer than a threshold, which is the larger one of either the contig length * `region_len_ratio`, '
            'or `min_region_cutoff`, default: %(default)s')
    correct_group.add_argument(
            '--min_region_cutoff', type=int, default=5000,
            help='minimum length cutoff for a "large" high-coverage region, default: %(default)s (bp)')

    # Parameters for preprocessing (contig / link filtration) before clustering
    filter_group = parser.add_argument_group('>>> Parameters for preprocessing (contig / Hi-C link filtration) before clustering')
    filter_group.add_argument(
            '--Nx', type=int, default=80,
            help='keep contigs longer than Nx, default: %(default)s. This parameter speeds up the process and reduces memory consumption by reducing "n" in the '
            'Markov cluster algorithm, which has a time and space complexity of O(n^3) in its straightforward implementation. Additionally, fragmented '
            'contigs often fail to elongate due to a high proportion of repeats and assembly errors, so ignoring them in the clustering step could be better. '
            'The contigs filtered out by `--Nx` can be rescued in the following reassignment steps. To disable it, set this parameter to 100')
    filter_group.add_argument(
            '--RE_site_cutoff', type=int, default=5,
            help='minimum restriction enzyme sites, default: %(default)s. Contigs with fewer RE sites than this value will be filtered out')
    filter_group.add_argument(
            '--density_lower', default='0.2X',
            help='lower limit for contig Hi-C link density, used to filter out contigs that provide minimal linking information, default: %(default)s. '
            'This parameter works in two modes. [Fraction mode] If `density_lower` is a numeric value within the range [0, 1) , contigs with a link density '
            'less than (1 - `density_lower`) * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `density_lower` is a string '
            'ending with "X", contigs with a link density lower than `density_lower` times the median density will be filtered out. To disable it, set this '
            'parameter to 0')
    filter_group.add_argument(
            '--density_upper', default='1.9X',
            help='upper limit for contig Hi-C link density, used to filter out potential collapsed contigs, default: %(default)s. This parameter works '
            'in two modes. [Fraction mode] If `density_upper` is a numeric value whitin the range (0, 1], contigs with a link density greater than '
            '`density_upper` * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `density_upper` is a string ending with "X", '
            'contigs with a link density greater than `density_upper` times the median density will be filtered out. To disable it, set this parameter to 1')
    filter_group.add_argument(
            '--read_depth_upper', default='1.5X',
            help='upper limit for contig read depth, used to filter out potential collapsed contigs, default: %(default)s. This parameter works '
            'in two modes. [Fraction mode] If `read_depth_upper` is a numeric value within the range (0, 1], contigs with a read depth greater than '
            '`read_depth_upper` * 100%% of all contigs will be filtered out before Markov clustering. [Multiple mode] If `read_depth_upper` is a string ending with "X", '
            'contigs with a read depth greater than Q3 + `read_depth_upper` times the IQR (interquartile range) will be filtered out. To disable it, set this '
            'parameter to 1')
    filter_group.add_argument(
            '--topN', type=int, default=10,
            help='top N nearest neighbors (contigs) used for rank-sum calculation, default: %(default)s')
    filter_group.add_argument(
            '--rank_sum_hard_cutoff', type=int, default=0,
            help='hard filtering cutoff for rank-sum value, default: %(default)s (disabled). This parameter may be helpful if you are experienced in '
            'estimating the range of rank sum values, but it can also be dangerous. Be aware that the rank-sum values are affected by variables such as `topN` '
            'and contig length')
    filter_group.add_argument(
            '--rank_sum_upper', default='1.5X',
            help='upper limit for rank-sum value, used to filter out potential chimeric and collapsed contigs, default: %(default)s. This parameter works '
            'in two modes. [Fraction mode] If `rank_sum_upper` is a numeric value within the range (0, 1], contigs with a rank-sum value greater than '
            '`rank_sum_upper` * 100%% of all contigs will be filtered out before Markov clustering. If `rank_sum_upper` is a string ending with "X", '
            'contigs with a rank-sum value greater than Q3 + `rank_sum_upper` times the IQR (interquartile range) will be filtered out. To disable it, set this '
            'parameter to 1')
    filter_group.add_argument(
            '--remove_allelic_links', type=int, default=0,
            help='This parameter identifies allelic contig pairs and removes the Hi-C links between them. The value should be the ploidy and must be >= 2. '
            'By default, this function is disabled. It is designed to handle haplotype-phased assemblies with high switch error rates. Allelic contig pairs are '
            'identified based on the pattern of Hi-C links')
    filter_group.add_argument(
            '--concordance_ratio_cutoff', type=float, default=0.2,
            help='concordance ratio cutoff for allelic contig pair identification, default: %(default)s. In most cases, the default cutoff works well. Increasing '
            'the cutoff will increase both the true and false positive rates (TPR and FPR). This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--nwindows', type=int, default=50,
            help='this parameter is used in the concordance ratio calculation to eliminate the interference of contig length. The range (2 times the distance) '
            'for counting the "concordant Hi-C links" is dynamically defined by dividing the length of the shorter contig from a pair by `nwindows`, default: %(default)s')
    filter_group.add_argument(
            '--remove_concentrated_links', default=False, action='store_true',
            help='remove concentrated links that are concentrated in a tiny region of the genome with a remarkably higher link density, default: %(default)s')
    filter_group.add_argument(
            '--max_read_pairs', type=int, default=200,
            help='maximum Hi-C read pairs for allelic contig pair identification, default: %(default)s. This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--min_read_pairs', type=int, default=20,
            help='minimum Hi-C read pairs for allelic contig pair identification, default: %(default)s. This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--phasing_weight', type=float, default=1.0, 
            help='weight of phasing information, default: %(default)s. When this parameter is set to 1.0, all Hi-C links between haplotypes will be removed. Set this'
            'parameter to 0 to completely ignore the phasing information in the GFA files. This parameter only works if `--gfa` is set and mulitple gfa files are input')

    # Parameters for parsing ultra-long reads
    ul_group = parser.add_argument_group('>>> Parameters for parsing ultra-long reads')
    ul_group.add_argument(
            '--min_ul_mapq', type=int, default=30,
            help='MAPQ cutoff, alignments with MAPQ lower than this value will be removed, default: %(default)s')
    ul_group.add_argument(
            '--min_ul_alignment_length', type=int, default=10000,
            help='alignment length cutoff, alignments shorter than this value will be removed, default: %(default)s (bp)')
    ul_group.add_argument(
            '--max_distance_to_end', type=int, default=100,
            help='the distance to the ends of a contig should be less than this value, default: %(default)s (bp)')
    ul_group.add_argument(
            '--max_overlap_ratio', type=float, default=0.5,
            help='the length ratio of overlap between primary and supplementary alignments should be less than this value, default: %(default)s')
    ul_group.add_argument(
            '--max_gap_len', type=int, default=10000,
            help='maximum gap length between primary and supplementary alignments, default: %(default)s (bp)')
    ul_group.add_argument(
            '--min_ul_support', type=int, default=2,
            help='a linkage between two contigs should be supported by more than this number of UL reads, default: %(default)s')

    # Parameters for adjacency matrix construction and Markov clustering
    mcl_group = parser.add_argument_group('>>> Parameters for adjacency matrix construction and Markov Clustering')
    mcl_group.add_argument(
            '--bin_size', type=int, default=-1,
            help="if a contig's is greater than this value (unit: Kbp), the contig will be split into bins in this size "
            "to construct an djacency matrix. By default, `bin_size` is the smaller of either average chromosome length divided by 30, or 2000 kbp, "
            "but it will not be shorter than 100 Kbp. You can manually designate this value or set it to 0 to disable the function")
    mcl_group.add_argument(
            '--flank', type=int, default=500,
            help="consider only the Hi-C links from the flanking regions (ends) of each contig to construct the adjacency matrix, default: %(default)s (Kbp). "
            "If a contig's length is less than or equal to 2 times `flank`, all links will be used. Set this parameter to 0 to use the whole contigs")
    mcl_group.add_argument(
            '--normalize_by_nlinks', default=False, action='store_true',
            help='normalize inter-contig Hi-C links by the number of links to other contigs or groups, default: %(default)s')
    mcl_group.add_argument(
            '--expansion', type=int, default=2,
            help='expansion for Markov clustering, default: %(default)s')
    mcl_group.add_argument(
            '--min_inflation', type=float, default=1.1,
            help='minimum inflation for parameter tuning, default: %(default)s')
    mcl_group.add_argument(
            '--max_inflation', type=float, default=3.0,
            help='maximum inflation for parameter tuning, default: %(default)s')
    mcl_group.add_argument(
            '--inflation_step', type=float, default=0.1,
            help='step length of inflation for parameter tuning, default: %(default)s')
    mcl_group.add_argument(
            '--max_iter', type=int, default=200,
            help='maximum number of iterations of Markov clustering for each inflation, default: %(default)s')
    mcl_group.add_argument(
            '--pruning', type=float, default=0.0001,
            help='pruning threshold for Markov clustering, default: %(default)s')
    mcl_group.add_argument(
            '--skip_clustering', default=False, action='store_true',
            help='skip Markov clustering, default: %(default)s')

    # Parameters for performance
    performance_group = parser.add_argument_group('>>> Parameters for performance')
    performance_group.add_argument(
            '--threads', type=int, default=8,
            help='threads for reading BAM file, default: %(default)s')
    performance_group.add_argument(
            '--dense_matrix', default=False, action='store_true',
            help='the Markov Clustering is optimized with Scipy sparse matrix and Intel Math Kernal Library (MKL) by default. Set this parameter to '
            'disable it using dense matrix, default: %(default)s')

    # Parameters for logging
    logging_group = parser.add_argument_group('>>> Parameters for logging')
    logging_group.add_argument(
            '--verbose', default=False, action='store_true',
            help='verbose logging, default: %(default)s')

    args = parser.parse_args()

    # check parameters
    check_param('--density_lower', args.density_lower, {'X', 'x'})
    check_param('--density_upper', args.density_upper, {'X', 'x'})
    check_param('--read_depth_upper', args.read_depth_upper, {'X', 'x'})
    check_param('--rank_sum_upper', args.rank_sum_upper, {'X', 'x'})

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

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # make sure Intel MKL and sparse_dot_mkl have been correctly configured
    if not INTEL_MKL:
        logger.warning('Module sparse_dot_mkl or Intel MKL is not correctly installed, HapHiC will be executed in dense matrix mode')
        args.dense_matrix = True
    elif args.dense_matrix:
        logger.warning('--dense_matrix is set, HapHiC will be executed in dense matrix mode')

    # check file format for Hi-C read alignments
    if args.aln_format == 'auto':
        detect_format(args)

    if args.correct_nrounds and args.ul:
        args.ul = None
        logger.warning('Ultra-long data are not supported now when assembly correction is enabled')

    # quick view mode
    if args.quick_view:
        args.bin_size = 0
        args.Nx = 100
        args.remove_allelic_links = 0
        args.remove_concentrated_links = False

    # read draft genome in FASTA format,
    # construct a dict to store sequence and length of each contig
    fa_dict = parse_fasta(args.fasta, RE=args.RE)

    # parse gfa file(s) to get the read depth and phasing information
    if args.gfa:
        gfa_list = args.gfa.split(',')
        read_depth_dict = parse_gfa(gfa_list, fa_dict)
    else:
        read_depth_dict = dict()

    if args.correct_nrounds:
        # performing assembly correction
        if args.aln_format == 'bam':
            ctg_cov_dict, ctg_link_pos_dict = parse_bam_for_correction(fa_dict, args)
        else:
            assert args.aln_format in {'pairs', 'bgzipped_pairs'}
            ctg_cov_dict, ctg_link_pos_dict = parse_pairs_for_correction(fa_dict, args)
        nbroken_ctgs, final_break_pos_dict, final_break_frag_dict  = correct_assembly(
                ctg_cov_dict, ctg_link_pos_dict, fa_dict, read_depth_dict, args)
        del ctg_cov_dict, ctg_link_pos_dict
        gc.collect()

    # the original read depth info of contigs after correction
    ctg_read_depth_dict = read_depth_dict.copy()

    if args.ul:
        whitelist = set()
        path_list = parse_ul_alignments(args)
        for path in path_list:
            for i in range(len(path) - 1):
                if i % 2 != 0:
                    node1, node2 = path[i], path[i+1]
                    ctg1 = node1.rsplit('_', 1)[0]
                    ctg2 = node2.rsplit('_', 1)[0]
                    whitelist.add(ctg1)
                    whitelist.add(ctg2)
    else:
        whitelist = set()
    # # record whitelist for reassignment
    args.whitelist = whitelist

    # statistics for fragments (contigs / bins)
    _, bin_set, bin_size, frag_len_dict, Nx_frag_set, RE_site_dict, split_ctg_set = stat_fragments(
            fa_dict, args.RE, read_depth_dict, whitelist, nchrs=args.nchrs, flank=args.flank, Nx=args.Nx, bin_size=args.bin_size)

    # construct linking matrix from Hi-C alignments

    if args.correct_nrounds and nbroken_ctgs:
        # need all read1 Hi-C links and contig conversion
        format_options = [b'filter=flag.read1']
        if split_ctg_set:
            if args.aln_format == 'bam':
                alignments = bam_generator_for_correction(
                        args.alignments, args.threads, format_options, final_break_pos_dict, final_break_frag_dict)
            else:
                alignments = pairs_generator_for_correction(
                        args.alignments, args.aln_format, final_break_pos_dict, final_break_frag_dict)
        else:
            if args.aln_format == 'bam':
                alignments = bam_generator_for_correction_ctg(
                        args.alignments, args.threads, format_options, final_break_pos_dict, final_break_frag_dict)
            else:
                alignments = pairs_generator_for_correction_ctg(
                        args.alignments, args.aln_format, final_break_pos_dict, final_break_frag_dict)
    elif split_ctg_set:
        # need all read1 Hi-C links, but contig conversion is not necessary
        if args.aln_format == 'bam':
            format_options = [b'filter=flag.read1']
            alignments = bam_generator(args.alignments, args.threads, format_options)
        else:
            alignments = pairs_generator(args.alignments, args.aln_format)
    else:
        # need read1 inter-contig Hi-C links only, and contig conversion is not necessary
        if args.aln_format == 'bam':
            format_options = [b'filter=flag.read1 && refid != mrefid']
            alignments = bam_generator(args.alignments, args.threads, format_options)
        else:
            alignments = pairs_generator_inter_ctgs(args.alignments, args.aln_format)

    # parse alignments
    if split_ctg_set:
        full_link_dict, flank_link_dict, HT_link_dict, clm_dict, frag_link_dict, ctg_coord_dict, ctg_pair_to_frag = parse_alignments(
                alignments, fa_dict, args, bin_size, frag_len_dict, Nx_frag_set, split_ctg_set)
    else:
        full_link_dict, flank_link_dict, HT_link_dict, clm_dict, frag_link_dict, ctg_coord_dict = parse_alignments_for_ctgs(
                alignments, fa_dict, args, frag_len_dict, Nx_frag_set)

    if args.ul:    
        # update HT_link_dict based on the contig paths supported by ultra-long reads 
        add_HT_links_based_on_ul(path_list, HT_link_dict)

    # output HT_link_dict.pkl for HapHiC sort
    output_pickle(HT_link_dict, 'HT_link_dict', 'HT_links.pkl')
    del HT_link_dict
    gc.collect()

    if args.quick_view:
        finished_time = time.time()
        logger.info('Program finished in {}s'.format(finished_time-start_time))
        return None

    # output clm file
    output_clm(clm_dict)
    del clm_dict
    gc.collect()

    # normalize flank_link_dict
    if args.normalize_by_nlinks:
        normalize_by_nlinks(flank_link_dict, frag_link_dict)

    # remove concentrated links in full_link_dict:
    if args.remove_concentrated_links:
        for ctg_name_pair, data in ctg_coord_dict.items():
            if isinstance(data, list):
                full_link_dict[ctg_name_pair] *= data[1]

    # filter fragments based on RE sites, link densities and ranksums
    filtered_frags = filter_fragments(
            Nx_frag_set, RE_site_dict, args.RE_site_cutoff, frag_link_dict,
            args.density_lower, args.density_upper, args.topN, args.rank_sum_upper,
            args.rank_sum_hard_cutoff, flank_link_dict, read_depth_dict, args.read_depth_upper, whitelist)

    # update flank_link_dict & full_link_dict by removing Hi-C links between alleic contig pairs
    if args.remove_allelic_links:
        if split_ctg_set:
            filtered_frags = remove_allelic_HiC_links(
                    fa_dict, ctg_coord_dict, full_link_dict, args, flank_link_dict, filtered_frags, ctg_pair_to_frag)
        else:
            filtered_frags = remove_allelic_HiC_links(
                    fa_dict, ctg_coord_dict, full_link_dict, args, flank_link_dict, filtered_frags)
        del ctg_coord_dict
        gc.collect()

    # update flank_link_dict and full_link_dict based on the contig paths supported by ultra-long reads
    if args.ul:
        add_flank_and_full_links_based_on_ul(path_list, flank_link_dict, full_link_dict, bin_set)

    # update flank_link_dict & full_link_dict based on the phasing information
    if args.gfa and len(gfa_list) >= 2 and args.phasing_weight:
        reduce_inter_hap_HiC_links(flank_link_dict, read_depth_dict, args.phasing_weight, target='flank_link_dict')
        reduce_inter_hap_HiC_links(full_link_dict, ctg_read_depth_dict, args.phasing_weight, target='full_link_dict')

    # output full_link_dict in pickle format after update
    output_pickle(full_link_dict, 'full_link_dict', 'full_links.pkl')

    # convert flank_link_dict to matrix (for clustering)
    flank_link_matrix, frag_index_dict = dict_to_matrix(
            flank_link_dict, filtered_frags, dense_matrix=args.dense_matrix, add_self_loops=True)

    del flank_link_dict
    gc.collect()

    matrix_time = time.time()
    logger.info('Hi-C linking matrix was constructed in {}s'.format(matrix_time-start_time))

    if not args.skip_clustering:
        # run Markov clustering, and output results
        result_clusters_list, mcl_nrounds = run_mcl_clustering(
                flank_link_matrix, bin_set, frag_len_dict, frag_index_dict, args.expansion,
                args.min_inflation, args.max_inflation, args.inflation_step,
                args.max_iter, args.pruning, fa_dict, args.nchrs, args.dense_matrix
                )

        clustering_time = time.time()
        logger.info('{} round(s) of Markov clustering finished in {}s, average {}s per round'.format(
            mcl_nrounds, clustering_time-matrix_time, (clustering_time-matrix_time)/mcl_nrounds))

        # statistics for RE sites, Hi-C links, link density, link density ratio
        output_statistics(fa_dict, full_link_dict, result_clusters_list)

    finished_time = time.time()
    logger.info('Program finished in {}s'.format(finished_time-start_time))


def main():

    # get arguments
    args = parse_arguments()

    run(args, 'HapHiC_cluster.log')


if __name__ == '__main__':
    main()

