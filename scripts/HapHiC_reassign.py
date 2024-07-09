#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Email: zengxf@sustech.edu.cn

import sys
import os
import pickle
import argparse
import logging
import time
import pysam
import gzip

from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from HapHiC_cluster import parse_fasta, parse_gfa, stat_fragments, remove_allelic_HiC_links, dict_to_matrix

from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

pysam.set_verbosity(0)


def parse_pickle(fa_dict, pickle_file):

    logger.info('Parsing input pickle file...')

    with open(pickle_file, 'rb') as f:
        full_link_dict = pickle.load(f)

    # sort by contig length
    sorted_ctg_list = sorted([(ctg, fa_dict[ctg][1]) for ctg in fa_dict], key=lambda x: x[1], reverse=True)

    RE_site_dict = {ctg: ctg_info[2] for ctg, ctg_info in fa_dict.items()}

    return full_link_dict, sorted_ctg_list, RE_site_dict


def parse_pairs_for_reassignment(fa_dict, args, logger=logger):

    """parse pairs file (for reassignment)"""
    logger.info('Parsing input pairs file...')

    # set default parameters
    pairs = args.links

    full_link_dict = defaultdict(int)

    # a dict storing coords of Hi-C links between contigs
    ctg_coord_dict = defaultdict(lambda: array('i'))

    if pairs.endswith('.pairs'):
        fopen = open
    else:
        assert pairs.endswith('.pairs.gz')
        fopen = gzip.open

    with fopen(pairs, 'rt') as f:

        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            ref, mref = cols[1], cols[3]

            # skip intra-contig Hi-C links and contigs not in the fasta file
            if ref == mref or ref not in fa_dict or mref not in fa_dict:
                continue

            # sort by contig name
            (ctg_i, coord_i), (ctg_j, coord_j) = sorted(((ref, int(cols[2])), (mref, int(cols[4]))))
            ctg_name_pair = (ctg_i, ctg_j)

            # update full_link_dict
            full_link_dict[ctg_name_pair] += 1

            # record coord pairs and calculate concordance ratios once there are enough coord pairs
            if args.remove_allelic_links:
                record_coord_pairs(ctg_coord_dict, ctg_name_pair, coord_i, coord_j, args.max_read_pairs, fa_dict, args.nwindows)

    return full_link_dict, ctg_coord_dict


def parse_bam_for_reassignment(fa_dict, args, logger=logger):

    """parse BAM file (for reassignment)"""
    logger.info('Parsing input BAM file...')

    # set default parameters
    bam = args.links

    full_link_dict = defaultdict(int)

    # a dict storing coords of Hi-C links between contigs
    ctg_coord_dict = defaultdict(lambda: array('i'))

    format_options = [b'filter=flag.read1 && refid != mrefid']

    with pysam.AlignmentFile(bam, mode='rb', threads=args.threads, format_options=format_options) as f:

        for aln in f:

            # may be slightly faster
            ref, mref = aln.reference_name, aln.next_reference_name

            # skip contigs not in the fasta file
            if ref not in fa_dict or mref not in fa_dict:
                continue

            # sort by contig name
            (ctg_i, coord_i), (ctg_j, coord_j) = sorted(((ref, aln.reference_start+1), (mref, aln.next_reference_start+1)))
            ctg_name_pair = (ctg_i, ctg_j)

            # update full_link_dict
            full_link_dict[ctg_name_pair] += 1

            # record coord pairs and calculate concordance ratios once there are enough coord pairs
            if args.remove_allelic_links:
                record_coord_pairs(ctg_coord_dict, ctg_name_pair, coord_i, coord_j, args.max_read_pairs, fa_dict, args.nwindows)

    return full_link_dict, ctg_coord_dict


def parse_clusters(clusters_file, RE_site_dict, fa_dict, min_group_len):

    logger.info('Parsing .clusters.txt file...')

    ctg_group_dict = dict()
    group_RE_dict = dict()

    with open(clusters_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            group = cols[0]
            # dismiss groups smaller than min_group_len
            if min_group_len and sum([fa_dict[ctg][1] for ctg in cols[2:]]) / 1000000 < min_group_len:
                continue
            # add pseudo-count of 1
            group_RE_dict[group] = 1
            for ctg in cols[2:]:
                ctg_group_dict[ctg] = group
                # minus the pseudo-count of 1 for each contig
                group_RE_dict[group] += RE_site_dict[ctg] - 1

    return ctg_group_dict, group_RE_dict


def parse_assembly(assembly_file, RE_site_dict, fa_dict, min_group_len):

    logger.info('Parsing .assembly file...')

    ctg_dict = dict()
    ctg_group_dict = dict()
    group_RE_dict = dict()
    n = 0
    with open(assembly_file) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if line.startswith('>'):
                ctg_dict[cols[1]] = cols[0][1:]
            else:
                n += 1
                group = 'group{}'.format(n)
                # dismiss groups smaller than min_group_len
                if min_group_len and sum([fa_dict[ctg_dict[num.strip('-')]][1] for num in cols]) / 1000000 < min_group_len:
                    continue
                # add pseudo-count of 1
                group_RE_dict[group] = 1
                for num in cols:
                    ctg = ctg_dict[num.strip('-')]
                    ctg_group_dict[ctg] = group
                    # minus the pseudo-count of 1 for each contig
                    group_RE_dict[group] += RE_site_dict[ctg] - 1

    return ctg_group_dict, group_RE_dict


def add_ungrouped_ctgs(fa_dict, ctg_group_dict):

    # In the agglomerative hierarchical clustering step, we only consider the filtered
    # (high-confidence) contigs used in the Markov clustering ('grouped', else 'un-').

    grouped_ctgs = set()
    for ctg in fa_dict:
        if ctg not in ctg_group_dict:
            ctg_group_dict[ctg] = 'ungrouped'
        else:
            grouped_ctgs.add(ctg)

    return grouped_ctgs


def parse_link_dict(link_dict, ctg_group_dict, normalize_by_nlinks=False):

    def add_ctg_group(ctg, group, links):
        if group != 'ungrouped':
            if group in ctg_group_link_dict[ctg]:
                ctg_group_link_dict[ctg][group] += links
            else:
                ctg_group_link_dict[ctg][group] = links

    ctg_group_link_dict = defaultdict(dict)
    linked_ctg_dict = defaultdict(set)
    ctg_link_dict = defaultdict(int)

    if normalize_by_nlinks:
        total_links = 0
        normalized_total_links = 0
        for (ctg_i, ctg_j), links in link_dict.items():
            ctg_link_dict[ctg_i] += links
            ctg_link_dict[ctg_j] += links
            total_links += links

    for (ctg_i, ctg_j), links in link_dict.items():
        if normalize_by_nlinks:
            # normalize by geometric mean
            links /= (ctg_link_dict[ctg_i] * ctg_link_dict[ctg_j]) ** 0.5
            link_dict[(ctg_i, ctg_j)] = links
            normalized_total_links += links
        else:
            group_i, group_j = ctg_group_dict[ctg_i], ctg_group_dict[ctg_j]
            add_ctg_group(ctg_i, group_j, links)
            add_ctg_group(ctg_j, group_i, links)
            linked_ctg_dict[ctg_i].add(ctg_j)
            linked_ctg_dict[ctg_j].add(ctg_i)

    # make sure the total links are the same
    if normalize_by_nlinks:
        scale_factor = total_links / normalized_total_links
        for (ctg_i, ctg_j), links in link_dict.items():
            links *= scale_factor
            link_dict[(ctg_i, ctg_j)] = links
            group_i, group_j = ctg_group_dict[ctg_i], ctg_group_dict[ctg_j]
            add_ctg_group(ctg_i, group_j, links)
            add_ctg_group(ctg_j, group_i, links)
            linked_ctg_dict[ctg_i].add(ctg_j)
            linked_ctg_dict[ctg_j].add(ctg_i)

    return ctg_group_link_dict, linked_ctg_dict


def run_reassignment(
        sorted_ctg_list, ctg_group_link_dict, ctg_group_dict, full_link_dict, linked_ctg_dict, fa_dict, RE_site_dict, gfa,
        group_RE_dict, max_ctg_len, min_RE_sites, min_links, min_link_density, min_density_ratio, ambiguous_cutoff, min_group_len, whitelist, nround):

    def update(ctg, max_group):

        # reassign new group (max_group) to the ctg
        former_group = ctg_group_dict[ctg]
        ctg_group_dict[ctg] = max_group
        # update ctg_group_link_dict
        for each_ctg in linked_ctg_dict[ctg]:
            group_links = ctg_group_link_dict[each_ctg]
            ctg_name_pair = tuple(sorted([ctg, each_ctg]))
            links = full_link_dict[ctg_name_pair]
            if former_group != 'ungrouped':
                group_links[former_group] -= links
            if max_group in group_links:
                group_links[max_group] += links
            elif max_group != 'ungrouped':
                group_links[max_group] = links


    def cal_link_density(ctg, group, former_group, links):

        group_RE_sites = group_RE_dict[group]
        if group == former_group:
            return links / group_RE_sites
        else:
            return links / (group_RE_sites + RE_site_dict[ctg] - 1)

    if nround:
        logger.info('Performing reassignment...')
        round_name = 'round{}'.format(nround)
    else:
        logger.info('Performing additional round of rescue...')
        round_name = 'additional_rescue'

    # dismiss groups smaller than min_group_len
    if min_group_len and nround > 1:

        # calculate group length
        group_len_dict = defaultdict(int)
        for ctg, group in ctg_group_dict.items():
            if group != 'ungrouped':
                group_len_dict[group] += fa_dict[ctg][1]

        # get dismissed groups
        dismissed_group_set = set()
        for group, group_len in group_len_dict.items():
            if group_len / 1000000 < min_group_len:
                dismissed_group_set.add(group)

        for ctg, group in ctg_group_dict.items():
            # ctgs in dismissed groups are assigned "ungrouped"
            if group in dismissed_group_set:
                ctg_group_dict[ctg] = 'ungrouped'
            # links to dismissed groups are set to 0
            for group in dismissed_group_set:
                ctg_group_link_dict[ctg][group] = 0

    # traverse all contigs from longer to shorter
    result_dict = defaultdict(int)
    for ctg, ctg_len in sorted_ctg_list:

        former_group = ctg_group_dict[ctg]
        group_links = ctg_group_link_dict[ctg]
        filtered = False

        # filter RE_site
        if (RE_site_dict[ctg] - 1 < min_RE_sites and ctg not in whitelist) or not group_links:
            max_group, max_links = 'ungrouped', 0
            logger.debug('[reassignment::{}] Contig {} is not rescued (too few RE sites, {})'.format(round_name, ctg, former_group))
            result_dict['not_rescued'] += 1
            filtered = True
        else:
            sorted_group_links = sorted(group_links.items(), key=lambda x: x[1], reverse=True)
            max_group, max_links = sorted_group_links[0]
            if len(sorted_group_links) > 1:
                second_links = sorted_group_links[1][1]
            else:
                second_links = 0
            # filter max_links
            if max_links < min_links and ctg not in whitelist:
                max_group, max_links = 'ungrouped', 0
                logger.debug('[reassignment::{}] Contig {} is not rescued (too few Hi-C links, {})'.format(round_name, ctg, former_group))
                result_dict['not_rescued'] += 1
                filtered = True
            # we define contigs with second_links / max_links >= ambiguous_cutoff as ambiguous,
            # these contigs will NOT be reassigned, but only rescued in the additional round of rescue
            elif nround and second_links / max_links >= ambiguous_cutoff and ctg not in whitelist:
                max_group, max_links = 'ungrouped', 0
                logger.debug('[reassignment::{}] Contig {} is not rescued (ambiguous, {})'.format(round_name, ctg, former_group))
                result_dict['not_rescued'] += 1
                filtered = True
            else:
                # calculate link densities
                max_link_density = cal_link_density(ctg, max_group, former_group, max_links)
                # filter max_link_density
                if max_link_density < min_link_density and ctg not in whitelist:
                    max_group, max_links = 'ungrouped', 0
                    logger.debug('[reassignment::{}] Contig {} is not rescued (too low Hi-C links density, {})'.format(round_name, ctg, former_group))
                    result_dict['not_rescued'] += 1
                    filtered = True
                else:
                    # calculate the average link density to other nonbest groups
                    # other_group_density_sum = sum([cal_link_density(ctg, group, former_group, links) for group, links in sorted_group_links[1:]])

                    # a patch for gfa file
                    if gfa:
                        other_group_density_sum = sum([cal_link_density(ctg, group, former_group, links) for group, links in sorted_group_links[1:] if links])
                    else:
                        other_group_density_sum = sum([cal_link_density(ctg, group, former_group, links) for group, links in sorted_group_links[1:]])

                    if other_group_density_sum:
                        if gfa:
                            avg_other_group_density = other_group_density_sum / len([group for group, links in sorted_group_links[1:] if links])
                        else:
                            avg_other_group_density = other_group_density_sum / (len(group_RE_dict) - 1)
                    else:
                        # assign a big number to prevent division by zero problem
                        avg_other_group_density = 1000000000

        if filtered:
            continue

        # if best group is "ungrouped", keep consistent
        if max_group == 'ungrouped':
            assert False
            # consistent
            logger.debug('[reassignment::{}] Contig {} is consistent ({})'.format(round_name, ctg, former_group))
            result_dict['consistent'] += 1
        else:
            # rescue
            if former_group == 'ungrouped':
                # rescued
                if max_link_density / avg_other_group_density >= min_density_ratio:
                    update(ctg, max_group)
                    group_RE_dict[max_group] += RE_site_dict[ctg] - 1
                    logger.debug('[reassignment::{}] Contig {} is rescued (ungrouped -> {})'.format(round_name, ctg, max_group))
                    result_dict['rescued'] += 1
                # not rescued
                else:
                    logger.debug('[reassignment::{}] Contig {} is not rescued (too low density ratio)'.format(round_name, ctg))
                    result_dict['not_rescued'] += 1
            # consistent
            elif former_group in group_links and group_links[former_group] == max_links:
                logger.debug('[reassignment::{}] Contig {} is consistent ({})'.format(round_name, ctg, former_group))
                result_dict['consistent'] += 1
            # reassignment, contigs longer than max_ctg_len * 1000 will NOT be reassgined
            elif nround and ctg_len <= max_ctg_len * 1000 and max_link_density / avg_other_group_density >= min_density_ratio:
                update(ctg, max_group)
                group_RE_dict[former_group] -= RE_site_dict[ctg] - 1
                group_RE_dict[max_group] += RE_site_dict[ctg] - 1
                logger.debug('[reassignment::{}] Contig {} is reassigned ({} -> {})'.format(round_name, ctg, former_group, max_group))
                result_dict['reassigned'] += 1
            # consistent
            else:
                logger.debug('[reassignment::{}] Contig {} is consistent ({})'.format(round_name, ctg, former_group))
                result_dict['consistent'] += 1

    logger.info('[result::{}] Total: {}, consistent: {}, rescued: {}, reassigned: {}, not rescued: {}'.format(
        round_name, len(sorted_ctg_list), result_dict['consistent'], result_dict['rescued'], result_dict['reassigned'], result_dict['not_rescued']))


def stat_clusters(ctg_group_dict, fa_dict, grouped_ctgs):

    group_ctg_dict = dict()
    group_hiconf_RE_dict = dict()

    for ctg, group in ctg_group_dict.items():
        if group == 'ungrouped':
            continue
        ctg_len = fa_dict[ctg][1]
        # 'RE_sites' is used to store the RE site sum of high-confidence contigs for the
        # link density calculation in the agglomerative hierarchical clustering step.
        RE_sites = fa_dict[ctg][2] - 1 if ctg in grouped_ctgs else 0
        if group in group_ctg_dict:
            group_ctg_dict[group][0].add((ctg, ctg_len))
            group_ctg_dict[group][1] += ctg_len
            group_hiconf_RE_dict[group] += RE_sites
        else:
            group_ctg_dict[group] = [{(ctg, ctg_len)}, ctg_len]
            group_hiconf_RE_dict[group] = RE_sites

    return group_ctg_dict, group_hiconf_RE_dict


def clusters_output(group_ctg_dict, fa_dict, out_prefix):

    new_ctg_group_dict = dict()
    new_old_group_dict = dict()
    new_group_ctg_dict = dict()

    sorted_clusters = sorted(group_ctg_dict.items(), key=lambda x: x[1][1], reverse=True)

    if out_prefix == 'reassigned':
        subdir = 'reassigned_groups'
    else:
        assert out_prefix == 'hc'
        subdir = 'hc_groups'

    with open('{}/{}_clusters.txt'.format(subdir, out_prefix), 'w') as fclusters:

        fclusters.write('#Group\tnContigs\tContigs\n')
        for n, (group, ctg_stat) in enumerate(sorted_clusters, 1):

            group_name = 'group{}_{}bp'.format(n, ctg_stat[1])
            new_old_group_dict[group_name] = group
            new_group_ctg_dict[group_name] = ctg_stat

            sorted_ctgs = [ctg for ctg, ctg_len in sorted(list(ctg_stat[0]), key=lambda x: x[1], reverse=True)]

            fclusters.write('{}\t{}\t{}\n'.format(group_name, len(ctg_stat[0]), ' '.join(sorted_ctgs)))

            with open('{}/{}_{}.txt'.format(subdir, out_prefix, group_name), 'w') as fgroup:
                fgroup.write('#Contig\tRECounts\tLength\n')
                for ctg in sorted_ctgs:
                    fgroup.write('{}\t{}\t{}\n'.format(ctg, fa_dict[ctg][2], fa_dict[ctg][1]))
                    new_ctg_group_dict[ctg] = group_name

    return new_ctg_group_dict, new_group_ctg_dict, new_old_group_dict


def agglomerative_hierarchical_clustering(full_link_dict, grouped_ctgs, new_ctg_group_dict, group_hiconf_RE_dict, new_old_group_dict, nclusters, normalize_by_nlinks=False):

    logger.info('Performing additional agglomerative hierarchical clustering...')

    # construct a Hi-C linking matrix between clusters
    group_link_dict = defaultdict(int)

    for (ctg_i, ctg_j), links in full_link_dict.items():
        if ctg_i not in grouped_ctgs or ctg_j not in grouped_ctgs:
            continue
        if ctg_i not in new_ctg_group_dict or ctg_j not in new_ctg_group_dict:
            continue
        group_i, group_j = new_ctg_group_dict[ctg_i], new_ctg_group_dict[ctg_j]
        if group_i == group_j:
            continue

        group_name_pair = tuple(sorted([group_i, group_j]))

        # group_link_dict only records links between high-confidence contigs from different groups
        group_link_dict[group_name_pair] += links

    # calculate link density
    with open('hc_groups/group_group_links.txt', 'w') as f:

        f.write('group1\tgroup2\tlinks\tlink_density\n')
        max_link_density = 0

        if normalize_by_nlinks:
            group_group_link_dict = defaultdict(int)
            for (group_i, group_j), links in group_link_dict.items():
                group_group_link_dict[group_i] += links
                group_group_link_dict[group_j] += links

        for group_name_pair, links in group_link_dict.items():
            group_i, group_j = group_name_pair
            group_RE_sites_i = group_hiconf_RE_dict[new_old_group_dict[group_i]]
            group_RE_sites_j = group_hiconf_RE_dict[new_old_group_dict[group_j]]
            if normalize_by_nlinks:
                link_density = links / (group_group_link_dict[group_i] * group_group_link_dict[group_j])
            else:
                link_density = links / (group_RE_sites_i * group_RE_sites_j)
            max_link_density = max(max_link_density, link_density)
            group_link_dict[group_name_pair] = link_density
            f.write('{}\t{}\t{}\t{}\n'.format(group_i, group_j, links, link_density))

    # for (group_i, group_j), links in group_link_dict.items():
    #     print(group_i, group_j, links)

    # convert group_link_dict to a matrix for agglomerative hierarchical clustering
    group_link_matrix, group_index_dict = dict_to_matrix(group_link_dict, set(new_old_group_dict.keys()))
    index_group_dict = {i : g for g, i in group_index_dict.items()}

    # use max - density to construct distance matrix
    group_link_matrix = max_link_density - group_link_matrix

    if 'affinity' in AgglomerativeClustering._get_param_names():
        clust = AgglomerativeClustering(n_clusters=nclusters, affinity="precomputed", linkage="average", distance_threshold=None)
    else:
        # add support for higher version of scikit-learn
        logger.info('The parameter "affinity" is not found in AgglomerativeClustering, use "metric" instead')
        assert 'metric' in AgglomerativeClustering._get_param_names()
        clust = AgglomerativeClustering(n_clusters=nclusters, metric="precomputed", linkage="average", distance_threshold=None)

    clusters = clust.fit_predict(group_link_matrix)

    hc_cluster_dict = defaultdict(list)

    for index, new_cluster in enumerate(clusters):
        old_cluster = index_group_dict[index]
        hc_cluster_dict[new_cluster].append(old_cluster)

    return hc_cluster_dict


def stat_hc_clusters(group_ctg_dict, hc_cluster_dict):

    # output Agglomerative Hierarchical Clustering result
    hc_group_ctg_dict = dict()

    with open('hc_groups/hc_result.txt', 'w') as fresult:
        fresult.write('hc_id\treassigned_groups\n')
        for new_cluster, old_clusters in hc_cluster_dict.items():
            fresult.write('{}\t{}\n'.format(new_cluster, ' '.join(old_clusters)))
            hc_group_ctg_dict[new_cluster] = [set(), 0]

            for old_cluster in old_clusters:
                hc_group_ctg_dict[new_cluster][0] |= group_ctg_dict[old_cluster][0]
                hc_group_ctg_dict[new_cluster][1] += group_ctg_dict[old_cluster][1]

    return hc_group_ctg_dict


def split_clm_file(clm_file, group_ctg_dict, ctg_group_dict, subdir):

    logger.info('Splitting clm file into subfiles by group...')

    # make directory for final groups
    final_dir = 'final_groups'
    os.mkdir(final_dir)

    # create symbolic links for final groups
    if subdir == 'reassigned_groups':
        prefix = 'reassigned'
    else:
        assert subdir == 'hc_groups'
        prefix = 'hc'

    for group in group_ctg_dict:
        # group files
        os.symlink('../{0}/{1}_{2}.txt'.format(subdir, prefix, group),
                   '{0}/{1}.txt'.format(final_dir, group))

    # clusters file
    os.symlink('../{0}/{1}_clusters.txt'.format(subdir, prefix),
               '{0}/final_clusters.txt'.format(final_dir))

    # make directory for clm splitting
    subdir = 'split_clms'
    os.mkdir(subdir)

    fp_dict = dict()

    for group in group_ctg_dict:
        fp_dict[group] = open('{}/{}.clm'.format(subdir, group), 'w')

    with open(clm_file) as f:
        for line in f:
            cols = line.split()
            ctg_1, ctg_2 = cols[0][:-1], cols[1][:-1]
            if ctg_1 in ctg_group_dict and ctg_2 in ctg_group_dict and ctg_group_dict[ctg_1] == ctg_group_dict[ctg_2]:
                fp_dict[ctg_group_dict[ctg_1]].write(line)

    for group, fp in fp_dict.items():
        fp.close()


def mock_clusters_file(fa_dicts, total_lens, final_dir):

    with open('{}/final_clusters.txt'.format(final_dir), 'w') as f:
        f.write('#Group\tnContigs\tContigs\n')
        for n, fa_dict in enumerate(fa_dicts):
            total_len = total_lens[n]
            f.write('group{}_{}bp\t{}\t{}\n'.format(n+1, total_len, len(fa_dict), ' '.join([ctg for ctg in fa_dict])))


def mock_group_file(fa_dicts, total_lens, final_dir):

    for n, fa_dict in enumerate(fa_dicts):
        total_len = total_lens[n]
        with open('{}/group{}_{}bp.txt'.format(final_dir, n+1, total_len), 'w') as f:
            f.write('#Contig\tRECounts\tLength\n')
            for ctg, ctg_info in fa_dict.items():
                f.write('{}\t{}\t{}\n'.format(ctg, ctg_info[2], ctg_info[1]))


def parse_arguments():

    parser = argparse.ArgumentParser(prog='haphic reassign')

    # Parameters for parsing input files and pipeline control
    input_group = parser.add_argument_group('>>> Parameters for parsing input files and pipeline control')
    input_group.add_argument(
            'fasta', help='draft genome in FASTA format. Use `corrected_asm.fa` generated in the clustering step when `--correct_nrounds` was set')
    input_group.add_argument(
            'links', help='`full_links.pkl` generated in the clustering step (much faster) or filtered Hi-C read alignments in BAM/pairs format (DO NOT sort it by coordinate)')
    input_group.add_argument(
            'clusters', help='`*.clusters.txt` file or `*.assembly` file')
    input_group.add_argument(
            'clm', help='`paired_links.clm` generated in the clustering step. After reassignment, the clm file will be split into subfiles by group for the following '
            'sorting (ordering and orientation) step')
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
            help='(experimental) GFA file(s) of the phased hifiasm assembly, separated with commas (--gfa "gfa1,gfa2"), default: %(default)s. In the reassignment step, '
            'this parameter only works with `--quick_view` to devide contigs into different haplotype groups based on the phasing information in the GFA files')

    # Parameters for reassignment and rescue
    reassign_group = parser.add_argument_group('>>> Parameters for reassignment and rescue')
    reassign_group.add_argument(
            '--min_group_len', type=float, default=5,
            help='minimum group length, default: %(default)s (Mbp). Groups smaller than this value will be dismissed and the contigs inside will be rescued and reassigned')
    reassign_group.add_argument(
            '--max_ctg_len', type=float, default=10000,
            help='maximum contig length for reassignment, default: %(default)s (Kbp). Contigs longer than this value will NOT be reassigned as they are error-prone. '
            'Set this parameter to 0 to rescue only, and no contig will be reassigned')
    reassign_group.add_argument(
            '--min_RE_sites', type=int, default=25,
            help='minimum restriction enzyme sites, default: %(default)s. Contigs with fewer RE sites than this value will be moved to "ungrouped"')
    reassign_group.add_argument(
            '--min_links', type=int, default=25,
            help='minimum Hi-C links, default: %(default)s. Contigs with fewer Hi-C links than this value will be moved to "ungrouped"')
    reassign_group.add_argument(
            '--min_link_density', type=float, default=0.0001,
            help='minimum Hi-C link density (= links / RE_sites), default: %(default)s. Contigs with Hi-C link density less than this value will be moved to "ungrouped"')
    reassign_group.add_argument(
            '--min_density_ratio', type=float, default=4,
            help='minimum Hi-C link density ratio (= link density to the best group / to the average of other groups), default: %(default)s. Contigs with link density ratio '
            'less than this value will NOT be reassigned and rescued')
    reassign_group.add_argument(
            '--ambiguous_cutoff', type=float, default=0.6,
            help='contigs with Hi-C links to the second-best group / to the best group >= this value are defined as ambiguous contigs, default: %(default)s. '
            'These contigs will NOT be reassigned, but only rescued in an additional round of rescue.')
    reassign_group.add_argument(
            '--reassign_nrounds', type=int, default=5,
            help='maximum number of rounds for reassignment, default: %(default)s. If the result converges, the iterations will be interrupted')
    reassign_group.add_argument(
            '--normalize_by_nlinks', default=False, action='store_true',
            help="normalize inter-contig and inter-group Hi-C links by the number of links to other contigs or groups, default: %(default)s")
    reassign_group.add_argument(
            '--nclusters', type=int, default=0,
            help='perform additional agglomerative hierarchical clustering to concatenate clusters after reassignment and rescue, '
            'the value is often set to the expected number of chromosomes to get chromosome-level scaffolds. Set the parameter to 0 to disable it, default: %(default)s')
    reassign_group.add_argument(
        '--no_additional_rescue', default=False, action='store_true',
        help='do not run the additional round of rescue, default: %(default)s')

    # Parameters for Hi-C link filtration
    filter_group = parser.add_argument_group('>>> Parameters for Hi-C link filtration before clustering, work only if the input "links" file is in BAM or pairs format')
    filter_group.add_argument(
            '--remove_allelic_links', type=int, default=0,
            help='This parameter identifies allelic contig pairs and removes the Hi-C links between them, the value should be the ploidy and must be >= 2. '
            'By default, this function is disabled. It is designed to handle haplotype-phased assemblies with high switch error rates. Allelic contig pairs are '
            'identified based on the pattern of Hi-C links')
    filter_group.add_argument(
            '--concordance_ratio_cutoff', type=float, default=0.2,
            help='concordance ratio cutoff for allelic contig pair identification, default: %(default)s. In most cases, the default cutoff works well. Increasing'
            'the cutoff will increase both the true and false positive rates (TPR and FPR). This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--nwindows', type=int, default=50,
            help='this parameter is used in the concordance ratio calculation to eliminate the interference of contig length. The range (2 times the distance) '
            'for counting the "concordant Hi-C links" is dynamically defined by dividing the length of the shorter contig from a pair by `nwindows`, default: %(default)s')
    filter_group.add_argument(
            '--max_read_pairs', type=int, default=200,
            help='maximum Hi-C read pairs for allelic contig pair identification, default: %(default)s. This parameter only works if `--remove_allelic_links` is set')
    filter_group.add_argument(
            '--min_read_pairs', type=int, default=20,
            help='minimum Hi-C read pairs for allelic contig pair identification, default: %(default)s. This parameter only works if `--remove_allelic_links` is set')

    # Parameters for performance
    performance_group = parser.add_argument_group('>>> Parameters for performance')
    performance_group.add_argument(
            '--threads', type=int, default=8,
            help='threads for reading BAM file, default: %(default)s')

    # Parameters for logging
    logging_group = parser.add_argument_group('>>> Parameters for logging')
    logging_group.add_argument(
            '--verbose', default=False, action='store_true',
            help='verbose logging, default: %(default)s')

    args = parser.parse_args()

    # check parameters
    if not (args.links.endswith('.bam') or args.links.endswith('.pkl') or args.links.endswith('.pairs') or args.links.endswith('.pairs.gz')):
        logger.error('The second positional argument "links" should end with .bam, .pkl, .pairs, or .pairs.gz')
        raise RuntimeError('Parameter check failed')

    if not args.clusters.endswith('.clusters.txt') and not args.clusters.endswith('.assembly'):
        logger.error('The third positional argument "clusters" should end with .clusters.txt or .assembly')
        raise RuntimeError('Parameter check failed')

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

    # read draft genome in FASTA format,
    # construct a dict to store sequence and length of each contig / scaffold
    fa_dict = parse_fasta(args.fasta, RE=args.RE, logger=logger)

    # quick view mode
    if args.quick_view:
        if args.gfa:
            gfa_list = args.gfa.split(',')
            read_depth_dict = parse_gfa(gfa_list, fa_dict)
        else:
            gfa_list = list()
            read_depth_dict = dict()

        final_dir = 'final_groups'
        os.mkdir(final_dir)
        if len(gfa_list) <= 1:
            total_len = sum([ctg_info[1] for ctg_info in fa_dict.values()])
            mock_clusters_file((fa_dict,), (total_len,), final_dir)
            mock_group_file((fa_dict,), (total_len,), final_dir)
        else:
            hap_ctg_dict = defaultdict(set)
            for ctg, (hap, _) in read_depth_dict.items():
                if ctg in fa_dict:
                    hap_ctg_dict[hap].add(ctg)
            fa_dicts, total_lens = list(), list()
            for hap, ctgs in hap_ctg_dict.items():
                hap_fa_dict = dict()
                for ctg in ctgs:
                    hap_fa_dict[ctg] = fa_dict[ctg]
                fa_dicts.append(hap_fa_dict)
                total_lens.append(sum([ctg_info[1] for ctg_info in hap_fa_dict.values()]))
            mock_clusters_file(fa_dicts, total_lens, final_dir)
            mock_group_file(fa_dicts, total_lens, final_dir)

        finished_time = time.time()
        logger.info('Program finished in {}s'.format(finished_time-start_time))
        return None

    if args.links.endswith('.pkl'):
        if args.remove_allelic_links:
            logger.info("A pickle file is input meanwhile --remove_allelic_links is set, "
                    "allelic Hi-C links will NOT be treated specially in the reassignment step")
        full_link_dict, sorted_ctg_list, RE_site_dict = parse_pickle(fa_dict, args.links)

    else:

        stat_output = stat_fragments(fa_dict, args.RE, {}, set(), logger=logger)
        sorted_ctg_list = stat_output[0]
        RE_site_dict = stat_output[-2]

        # parse bam
        if args.links.endswith('.bam'):
            full_link_dict, ctg_coord_dict = parse_bam_for_reassignment(fa_dict, args, logger=logger)
        # parse pairs
        else:
            full_link_dict, ctg_coord_dict = parse_pairs_for_reassignment(fa_dict, args, logger=logger)


        # update full_link_dict by removing Hi-C links between alleic contig pairs
        if args.remove_allelic_links:
            remove_allelic_HiC_links(fa_dict, ctg_coord_dict, full_link_dict, args, logger=logger)


    # parse clusters file
    if args.clusters.endswith('.clusters.txt'):
        ctg_group_dict, group_RE_dict = parse_clusters(args.clusters, RE_site_dict, fa_dict, args.min_group_len)
    else:
        ctg_group_dict, group_RE_dict = parse_assembly(args.clusters, RE_site_dict, fa_dict, args.min_group_len)

    # add ungrouped contigs
    grouped_ctgs = add_ungrouped_ctgs(fa_dict, ctg_group_dict)

    # reassignment
    ctg_group_link_dict, linked_ctg_dict = parse_link_dict(full_link_dict, ctg_group_dict, normalize_by_nlinks=args.normalize_by_nlinks)

    preparation_time = time.time()
    logger.info('File parsing and data preparation finished in {}s'.format(preparation_time-start_time))

    if 'whitelist' in args:
        whitelist = args.whitelist
    else:
        whitelist = set()

    for n in range(args.reassign_nrounds):
        run_reassignment(
                sorted_ctg_list, ctg_group_link_dict, ctg_group_dict, full_link_dict, linked_ctg_dict,
                fa_dict, RE_site_dict, args.gfa, group_RE_dict, args.max_ctg_len, args.min_RE_sites, args.min_links,
                args.min_link_density, args.min_density_ratio, args.ambiguous_cutoff, args.min_group_len, whitelist, n+1)
        if n > 0 and last_round == ctg_group_dict:
            logger.info('[result::round{}] Result has converged after {} rounds of reassignment, break'.format(n+1, n))
            break
        last_round = ctg_group_dict.copy()

    if not args.no_additional_rescue:
        # additional rescue round
        run_reassignment(
                sorted_ctg_list, ctg_group_link_dict, ctg_group_dict, full_link_dict, linked_ctg_dict,
                fa_dict, RE_site_dict, args.gfa, group_RE_dict, args.max_ctg_len, args.min_RE_sites, args.min_links,
                args.min_link_density, args.min_density_ratio, args.ambiguous_cutoff, args.min_group_len, whitelist, 0)

    # cluster output
    os.mkdir('reassigned_groups')
    group_ctg_dict, group_hiconf_RE_dict = stat_clusters(ctg_group_dict, fa_dict, grouped_ctgs)
    new_ctg_group_dict, new_group_ctg_dict, new_old_group_dict = clusters_output(group_ctg_dict, fa_dict, 'reassigned')

    reassignment_time = time.time()
    logger.info('{} round(s) of reassignment finished in {}s, average {}s per round'.format(
        n+1, reassignment_time-preparation_time, (reassignment_time-preparation_time)/(n+1)))

    run_hc = False
    # Agglomerative Hierarchical Clustering
    if args.nclusters and args.nclusters < len(new_group_ctg_dict):
        os.mkdir('hc_groups')
        hc_cluster_dict = agglomerative_hierarchical_clustering(
                full_link_dict, grouped_ctgs, new_ctg_group_dict, group_hiconf_RE_dict, new_old_group_dict, args.nclusters, normalize_by_nlinks=args.normalize_by_nlinks)
        hc_group_ctg_dict = stat_hc_clusters(new_group_ctg_dict, hc_cluster_dict)
        new_ctg_group_dict, new_group_ctg_dict, _ = clusters_output(hc_group_ctg_dict, fa_dict, 'hc')
        run_hc = True
    elif args.nclusters == 0:
        logger.info('Parameter --nclusters is set to 0, skip Agglomerative Hierarchical Clustering')
    elif args.nclusters == len(new_group_ctg_dict):
        logger.info('Parameter --nclusters ({}) is equal to the number of clusters ({}) after reassignment, '
                    'skip Agglomerative Hierarchical Clustering'.format(args.nclusters, len(new_group_ctg_dict)))
    else:
        logger.info('Parameter --nclusters ({}) is greater than the number of clusters ({}) after reassignment, '
                    'try higher inflations'.format(args.nclusters, len(new_group_ctg_dict)))

    # split clm file into subfiles
    if run_hc:
        split_clm_file(args.clm, new_group_ctg_dict, new_ctg_group_dict, 'hc_groups')
    else:
        split_clm_file(args.clm, new_group_ctg_dict, new_ctg_group_dict, 'reassigned_groups')

    finish_time = time.time()
    logger.info('Program finished in {}s'.format(finish_time-start_time))


def main():

    # get arguments
    args = parse_arguments()

    run(args, 'HapHiC_reassign.log')


if __name__ == '__main__':
    main()

