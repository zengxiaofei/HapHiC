#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-12-01 12:17


import argparse
import pickle
import logging
import time
import sys
import gzip

import pysam
from portion import closed
import collections
import numpy as np
from math import ceil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

pysam.set_verbosity(0)


def parse_agp(agp, bin_size):

    logger.info('Parsing input AGP file...')

    ctg_dict = collections.defaultdict(dict)
    ctg_aln_dict = collections.defaultdict(dict)
    group_size_dict = collections.OrderedDict()

    frag_set = set()
    group_frag_dict = collections.defaultdict(set)

    with open(agp) as f:

        for line in f:

            if line.startswith('#') or not line.strip():
                continue

            cols = line.split()

            if cols[4] == 'W':
                group = cols[0]
                group_start, group_end = int(cols[1]), int(cols[2])
                group_start_bin = (group_start - 1) // bin_size
                group_end_bin = (group_end - 1) // bin_size

                ctg = cols[5]
                ctg_raw_start, ctg_raw_end = int(cols[6]), int(cols[7])

                group_size_dict[group] = group_end

                frag_id = (ctg, ctg_raw_start, ctg_raw_end)
                frag_set.add(frag_id)
                group_frag_dict[group].add(frag_id)

                for group_bin in range(group_start_bin, group_end_bin + 1):

                    group_bin_start = group_bin * bin_size + 1
                    group_bin_end = (group_bin + 1) * bin_size
                    ctg_group_range = closed(group_bin_start, group_bin_end) & closed(group_start, group_end)

                    # calculate bin_start and bin_end on raw contigs
                    if cols[8] == '+':
                        ctg_raw_bin_start = ctg_group_range.lower - group_start + ctg_raw_start
                        ctg_raw_bin_end = ctg_group_range.upper - group_start + ctg_raw_start
                    else:
                        assert cols[8] == '-'
                        ctg_raw_bin_start = ctg_raw_end - (ctg_group_range.upper - group_start)
                        ctg_raw_bin_end =  ctg_raw_end - (ctg_group_range.lower - group_start)

                    ctg_raw_bin_range = closed(ctg_raw_bin_start, ctg_raw_bin_end)
                    ctg_dict[ctg][ctg_raw_bin_range] = (group, group_bin)

                    ctg_aln_start_bin = (ctg_raw_bin_start-1)//bin_size
                    ctg_aln_end_bin = (ctg_raw_bin_end-1)//bin_size

                    for aln_bin in range(ctg_aln_start_bin, ctg_aln_end_bin+1):
                        if aln_bin in ctg_aln_dict[ctg]:
                            ctg_aln_dict[ctg][aln_bin].append(ctg_raw_bin_range)
                        else:
                            ctg_aln_dict[ctg][aln_bin] = [ctg_raw_bin_range]

    return ctg_dict, ctg_aln_dict, group_size_dict, frag_set, group_frag_dict


def generate_contact_matrix(group_size_dict, frag_set, group_frag_dict, bin_size, min_len, specified_scaffolds):


    logger.info('Generating an empty contact matrix...')

    total_n_bins = 0
    min_len *= 1000000
    
    if specified_scaffolds:
        scaffolds = specified_scaffolds.split(',')

    # group bin ID to total bin ID
    group_to_total_bin_dict = dict()
    group_list = list()

    for group, size in group_size_dict.items():
        if size >= min_len and ((specified_scaffolds and group in scaffolds) or not specified_scaffolds):
            group_n_bins = size // bin_size + 1
            for n in range(group_n_bins):
                group_to_total_bin_dict[(group, n)] = total_n_bins + n
            total_n_bins += group_n_bins
            group_list.append(group)
        else:
            frag_set -= group_frag_dict[group]

    contact_matrix = np.zeros((total_n_bins, total_n_bins), dtype=int)

    ctg_set = set()
    for ctg, _, __ in frag_set:
        ctg_set.add(ctg)

    return contact_matrix, group_to_total_bin_dict, group_list, ctg_set


def parse_pairs(pairs, ctg_dict, ctg_aln_dict, bin_size, contact_matrix, group_to_total_bin_dict, group_list, ctg_set):

    def convert_group_bin_id(ctg, pos):

        for ctg_raw_bin_range in ctg_aln_dict[ctg][(pos-1)//bin_size]:
            group_and_bin = ctg_dict[ctg][ctg_raw_bin_range]
            if pos in ctg_raw_bin_range:
                if group_and_bin[0] not in group_list:
                    return None
                return group_to_total_bin_dict[group_and_bin]

    logger.info('Parsing input pairs file...')

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
            ref, pos, mref, mpos = cols[1], int(cols[2]), cols[3], int(cols[4])

            # only consider contigs in groups with a length equal or greater than `min_len`
            if ref not in ctg_set or mref not in ctg_set:
                continue

            ref_group_bin_id = convert_group_bin_id(ref, pos)
            if ref_group_bin_id is None:
                continue

            mref_group_bin_id = convert_group_bin_id(mref, mpos)
            if mref_group_bin_id is None:
                continue

            contact_matrix[ref_group_bin_id, mref_group_bin_id] += 1

    return contact_matrix


def parse_bam(bam, ctg_dict, ctg_aln_dict, bin_size, contact_matrix, group_to_total_bin_dict, group_list, ctg_set, threads):

    def convert_group_bin_id(ctg, pos):

        for ctg_raw_bin_range in ctg_aln_dict[ctg][(pos-1)//bin_size]:
            group_and_bin = ctg_dict[ctg][ctg_raw_bin_range]
            if pos in ctg_raw_bin_range:
                if group_and_bin[0] not in group_list:
                    return None
                return group_to_total_bin_dict[group_and_bin]


    logger.info('Parsing input BAM file...')

    format_options = [b'filter=flag.read1']

    with pysam.AlignmentFile(bam, mode='rb', threads=threads, format_options=format_options) as f:

        for aln in f:

            ref, mref = aln.reference_name, aln.next_reference_name

            # only consider contigs in groups with a length equal or greater than `min_len`
            if ref not in ctg_set or mref not in ctg_set:
                continue

            ref_group_bin_id = convert_group_bin_id(ref, aln.reference_start + 1)
            if ref_group_bin_id is None:
                continue

            mref_group_bin_id = convert_group_bin_id(mref, aln.next_reference_start + 1)
            if mref_group_bin_id is None:
                continue

            contact_matrix[ref_group_bin_id, mref_group_bin_id] += 1

    return contact_matrix


def output_pickle(contact_matrix, args):

    logger.info('Writing raw contact matrix to a pickle file...')

    with open('contact_matrix.pkl', 'wb') as fpkl:
        pickle.dump((contact_matrix, args), fpkl)


def load_pickle(pickle_file, args):

    logger.info('Reading raw contact matrix from a previously generated pickle file...')

    with open(pickle_file, 'rb') as fpkl:
        contact_matrix, old_args = pickle.load(fpkl)
        # check parameters
        if old_args.bin_size != args.bin_size or old_args.min_len != args.min_len or old_args.specified_scaffolds != args.specified_scaffolds:
            error_message = ('The input parameters (--bin_size {} --min_len {} --specified_scaffolds {}) are not consistent with '
                             'those used to generate `contact_map.pkl` (--bin_size {} --min_len {} --specified_scaffolds {})'.format(
                                args.bin_size, args.min_len, args.specified_scaffolds, 
                                old_args.bin_size, old_args.min_len, old_args.specified_scaffolds))
            logger.error(error_message)
            raise Exception(error_message)

        return contact_matrix


def bnewt(A, tol=1e-6, x0=None, delta=0.1, Delta=3, fl=0):

    # A python implemention of the Knight-Ruiz (KR) normalization algorithm described in:
    # https://academic.oup.com/imajna/article/33/3/1029/659457

    n = A.shape[0]
    e = np.ones(n)
    res = []

    if x0 is None:
        x0 = e

    g = 0.9
    etamax = 0.1
    eta = etamax
    stop_tol = tol * 0.5

    x = x0
    rt = tol**2
    v = x * (A @ x)
    rk = 1 - v
    rho_km1 = rk @ rk
    rout = rho_km1
    rold = rout

    MVP = 0
    i = 0

    if fl == 1:
        # print('it in. it res')
        pass

    nn = 0
    max_nn, max_mm = 1000, 10000
    error_message = (
            'Unable to converge. Maybe the matrix is too sparse (too few Hi-C links). '
            'You can try another normalization method.')

    while rout > rt:

        # to avoid endless loop
        nn += 1
        mm = 0
        if nn > max_nn:
            logger.info(error_message)
            raise Exception(error_message)

        i += 1
        k = 0
        y = e
        innertol = max([eta**2 * rout, rt])

        while rho_km1 > innertol:

            # to avoid endless loop
            mm += 1
            if mm > max_mm:
                logger.info(error_message)
                raise Exception(error_message)

            k += 1
            if k == 1:
                Z = rk / v
                p = Z
                rho_km1 = rk @ Z
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            w = x * (A @ (x * p)) + v * p
            alpha = rho_km1 / (p @ w)
            ap = alpha * p

            ynew = y + ap
            if min(ynew) <= delta:
                if delta == 0:
                    break
                ind = np.where(ap < 0)
                gamma = min((delta - y[ind]) / ap[ind])
                y = y + gamma * ap
                break
            if max(ynew) >= Delta:
                ind = np.where(ynew > Delta)
                gamma = min((Delta - y[ind]) / ap[ind])
                y = y + gamma * ap
                break
            y = ynew
            rk = rk - alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = rk @ Z

        x = x * y
        v = x * (A @ x)
        rk = 1 - v
        rho_km1 = rk @ rk
        rout = rho_km1
        MVP += k + 1

        rat = rout / rold
        rold = rout
        res_norm = np.sqrt(rout)
        eta_o = eta
        eta = g * rat
        if g * eta_o**2 > 0.1:
            eta = max([eta, g * eta_o**2])
        eta = max([min([eta, etamax]), stop_tol / res_norm])

        if fl == 1:
            # print(f'{i:3d} {k:6d} {res_norm:.3e}')
            res.append(res_norm)

    # print(f'Matrix-vector products = {MVP:6d}')
    return x, res


def normalize_matrix(contact_matrix, group_list, group_size_dict, bin_size, normalization, vmax_coef, manual_vmax):

    if normalization == 'KR':

        # a dict used to save intra-scaffold matrix for each scaffold
        normalized_intra_matrix_dict = dict()

        logger.info('Normalizing contact mattrix using the Knight-Ruiz (KR) balancing algorithm')
        group_start_bin = 0

        # save indices for zero-value elements
        zero_indices = np.argwhere(contact_matrix == 0)
        # to avoid divide-by-zero warning/error
        contact_matrix = contact_matrix + 0.00001

        for group in group_list:
            group_bin_num = ceil(group_size_dict[group]/bin_size)
            group_end_bin = group_start_bin + group_bin_num
            intra_matrix = contact_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin]
            group_start_bin += group_bin_num

            # intra-scaffold KR normalization for each scaffold
            x, _ = bnewt(intra_matrix)
            d = np.diag(x)
            normalized_intra_matrix = d @ intra_matrix @ d
            normalized_intra_matrix_dict[group] = normalized_intra_matrix

        # inter-scaffold KR normalization
        x, _ = bnewt(contact_matrix)
        d = np.diag(x)
        normalized_inter_matrix = d @ contact_matrix @ d

        # combine inter- and intra- matrices, and calculate vmax
        non_diagonal_list = list()
        group_start_bin = 0

        for group in group_list:
            group_bin_num = ceil(group_size_dict[group]/bin_size)
            group_end_bin = group_start_bin + group_bin_num
            normalized_inter_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin] = normalized_intra_matrix_dict[group]
            for n, row in enumerate(normalized_intra_matrix_dict[group]):
                for m, i in enumerate(row):
                    if n != m:
                        non_diagonal_list.append(i)
            group_start_bin += group_bin_num

        # retrieve zero values
        normalized_inter_matrix[zero_indices[:, 0], zero_indices[:, 1]] = 0

        if manual_vmax < 0:
            vmax = np.median(non_diagonal_list) * vmax_coef
            logger.info('The vmax for the KR-normalized matrix is calculated to be {} ({} * median)'.format(vmax, vmax_coef))
        else:
            vmax = manual_vmax
            logger.info('The vmax for the KR-normalized matrix is manually designated as {})'.format(vmax))

        return normalized_inter_matrix, vmax

    else:

        if normalization == 'log10':
            logger.info('Normalizing contact matrix using log10...')
            normalized_matrix = np.log10(contact_matrix + 1)
        else:
            logger.info('Normalization is disabled')
            normalized_matrix = contact_matrix

        # calculate vmax
        non_diagonal_list = list()
        group_start_bin = 0

        for group in group_list:
            group_bin_num = ceil(group_size_dict[group]/bin_size)
            group_end_bin = group_start_bin + group_bin_num
            intra_matrix = normalized_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin]
            for n, row in enumerate(intra_matrix):
                for m, i in enumerate(row):
                    if n != m:
                        non_diagonal_list.append(i)
            group_start_bin += group_bin_num

        if manual_vmax < 0:
            vmax = np.median(non_diagonal_list) * vmax_coef
        else:
            vmax = manual_vmax

        if normalization == 'log10':
            if manual_vmax < 0:
                logger.info('The vmax for the log-normalized matrix is calculated to be {} ({} * median)'.format(vmax, vmax_coef))
            else:
                logger.info('The vmax for the log-normalized matrix is manually designated as {}'.format(vmax))
        else:
            if manual_vmax < 0:
                logger.info('The vmax for the raw matrix is calculated to be {} ({} * median)'.format(vmax, vmax_coef))
            else:
                logger.info('The vmax for the raw matrix is manually designated as {}'.format(vmax))

        return normalized_matrix, vmax


def get_resolution(length, n):

    for resolution in (0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000):
        n_xticks = length // int(resolution * 1000000) + 1
        if n_xticks < n:
            return resolution


def get_xticks(length, resolution, bin_size):

    xtick_list, xtick_values = list(), list()
    for x in range(0, length+1, int(resolution * 1000000)):
        xtick_list.append(x // bin_size)
        xtick_values.append(x / 1000000)

    if any([True for v in xtick_values if v != int(v)]):
        xticklabel_list = [str(v) for v in xtick_values]
    else:
        xticklabel_list = [str(int(v)) for v in xtick_values]

    return xtick_list, xticklabel_list


def set_origin(ax, n_bins, origin):

    if origin == 'bottom_left':
        ax.axis([0, n_bins, 0, n_bins])
    elif origin == 'top_left':
        ax.axis([0, n_bins, n_bins, 0])
    elif origin == 'bottom_right':
        ax.axis([n_bins, 0, 0, n_bins])
    else:
        ax.axis([n_bins, 0, n_bins, 0])


def get_cmap(cmap):

    if ',' in cmap:
        color_list = cmap.split(',')
        return colors.LinearSegmentedColormap.from_list('my_cmap', color_list)
    else:
        return cmap


def normalized_imshow(ax, contact_matrix, cmap, vmax):

    heatmap = ax.imshow(contact_matrix, cmap=cmap)
    norm = colors.Normalize(vmin=0, vmax=vmax)
    heatmap.set_norm(norm)

    return heatmap

def get_line_style(line_style):

        if line_style == 'dashed':
            return (0, (10, 20))
        else:
            return 'solid'


def draw_heatmap(contact_matrix, group_list, group_size_dict, bin_size, vmax, args):

    logger.info('Drawing heatmap...')

    plt.rcParams['pdf.fonttype'] = 42
    # convert inch to cm
    fig, ax = plt.subplots(figsize=(args.figure_width/2.54, args.figure_height/2.54), dpi=1000)
    plt.subplots_adjust(bottom=0.05, left=0.15, right=1, top=0.95)

    # set yticks
    ytick = 0
    ytick_list = list()
    group_edge_list = list()
    for group in group_list:
        group_bin_num = ceil(group_size_dict[group]/bin_size)
        ytick_list.append(ytick + group_bin_num / 2)
        ytick += group_bin_num
        group_edge_list.append(ytick - 0.5)
    ax.set_yticks(ytick_list)
    ax.set_yticklabels(group_list, size=6)

    # set xticks
    total_n_bins = contact_matrix.shape[0]
    total_len = total_n_bins * bin_size

    resolution = get_resolution(total_len, 10)

    xtick_list, xticklabel_list = get_xticks(total_len, resolution, bin_size)

    ax.set_xticks(xtick_list)
    ax.set_xticklabels(xticklabel_list, size=6)

    # set the origin for the plot
    set_origin(ax, total_n_bins, args.origin)

    # get cmap for heatmap
    cmap = get_cmap(args.cmap)

    # draw heatmap with normalization
    heatmap = normalized_imshow(ax, contact_matrix, cmap, vmax)

    ax.set_title('Hi-C contact map (bin size: {} Kb, x-axis unit: Mb)'.format(args.bin_size), fontsize=7)

    if args.border_style == 'grid':
        # get gridline style
        gridline_style = get_line_style(args.gridline_style)
        # add gridlines
        for edge in group_edge_list[:-1]:
            ax.vlines(edge, 0, total_n_bins, color=args.gridline_color, linestyle=gridline_style, linewidth=args.gridline_width)
            ax.hlines(edge, 0, total_n_bins, color=args.gridline_color, linestyle=gridline_style, linewidth=args.gridline_width)
    else:
        # get outline style
        outline_style = get_line_style(args.outline_style)
        # add outlines
        last_edge = 0
        for edge in group_edge_list:
            ax.vlines(edge, last_edge, edge, color=args.outline_color, linestyle=outline_style, linewidth=args.outline_width)
            ax.hlines(edge, last_edge, edge, color=args.outline_color, linestyle=outline_style, linewidth=args.outline_width)
            ax.vlines(last_edge, last_edge, edge, color=args.outline_color, linestyle=outline_style, linewidth=args.outline_width)
            ax.hlines(last_edge, last_edge, edge, color=args.outline_color, linestyle=outline_style, linewidth=args.outline_width)
            last_edge = edge

    cb = fig.colorbar(heatmap, shrink=0.5)
    if args.normalization== 'KR':
        cb.set_label('KR normalized counts', fontsize=7)
    elif args.normalization == 'log10':
        cb.set_label('Log$_{10}$(counts+1)', fontsize=7)
    else:
        cb.set_label('counts', fontsize=7)
    cb.ax.tick_params(labelsize=6)

    plt.savefig('contact_map.pdf')
    plt.close()


def draw_separate_heatmaps(contact_matrix, group_list, group_size_dict, bin_size, vmax, args):

    logger.info('Drawing heatmap for each scaffold...')

    plt.rcParams['pdf.fonttype'] = 42
    group_num = len(group_list)

    # ratio = height / width
    nrows = ceil(group_num/args.ncols)
    ratio = nrows / args.ncols

    # convert inch to cm
    fig, axes = plt.subplots(
            ceil(len(group_list)/args.ncols), args.ncols,
            figsize=(args.figure_width/2.54, ((args.figure_width-2)*ratio*1.2)/2.54), dpi=1000)

    ax_list = list()
    for n, ax in enumerate(axes.flat):
        if n < group_num:
            ax_list.append(ax)
        elif group_num > args.ncols:
            fig.delaxes(axes[n//args.ncols, n%args.ncols])
        else:
            fig.delaxes(axes[n%args.ncols])

    group_start_bin = 0
    for n, group in enumerate(group_list):
        group_bin_num = ceil(group_size_dict[group]/bin_size)
        group_end_bin = group_start_bin + group_bin_num
        group_matrix = contact_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin]
        group_start_bin += group_bin_num
        ax = ax_list[n]

        # set x and yticks
        group_n_bins = group_matrix.shape[0]
        group_len = group_n_bins * bin_size

        resolution = get_resolution(group_len, 6)

        xtick_list, xticklabel_list = get_xticks(group_len, resolution, bin_size)

        ax.set_xticks(xtick_list)
        ax.set_xticklabels(xticklabel_list, size=5)
        ax.set_yticks(xtick_list)
        ax.set_yticklabels(xticklabel_list, size=5)

        # set the origin for each plot
        set_origin(ax, group_n_bins, args.origin)

        # get cmap for heatmap
        cmap = get_cmap(args.cmap)

        # draw heatmap with normalization
        normalized_imshow(ax, group_matrix, cmap, vmax)

        ax.set_title(group, fontsize=6)

    fig.tight_layout()

    plt.savefig('separate_plots.pdf')
    plt.close()


def parse_arguments():

    parser = argparse.ArgumentParser(prog='haphic plot')

    parser.add_argument(
            'agp', help='scaffolding result in AGP format. The IDs in this file should match those in the BAM file')
    parser.add_argument(
            'alignments', help='filtered Hi-C read alignments in BAM/pairs format (slow) or previously generated `contact_matrix.pkl` (much faster)')
    parser.add_argument(
            '--bin_size', type=int, default=500, 
            help='bin size for generating contact matrix, default: %(default)s (kbp)')
    parser.add_argument(
            '--specified_scaffolds', default=None,
            help='specify scaffolds to visualize, separated with commas, default: %(default)s. When this parameter is set, the `--min_len` parameter will be disabled')
    parser.add_argument(
            '--min_len', type=int, default=1, 
            help='minimum scaffold length for visualization, default: %(default)s (Mbp)')
    parser.add_argument(
            '--cmap', default='white,red', 
            help='define the colormap for the heatmap, default: %(default)s. It can be any built-in sequential colormap from Matplotlib '
            '(refer to: https://matplotlib.org/stable/users/explain/colors/colormaps.html). You can also create a custom colormap by listing '
            'colors separated by commas')
    parser.add_argument(
            '--normalization', choices=('KR', 'log10', 'none'), default='KR', 
            help='method for matrix normalization, default: %(default)s')
    parser.add_argument(
            '--vmax_coef', type=float, default=4.0, 
            help='for contact matrix visualization, values greater than `vmax_coef` times the median of the intra-scaffold matrices '
            'are displayed in the same color, default: %(default)s')
    parser.add_argument(
            '--manual_vmax', type=float, default=-1,
            help='manually designate the maximum value of the colorbar (vmax) using a positive float number, default: disabled')
    parser.add_argument(
            '--separate_plots', default=False, action='store_true', 
            help='generate `separate_plots.pdf`, depicting the heatmap for each scaffold individually, default: %(default)s')
    parser.add_argument(
            '--ncols', type=int, default=5, 
            help='number of scaffolds per row in `separate_plots.pdf`, default: %(default)s')
    parser.add_argument(
            '--origin', choices=('bottom_left', 'top_left', 'bottom_right', 'top_right'), default='bottom_left', 
            help='set the origin of each heatmap, default: %(default)s')
    parser.add_argument(
            '--border_style', choices=('grid', 'outline'), default='grid', 
            help='border style for scaffolds, default: %(default)s')
    parser.add_argument(
            '--gridline_color', default='grey', 
            help='color for gridlines, default: %(default)s')
    parser.add_argument(
            '--gridline_style', choices=('solid', 'dashed'), default='solid', 
            help='style for gridlines, default: %(default)s')
    parser.add_argument(
            '--gridline_width', type=float, default=0.2, 
            help='width for gridlines, default: %(default)s')
    parser.add_argument(
            '--outline_color', default='blue', 
            help='color for outlines, default: %(default)s')
    parser.add_argument(
            '--outline_style', choices=('solid', 'dashed'), default='solid', 
            help='style for outlines, default: %(default)s')
    parser.add_argument(
            '--outline_width', type=float, default=0.2, 
            help='width for outlines, default: %(default)s')
    parser.add_argument(
            '--figure_width', type=int, default=15, 
            help='figure width, default: %(default)s (cm)')
    parser.add_argument(
            '--figure_height', type=int, default=12, 
            help='figure height, default: %(default)s (cm)')
    parser.add_argument(
            '--threads', type=int, default=8, 
            help='number of threads for reading BAM file, default: %(default)s')

    args = parser.parse_args()
    
    # When this parameter is set, the `--min_len` parameter will be disabled
    if args.specified_scaffolds:
        args.min_len = 0

    return args


def main():

    args = parse_arguments()

    log_file = 'HapHiC_plot.log'
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

    bin_size = args.bin_size * 1000

    ctg_dict, ctg_aln_dict, group_size_dict, frag_set, group_frag_dict = parse_agp(args.agp, bin_size)

    contact_matrix, group_to_total_bin_dict, group_list, ctg_set = generate_contact_matrix(
            group_size_dict, frag_set, group_frag_dict, bin_size, args.min_len, args.specified_scaffolds)

    if args.alignments.endswith('.pkl'):
        contact_matrix = load_pickle(args.alignments, args)
    else:
        if args.alignments.endswith('.bam'):
            contact_matrix = parse_bam(
                    args.alignments, ctg_dict, ctg_aln_dict, bin_size, contact_matrix, group_to_total_bin_dict, group_list, ctg_set, args.threads)
        else:
            assert args.alignments.endswith('.pairs') or args.alignments.endswith('.pairs.gz')
            contact_matrix = parse_pairs(
                    args.alignments, ctg_dict, ctg_aln_dict, bin_size, contact_matrix, group_to_total_bin_dict, group_list, ctg_set)
        contact_matrix = contact_matrix + np.transpose(contact_matrix)
        for i in range(contact_matrix.shape[0]):
            contact_matrix[i,i] /= 2
        output_pickle(contact_matrix, args)

    normalized_contact_matrix, vmax = normalize_matrix(contact_matrix, group_list, group_size_dict, bin_size, args.normalization, args.vmax_coef, args.manual_vmax)

    draw_heatmap(normalized_contact_matrix, group_list, group_size_dict, bin_size, vmax, args)

    if args.separate_plots:
        draw_separate_heatmaps(normalized_contact_matrix, group_list, group_size_dict, bin_size, vmax, args)

    finished_time = time.time()
    logger.info('Program finished in {}s'.format(finished_time-start_time))


if __name__ == '__main__':
    main()

