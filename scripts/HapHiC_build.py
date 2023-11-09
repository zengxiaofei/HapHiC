#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-07-11 15:39


import sys
import argparse
import os
import logging
import time

from collections import OrderedDict

from HapHiC_cluster import parse_fasta
from _version import __version__, __update_time__

logging.basicConfig(
        format='%(asctime)s <%(filename)s> [%(funcName)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
        )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_tours(tour_files, fa_dict):

    logger.info('Parsing tour files...')

    output_ctgs = set()
    tour_dict = OrderedDict()

    for tour_file in tour_files:
        basename = os.path.basename(tour_file)
        group = os.path.splitext(basename)[0].split('_')[0]
        tour_dict[group] = list()
        last_line  = ''
        with open(tour_file) as f:
            for line in f:
                if line.strip():
                    last_line = line.strip()

        for ctg_ori in last_line.split():
            ctg = ctg_ori[:-1]
            ori = ctg_ori[-1]
            if ctg not in fa_dict:
                raise RuntimeError('CANNOT find ctg {} in FASTA file'.format(ctg))
            elif ctg in output_ctgs:
                raise RuntimeError('Contig {} is repeated'.format(ctg))
            else:
                output_ctgs.add(ctg)
            tour_dict[group].append((ctg, ori))

    return tour_dict, output_ctgs


def build_final_scaffolds(tour_dict, fa_dict, output_ctgs, args):

    def get_output_seq(group):

        seq_list = list()
        for ctg, ori in tour_dict[group]:
            if ori == '+':
                seq_list.append(fa_dict[ctg][0])
            else:
                seq_list.append(fa_dict[ctg][0].translate(com_tab)[::-1])
        output_seq = ('N'*args.Ns).join(seq_list)

        return output_seq


    def write_agp(group, agp_out):

        n = 0
        accumulated_len = 0
        for ctg, ori in tour_dict[group]:
            n += 1
            ctg_len = fa_dict[ctg][1]
            start = accumulated_len + 1
            end = accumulated_len + ctg_len
            accumulated_len += ctg_len
            agp_out.write('{}\t{}\t{}\t{}\tW\t{}\t1\t{}\t{}\n'.format(
                group, start, end, n, ctg, ctg_len, ori))
            if n < len(tour_dict[group]*2)-1:
                n += 1
                start = accumulated_len + 1
                end = accumulated_len + args.Ns
                accumulated_len += args.Ns
                agp_out.write('{}\t{}\t{}\t{}\tU\t{}\tcontig\tyes\tmap\n'.format(
                    group, start, end, n, args.Ns))

    logger.info('Building final scaffolds...')

    # get complementary base
    base_plus = 'ATCGNatcgn'
    base_minus = 'TAGCNtagcn'
    com_tab = str.maketrans(base_plus, base_minus)

    # contigs in tour files
    # get the output order of final scaffolds
    order_list = list()
    # sort by length
    if not args.sort_by_input:
        for group, ctgs in tour_dict.items():
            length = sum([fa_dict[ctg][1] for ctg, ori in ctgs]) + (len(ctgs)-1)*args.Ns
            order_list.append((group, length))
        order_list.sort(key=lambda x: x[1], reverse=True)
        order_list = [group for group, length in order_list]
    # sort by input order of tour files
    else:
        order_list = tour_dict.keys()

    # contigs that are not in tour files
    unanchored_ctg_list = list()
    for ctg in fa_dict:
        if ctg not in output_ctgs:
            unanchored_ctg_list.append((ctg, fa_dict[ctg][1]))
    unanchored_ctg_list.sort(key=lambda x: x[1], reverse=True)

    # output FASTA file & AGP file
    with open('{}.fa'.format(args.prefix), 'w') as fa_out, open('{}.agp'.format(args.prefix), 'w') as agp_out:
        # the contigs in tour files
        for group in order_list:
            # output FASTA file
            output_seq = get_output_seq(group)
            fa_out.write('>{}\n{}\n'.format(group, output_seq))
            # output AGP file
            write_agp(group, agp_out)
        # contigs that are not in tour files
        for ctg, ctg_len in unanchored_ctg_list:
            # output FASTA file
            fa_out.write('>{}\n{}\n'.format(ctg, fa_dict[ctg][0]))
            # output AGP file
            agp_out.write('{0}\t1\t{1}\t1\tW\t{0}\t1\t{1}\t+\n'.format(
                ctg, ctg_len))


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(
            'fasta', help='draft genome in FASTA format')
    parser.add_argument(
            'tours', nargs='+', help='`*.tour` files generated in the sorting (ordering and orientation) step')
    parser.add_argument(
            '--Ns', type=int, default=100,
            help='number of Ns representing gaps, default: %(default)s')
    parser.add_argument(
            '--sort_by_input',
            default=False, action='store_true',
            help='sort output scaffolds by the order of input `*.tour files`, otherwise by group length, default: %(default)s')
    parser.add_argument(
            '--keep_letter_case',
            default=False, action='store_true',
            help='keep the letter case of bases in the draft genome (usually lowercase represents softmaksing), default: %(default)s')
    parser.add_argument(
            '--prefix', default='scaffolds',
            help='prefix of output FASTA and AGP files, default: %(default)s')

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

    # Using default RE here is ok. Because in the building step,
    # we don't care about the restriction sites.
    fa_dict = parse_fasta(args.fasta, keep_letter_case=args.keep_letter_case, logger=logger)

    tour_dict, output_ctgs = parse_tours(args.tours, fa_dict)

    build_final_scaffolds(tour_dict, fa_dict, output_ctgs, args)

    finish_time = time.time()
    logger.info('Program finished in {}s'.format(finish_time-start_time))


def main():

    # get arguments
    args = parse_arguments()

    run(args, 'HapHiC_build.log')


if __name__ == '__main__':
    main()

