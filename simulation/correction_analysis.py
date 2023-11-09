#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-05-08 17:46

import argparse
import pysam
from numpy import array, zeros, int32
import collections
import sys

pysam.set_verbosity(0)

def parse_raw_asm(raw):

    len_dict = dict()

    with open(raw) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                len_dict[ID] = 0
            else:
                len_dict[ID] += len(line.strip())

    return len_dict


def parse_ctg_anno(ctg_anno, len_dict):
    
    ctg_anno_dict = dict()
    chimeric_set, nonchimeric_set = set(), set()
    with open(ctg_anno) as f:
        for line in f:
            line = line.strip()
            if line.endswith(':'):
                ctg_type = line[:-1]
            elif line.startswith('CM'):
                try:
                    assert line in len_dict
                except:
                    print(line)
                ctg_anno_dict[line] = ctg_type
                if ctg_type  == 'nonchimeric':
                    nonchimeric_set.add(line)
                else:
                    chimeric_set.add(line)

    return ctg_anno_dict, chimeric_set, nonchimeric_set


def parse_autohic_agp(autohic_agp, ctg_anno_dict):
    
    autohic_correction_dict = dict()
    autohic_broken_set, autohic_unbroken_set = set(), set()
    with open(autohic_agp) as f:
        for line in f:
            if 'CM039579.1' in line:
                cols = line.split()
                ctg = cols[5]
                # broken
                if 'break' in ctg:
                    # attention: I have an a priori knowledge which all contigs 
                    # in CM039579.1 were broken by AutoHiC at most once, so the
                    # code below is rational
                    ctg, breakn = ctg.rsplit('_', 1)
                    if ctg in ctg_anno_dict:
                        if breakn == 'break1':
                            autohic_correction_dict[ctg] = [int(cols[7])]
                            autohic_broken_set.add(ctg)
                # unbroken
                elif ctg in ctg_anno_dict:
                    autohic_correction_dict[ctg] = 'unbroken'
                    autohic_unbroken_set.add(ctg)

    return autohic_correction_dict, autohic_broken_set, autohic_unbroken_set


def parse_allhic_agp(allhic_agp, ctg_anno_dict, len_dict):

    allhic_correction_dict = dict()
    allhic_broken_set, allhic_unbroken_set = set(), set()
    with open(allhic_agp) as f:
        for line in f:
            if 'CM039579.1' in line:
                cols = line.split()
                ctg = cols[5]
                # broken
                if ctg[-1] not in {'+', '-'}:
                    ctg, start, end = ctg.rsplit('_', 2)
                    if ctg in ctg_anno_dict:
                        end = int(end)
                        if end < len_dict[ctg]:
                            if ctg not in allhic_correction_dict:
                                allhic_correction_dict[ctg] = [end]
                            else:
                                allhic_correction_dict[ctg].append(end)
                        allhic_broken_set.add(ctg)
                elif ctg in ctg_anno_dict:
                    allhic_correction_dict[ctg] = 'unbroken'
                    allhic_unbroken_set.add(ctg)

    return allhic_correction_dict, allhic_broken_set, allhic_unbroken_set


def parse_bam(bam, ctg_anno_dict, len_dict):
    
    resolution = 10000
    ctg_spanning_cov_dict = dict()
    ctg_HiC_link_dict = dict()
    
    for ctg in ctg_anno_dict:
        ctg_len = len_dict[ctg]
        nbins = ctg_len // resolution + 1
        ctg_spanning_cov_dict[ctg] = array([0] * nbins, dtype=int32)
        # hap1, hap2, hap3, hap4, other_chrom
        ctg_HiC_link_dict[ctg] = zeros((5, nbins), dtype=int32) 
    
    format_options = [b'filter=flag.read1']
    with pysam.AlignmentFile(bam, mode='rb', threads=12, format_options=format_options) as f:
        for aln in f:
            ref, mref = aln.reference_name, aln.next_reference_name
            if ref == mref and ref in ctg_anno_dict:
                # spanning coverage
                spanning_region = sorted([aln.reference_start, aln.next_reference_start])
                spanning_start_bin = spanning_region[0]//resolution
                spanning_end_bin = spanning_region[1]//resolution
                ctg_spanning_cov_dict[ref][spanning_start_bin:spanning_end_bin+1] += 1
                # intra-contig links are excluded
                # # Hi-C links
                # bin_n = aln.reference_start // resolution
                # mbin_n = aln.next_reference_start // resolution
                # ctg_HiC_link_dict[ref][0, bin_n] += 1
                # ctg_HiC_link_dict[mref][0, mbin_n] += 1
            else:
                if ref in ctg_anno_dict:
                    bin_n = aln.reference_start // resolution
                    # hap1
                    if mref.startswith('CM039579.1'):
                        ctg_HiC_link_dict[ref][0, bin_n] += 1
                    # hap2
                    elif mref.startswith('CM039580.1'):
                        ctg_HiC_link_dict[ref][1, bin_n] += 1
                        
                    # hap3
                    elif mref.startswith('CM039581.1'):
                        ctg_HiC_link_dict[ref][2, bin_n] += 1
                    # hap4
                    elif mref.startswith('CM039582.1'):
                        ctg_HiC_link_dict[ref][3, bin_n] += 1
                    # other_chrom
                    else:
                        ctg_HiC_link_dict[ref][4, bin_n] += 1
                if mref in ctg_anno_dict:
                    bin_n = aln.next_reference_start // resolution
                    # hap1
                    if ref.startswith('CM039579.1'):
                        ctg_HiC_link_dict[mref][0, bin_n] += 1
                    # hap2
                    elif ref.startswith('CM039580.1'):
                        ctg_HiC_link_dict[mref][1, bin_n] += 1
                    # hap3
                    elif ref.startswith('CM039581.1'):
                        ctg_HiC_link_dict[mref][2, bin_n] += 1
                    # hap4
                    elif ref.startswith('CM039582.1'):
                        ctg_HiC_link_dict[mref][3, bin_n] += 1
                    # other_chrom
                    else:
                        ctg_HiC_link_dict[mref][4, bin_n] += 1

    return ctg_spanning_cov_dict, ctg_HiC_link_dict


def output_statistics(
            ctg_anno_dict, ctg_spanning_cov_dict, ctg_HiC_link_dict, autohic_correction_dict, allhic_correction_dict, 
            chimeric_set, nonchimeric_set, autohic_broken_set, autohic_unbroken_set, allhic_broken_set, allhic_unbroken_set):
    
    resolution = 10000

    # statistics for plotting
    max_cov_dict = collections.defaultdict(int)
    for ctg, spanning_cov_list in ctg_spanning_cov_dict.items():
        for n, cov in enumerate(spanning_cov_list):
            if cov > max_cov_dict[ctg]:
                max_cov_dict[ctg] = cov
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, n*resolution+1, (n+1)*resolution, 0, cov, 'Spanning_coverage', 'Spanning_coverage'))

    for ctg, link_array in ctg_HiC_link_dict.items():
        for n in range(link_array.shape[1]):
            coefficient = sum(link_array[:,n]) / max_cov_dict[ctg]
            if coefficient:
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, n*resolution+1, (n+1)*resolution, 0, link_array[0,n]/coefficient, 'nHiC_links', 'Hap1'))
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, n*resolution+1, (n+1)*resolution, link_array[0,n]/coefficient, (link_array[0,n]+link_array[1,n])/coefficient, 'nHiC_links', 'Hap2'))
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, n*resolution+1, (n+1)*resolution, (link_array[0,n]+link_array[1,n])/coefficient, (link_array[0,n]+link_array[1,n]+link_array[2,n])/coefficient, 'nHiC_links', 'Hap3'))
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, n*resolution+1, (n+1)*resolution, (link_array[0,n]+link_array[1,n]+link_array[2,n])/coefficient, (link_array[0,n]+link_array[1,n]+link_array[2,n]+link_array[3,n])/coefficient, 'nHiC_links', 'Hap4'))
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(ctg, n*resolution+1, (n+1)*resolution, (link_array[0,n]+link_array[1,n]+link_array[2,n]+link_array[3,n])/coefficient, max_cov_dict[ctg], 'nHiC_links', 'Other_chrom'))
    for ctg, break_points in autohic_correction_dict.items():
        if break_points == 'unbroken':
            continue
        for break_point in break_points:
            print('{}\t{}\tNA\t{}\tNA\tBreakpoint\tAutoHiC'.format(ctg, break_point, max_cov_dict[ctg]))
    
    for ctg, break_points in allhic_correction_dict.items():
        if break_points == 'unbroken':
            continue
        for break_point in break_points:
            print('{}\t{}\tNA\t{}\tNA\tBreakpoint\tALLHiC'.format(ctg, break_point, max_cov_dict[ctg]))

    # for venn diagram
    with open('venn.txt', 'w') as f:
        f.write('Contig\tType\tAutoHiC\tALLHiC\n')
        for ctg, ctg_type in ctg_anno_dict.items():
            if ctg in autohic_broken_set:
                autohic_broken = 'yes'
            else:
                autohic_broken = 'no'
            if ctg in allhic_broken_set:
                allhic_broken = 'yes'
            else:
                allhic_broken = 'no'
            f.write('{}\t{}\t{}\t{}\n'.format(ctg, ctg_type, autohic_broken, allhic_broken))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw', help='raw assembly in FASTA format')
    parser.add_argument('bam', help='Hi-C mapping in BAM format (on raw assembly)')
    parser.add_argument('ctg_anno', help='contig annotation list')
    parser.add_argument('autohic_agp', help='agp file generated by AutoHiC (after correction)')
    parser.add_argument('allhic_agp', help='agp file generated by ALLHiC (after correction)')
    args = parser.parse_args()
    
    # parse ctg annotation list
    len_dict = parse_raw_asm(args.raw)
    ctg_anno_dict, chimeric_set, nonchimeric_set = parse_ctg_anno(args.ctg_anno, len_dict)
    
    # parse autohic agp
    autohic_correction_dict, autohic_broken_set, autohic_unbroken_set = parse_autohic_agp(args.autohic_agp, ctg_anno_dict)
    
    # parse allhic agp
    allhic_correction_dict, allhic_broken_set, allhic_unbroken_set = parse_allhic_agp(args.allhic_agp, ctg_anno_dict, len_dict)
    
    # read bam file
    ctg_spanning_cov_dict, ctg_HiC_link_dict = parse_bam(args.bam, ctg_anno_dict, len_dict)
    
    # statistics
    output_statistics(
            ctg_anno_dict, ctg_spanning_cov_dict, ctg_HiC_link_dict, autohic_correction_dict, allhic_correction_dict, 
            chimeric_set, nonchimeric_set, autohic_broken_set, autohic_unbroken_set, allhic_broken_set, allhic_unbroken_set)

if __name__ == '__main__':
    main()
