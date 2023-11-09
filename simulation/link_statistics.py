#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-09-20 15:04


import sys
import pysam
import argparse


def parse_fasta(fasta):
    
    fa_dict = dict()
    
    with open(fasta) as f:
        
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                # intra_chrom, inter_homo, inter_nonhomo
                fa_dict[ID] = [0, 0, 0]
    
    return fa_dict


def parse_bam(bam, fa_dict, threads):
    
    format_options = [b'filter=flag.read1 && refid != mrefid']
    
    with pysam.AlignmentFile(bam, threads=threads, format_options=format_options) as f:
        
        for aln in f:

            reference_name, next_reference_name = aln.reference_name, aln.next_reference_name
            if reference_name.split('_')[0] == next_reference_name.split('_')[0]:
                if reference_name.split('_')[1] != next_reference_name.split('_')[1]:
                    ctg_pair_type = 'inter_homo'
                    fa_dict[reference_name][1] += 1
                    fa_dict[next_reference_name][1] += 1
                else:
                    ctg_pair_type = 'intra_chrom'
                    fa_dict[reference_name][0] += 1
                    fa_dict[next_reference_name][0] += 1
            else:
                ctg_pair_type = 'inter_nonhomo'
                fa_dict[reference_name][2] += 1
                fa_dict[next_reference_name][2] += 1

            # print('{}\t{}\t{}'.format(reference_name, next_reference_name, ctg_pair_type), file=sys.stderr)


def output_statistics(fa_dict, tag):

    with open('{}_HiC_links.txt'.format(tag), 'w') as f:
        
        for ctg, link_list in fa_dict.items():
            
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                ctg, link_list[0], link_list[1], link_list[2],
                sum(link_list), tag))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('bam', help='input bam file')
    parser.add_argument('tag', help='a specified tag for output file')
    parser.add_argument('--threads', type=int, default=8, help='threads for bam reading, default: %(default)s')
    args = parser.parse_args()

    fa_dict = parse_fasta(args.fasta)

    parse_bam(args.bam, fa_dict, args.threads)

    output_statistics(fa_dict, args.tag)

if __name__ == '__main__':
    main()

