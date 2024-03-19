#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-03-08 10:26


import argparse
import collections
import re
import os

def parse_agp(query_agp):
    
    alignments_dict = collections.defaultdict(list)
    ref_chrom_dict = collections.defaultdict(list)
    with open(query_agp) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            scaffold_name = cols[0]
            if cols[4] == 'W':
                if scaffold_name not in alignments_dict:
                    alignments_dict[scaffold_name] = [0]
                offset = alignments_dict[scaffold_name][0]
                scaffold_start, scaffold_end = int(cols[1]) - offset, int(cols[2]) - offset
                ctg_name, ctg_len, aln_orient = cols[5], int(cols[7]), cols[8]
                ref_chrom, ctg_order, ctg_orient = ctg_name.rsplit('_', 2)
                ctg_order = int(ctg_order[3:])
                if aln_orient == ctg_orient:
                    orient = '+'
                else:
                    orient = '-'
                alignments_dict[scaffold_name].append([scaffold_start, scaffold_end, ctg_name, orient])
                ref_chrom_dict[ref_chrom].append([ctg_name, ctg_len, ctg_order])
            elif cols[4] in {'U', 'N'}:
                alignments_dict[scaffold_name][0] += int(cols[5])

    for ref_chrom, ctg_list in ref_chrom_dict.items():
        ctg_list.sort(key=lambda x: x[-1])
    
    return alignments_dict, ref_chrom_dict


def parse_fasta(fasta):
    
    splitext = os.path.splitext(os.path.basename(fasta))
    fasta_chrs = '{}.chrs.nogaps{}'.format(splitext[0], splitext[1])
    
    seq_dict = dict()
    with open(fasta) as fin:
        for line in fin:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                seq_dict[ID] = []
            else:
                seq_dict[ID] += line.strip()
    
    with open(fasta_chrs, 'w') as fout:
        for ID, seq_list in seq_dict.items():
            seq = ''.join(seq_list).upper().replace('N', '')
            fout.write('>{}\n{}\n'.format(ID, seq))
            seq_dict[ID] = len(seq)

    return os.path.abspath(fasta_chrs), seq_dict


def mock_delta_file(alignments_dict, ref_chrom_dict, ref_seq_len_dict, query_seq_len_dict, ref_fasta_chrs, query_fasta_chrs):

    
    ctg_pos_dict = dict()
    for ref_chrom, ctg_list in ref_chrom_dict.items():
        start = 1
        for n, (ctg_name, ctg_len, _) in enumerate(ctg_list):
            end = start + ctg_len - 1
            ctg_pos_dict[ctg_name] = (ref_chrom, start, end)
            start = end + 1
    
    delta_file = 'mock.delta'

    with open(delta_file, 'w') as f:
        
        f.write('{} {}\nNUCMER\n'.format(ref_fasta_chrs, query_fasta_chrs))

        for scaffold_name, alignments in alignments_dict.items():
            last_scaffold_end = 0
            last_ctg_end = 0
            last_ref_chrom = ''
            last_orient = ''
            for scaffold_start, scaffold_end, ctg_name, orient in alignments[1:]:
                ref_chrom, start, end = ctg_pos_dict[ctg_name]
                if orient == '-':
                    start, end = end, start
                    ctg_match = (last_ctg_end == start + 1)
                else:
                    ctg_match = (last_ctg_end == start - 1)
                # merge
                if last_ref_chrom == '' or (ref_chrom == last_ref_chrom and ctg_match and last_scaffold_end == scaffold_start - 1 and last_orient == orient):
                    # the first contig
                    if last_ref_chrom == '':
                        merge_ctg_start = start
                        merge_scaffold_start = scaffold_start
                        last_ref_chrom = ref_chrom
                        last_orient = orient
                    # always update:
                    last_ctg_end = end
                    last_scaffold_end = scaffold_end
                    continue
                # not merge
                else:
                    f.write('>{} {} {} {}\n'.format(last_ref_chrom, scaffold_name, ref_seq_len_dict[last_ref_chrom], query_seq_len_dict[scaffold_name]))
                    if last_orient == '+':
                        f.write('{} {} {} {} 0 0 0\n0\n'.format(merge_ctg_start, last_ctg_end, merge_scaffold_start, last_scaffold_end))
                    else:
                        f.write('{} {} {} {} 0 0 0\n0\n'.format(last_ctg_end, merge_ctg_start, last_scaffold_end, merge_scaffold_start))
                    merge_ctg_start = start
                    merge_scaffold_start = scaffold_start
                    last_ctg_end = end
                    last_scaffold_end = scaffold_end
                    last_ref_chrom = ref_chrom
                    last_orient = orient

            f.write('>{} {} {} {}\n'.format(last_ref_chrom, scaffold_name, ref_seq_len_dict[last_ref_chrom], query_seq_len_dict[scaffold_name]))
            if last_orient == '+':
                f.write('{} {} {} {} 0 0 0\n0\n'.format(merge_ctg_start, last_ctg_end, merge_scaffold_start, last_scaffold_end))
            else:
                f.write('{} {} {} {} 0 0 0\n0\n'.format(last_ctg_end, merge_ctg_start, last_scaffold_end, merge_scaffold_start))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('query_agp', help='input agp file for query assembly')
    parser.add_argument('query_fasta', help='fasta file for query assembly')
    parser.add_argument('ref_fasta', help='fasta file for reference genome')
    args = parser.parse_args()

    alignments_dict, ref_chrom_dict = parse_agp(args.query_agp)
    ref_fasta_chrs, ref_seq_len_dict = parse_fasta(args.ref_fasta)
    query_fasta_chrs, query_seq_len_dict = parse_fasta(args.query_fasta)
    mock_delta_file(alignments_dict, ref_chrom_dict, ref_seq_len_dict, query_seq_len_dict, ref_fasta_chrs, query_fasta_chrs)

if __name__ == '__main__':
    main()

