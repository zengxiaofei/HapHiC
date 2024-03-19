#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-03-08 15:34


import argparse

comp_table = str.maketrans('ATCG', 'TAGC')

def parse_fasta(fasta):

    fa_dict = dict()

    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                fa_dict[ID] = []
            else:
                fa_dict[ID].append(line.strip().upper())

    for ID, seq_list in fa_dict.items():
        seq = ''.join(seq_list)
        fa_dict[ID] = seq

    return fa_dict


def parse_delta(mock_delta, ref_fa_dict, query_fa_dict):

    def get_comrev(seq):
        return seq.translate(comp_table)[::-1]

    with open(mock_delta) as f:
        for line in f:
            if line.startswith('>'):
                cols = line.split()
                ref_chrom = cols[0][1:]
                scaffold_name = cols[1]
                ref_seq, query_seq = ref_fa_dict[ref_chrom], query_fa_dict[scaffold_name]
                assert len(ref_seq) == int(cols[2])
                assert len(query_seq) == int(cols[3])
                aln_info = f.readline().split()
                ref_start, ref_end, scaffold_start, scaffold_end = [int(v) for v in aln_info[:4]]
                assert 'N' not in ref_seq[ref_start-1:ref_end]
                if scaffold_start < scaffold_end:
                    assert ref_seq[ref_start-1:ref_end] == query_seq[scaffold_start-1:scaffold_end]
                else:
                    assert ref_seq[ref_start-1:ref_end] == get_comrev(query_seq[scaffold_end-1:scaffold_start])
    print('Check Passed...')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_fa', help='fasta for reference genome')
    parser.add_argument('query_fa', help='fasta for query assembly')
    parser.add_argument('mock_delta', help='mocked delta file')
    args = parser.parse_args()

    ref_fa_dict = parse_fasta(args.ref_fa)
    query_fa_dict = parse_fasta(args.query_fa)
    parse_delta(args.mock_delta, ref_fa_dict, query_fa_dict)

if __name__ == '__main__':
    main()
