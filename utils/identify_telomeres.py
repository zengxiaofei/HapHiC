#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn

import argparse
import sys


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def parse_fasta(fasta):
    seq_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                seq_dict[ID] = []
            else:
                seq_dict[ID].append(line.strip())

    for ID, seq_list in seq_dict.items():
        seq_dict[ID] = ''.join(seq_list)

    return seq_dict


def identify_telomeres(seq_dict, repeat, window, min_count, min_length):
    repeat_rev_comp = revcomp(repeat)
    results = []

    for ctg_id, seq in seq_dict.items():
        seq_len = len(seq)
        if seq_len < min_length:
            continue

        types = []
        if seq[:window].count(repeat) >= min_count:
            types.append('telomeric_at_5')
        if seq[-window:].count(repeat_rev_comp) >= min_count:
            types.append('telomeric_at_3')

        results.append((ctg_id, seq_len, types))

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Identify telomeric repeats at the ends of sequences in a FASTA file.')

    parser.add_argument(
        'fasta', help='input FASTA file')
    parser.add_argument(
        '--repeat', default='CCCTAAA',
        help='telomeric repeat motif, default: %(default)s')
    parser.add_argument(
        '--window', type=int, default=1000,
        help='window size (bp) at each end to search for repeats, default: %(default)s')
    parser.add_argument(
        '--min-count', type=int, default=5,
        help='minimum number of repeat occurrences within the window, default: %(default)s')
    parser.add_argument(
        '--min-length', type=int, default=0,
        help='minimum sequence length to analyze, default: %(default)s')
    parser.add_argument(
        '--output', '-o', help='output file for telomere annotation (optional)')

    args = parser.parse_args()

    seq_dict = parse_fasta(args.fasta)

    results = identify_telomeres(seq_dict, args.repeat, args.window, args.min_count, args.min_length)

    n_total = len(seq_dict)
    n_analyzed = len(results)
    n_skipped = n_total - n_analyzed

    n_telomeric_at_5 = 0
    n_telomeric_at_3 = 0
    n_both = 0
    n_at_5_only = 0
    n_at_3_only = 0

    if args.output:
        fout = open(args.output, 'w')
    else:
        fout = None

    print('Per-sequence results:')
    print('-' * 60)
    print('{:<30} {:>10} {:>20}'.format('Sequence ID', 'Length(bp)', 'Telomere Type'))
    print('-' * 60)

    for ctg_id, seq_len, types in results:
        type_str = ','.join(types) if types else 'None'
        print('{:<30} {:>10} {:>20}'.format(ctg_id, seq_len, type_str))

        if fout:
            fout.write('{}\t{}\n'.format(ctg_id, ','.join(types) if types else 'None'))

        if types:
            if 'telomeric_at_5' in types and 'telomeric_at_3' in types:
                n_both += 1
                n_at_5_only += 1
                n_at_3_only += 1
            elif 'telomeric_at_5' in types:
                n_telomeric_at_5 += 1
                n_at_5_only += 1
            elif 'telomeric_at_3' in types:
                n_telomeric_at_3 += 1
                n_at_3_only += 1

    if fout:
        fout.close()

    n_any = n_telomeric_at_5 + n_telomeric_at_3 + n_both

    print('-' * 60)
    print()
    print('Summary:')
    print('  Total sequences in FASTA:          {}'.format(n_total))
    if args.min_length > 0:
        print('  Sequences >= {} bp (analyzed):     {}'.format(args.min_length, n_analyzed))
        print('  Sequences <  {} bp (skipped):      {}'.format(args.min_length, n_skipped))
    else:
        print('  Sequences analyzed:                {}'.format(n_analyzed))
    print('  Sequences with any telomere:       {}'.format(n_any))
    print('    telomeric_at_5 (only):           {}'.format(n_telomeric_at_5))
    print('    telomeric_at_3 (only):           {}'.format(n_telomeric_at_3))
    print('    both telomeric_at_5 and _at_3:   {}'.format(n_both))
    print('  Total with telomeric_at_5:         {}'.format(n_at_5_only))
    print('  Total with telomeric_at_3:         {}'.format(n_at_3_only))


if __name__ == '__main__':
    main()
