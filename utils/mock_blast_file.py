#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-07-01 18:08


import os
import argparse
import collections


def parse_fasta(fasta_file):
    
    ref_len_dict =  collections.defaultdict(int)
    with open(fasta_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                source_chr = ID.split('_')[0]
            else:
                ref_len_dict[source_chr] += len(line.strip())

    return ref_len_dict


def parse_tour(tour_file):
    
    last_line = ''
    with open(tour_file) as f:
        for line in f:
            last_line = line
    
    return last_line.split()


def parse_pre_sorted_ctgs(ctg_list, qname):

    total_len = 0

    # get dominant chromosome name (qname)
    len_dict = collections.defaultdict(int)
    
    for ctg in ctg_list:
        elements = ctg.split('_')
        source_chr = elements[0]
        len_dict[source_chr] += int(elements[-1][:-1])
    
    sname = sorted(len_dict.items(), key=lambda x: x[1])[-1][0]

    with open('blast_{}_{}.out'.format(qname, sname), 'w') as fout:
        for ctg in ctg_list:
            elements = ctg.split('_')
            if elements[-2] == elements[-1][-1]:
                strand = '+'
            else:
                strand = '-'
            length = int(elements[-1][:-1])
            
            if elements[0] == sname:
                # query
                qstart = total_len + 1
                qend = total_len + length

                # subject
                sstart, send = int(elements[2]), int(elements[3])
                if strand == '-':
                    sstart, send = send, sstart
                fout.write('{}\t{}\t100\t{}\t0\t0\t{}\t{}\t{}\t{}\t0\t10000\n'.format(
                    qname, sname, length, qstart, qend, sstart, send))
            total_len += length

    return sname, total_len


def generate_sizes(length, name, file_name):
    
    with open('{}.sizes'.format(file_name), 'w') as fout:
        fout.write('{}\t{}\n'.format(name, length))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('tour', help='input tour file')

    args = parser.parse_args()

    ref_len_dict = parse_fasta(args.fasta)

    pre_sorted_ctg_list = parse_tour(args.tour)

    qname = os.path.splitext(os.path.basename(args.tour))[0]

    sname, total_len = parse_pre_sorted_ctgs(pre_sorted_ctg_list, qname)

    generate_sizes(ref_len_dict[sname], sname, 'subject')

    generate_sizes(total_len, qname, 'query')

    os.system('python3 -m jcvi.graphics.blastplot blast_{}_{}.out --qsizes query.sizes --ssizes subject.sizes --style whitegrid'.format(qname, sname))

if __name__ == '__main__':
    main()

