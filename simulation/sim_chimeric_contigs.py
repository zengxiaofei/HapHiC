#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-03-27 19:37


import argparse
import scipy.stats
from itertools import combinations
import random
import scipy.stats

def parse_fasta(fasta):

    contig_dict = dict()

    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                ID = line.split()[0][1:]
                contig_dict[ID] = list()
            else:
                contig_dict[ID].append(line.strip())
    
    for ID, seq_list in contig_dict.items():
        contig_dict[ID] = ''.join(seq_list)

    return contig_dict


def sim_chimeric_ctgs(contig_dict, chimeric_ratio, inner_chrom_weight, inter_homo_weight, inter_nonhomo_weight, seed):
    
    total_ctgs = len(contig_dict)
    
    # define chimeric_ratio = sim_n_pairs / total_ctgs
    sim_n_pairs = int(total_ctgs * chimeric_ratio)

    weight_sum = inner_chrom_weight + inter_homo_weight + inter_nonhomo_weight
    
    normalized_homo_weight = inter_homo_weight / weight_sum
    normalized_nonhomo_weight = inter_nonhomo_weight / weight_sum
    normalized_inner_weight = inner_chrom_weight / weight_sum

    homo_n_pairs = int(sim_n_pairs * normalized_homo_weight)
    nonhomo_n_pairs = int(sim_n_pairs * normalized_nonhomo_weight)
    inner_n_pairs = int(sim_n_pairs * normalized_inner_weight)

    total_n_pairs = homo_n_pairs + nonhomo_n_pairs + inner_n_pairs
    
    ctg_list = list(contig_dict.keys())
    random.seed(seed)
    random.shuffle(ctg_list)

    modified_ctg_set = set()
    
    sampled_inter_homo = list()
    sampled_inter_nonhomo = list()
    sampled_inner_chrom = list()

    random.seed(seed)
    random_ctg_list1 = list(contig_dict.keys())
    random.shuffle(random_ctg_list1)
    random.seed(seed*2)
    random_ctg_list2 = list(contig_dict.keys())
    random.shuffle(random_ctg_list2)

    for ctg_1 in random_ctg_list1:
            
        if homo_n_pairs + nonhomo_n_pairs + inner_n_pairs == 0:
            break

        if ctg_1 in modified_ctg_set or 'collapsed' in ctg_1:
            continue
        
        for ctg_2 in random_ctg_list2:
        
            if ctg_1 == ctg_2 or ctg_2 in modified_ctg_set or 'collapsed' in ctg_2:
                continue
            
            ctg_1_source = ctg_1.split('_')[:2]
            ctg_2_source = ctg_2.split('_')[:2]
            
            if ctg_1_source != ctg_2_source and ctg_1_source[0] == ctg_2_source[0]:
                if homo_n_pairs:
                    sampled_inter_homo.append((ctg_1, ctg_2))
                    modified_ctg_set.add(ctg_1)
                    modified_ctg_set.add(ctg_2)
                    homo_n_pairs -= 1
                    break
            elif ctg_1_source != ctg_2_source:
                if nonhomo_n_pairs:
                    sampled_inter_nonhomo.append((ctg_1, ctg_2))
                    modified_ctg_set.add(ctg_1)
                    modified_ctg_set.add(ctg_2)
                    nonhomo_n_pairs -= 1
                    break
            else:
                if inner_n_pairs:
                    sampled_inner_chrom.append((ctg_1, ctg_2))
                    modified_ctg_set.add(ctg_1)
                    modified_ctg_set.add(ctg_2)
                    inner_n_pairs -= 1
                    break

    # (0) ctg_1 left/right, (1) ctg_2 left/right,
    # (2) ctg_1_half revcom/not, (3) ctg_2_half revcom/not
    # (4) ctg_1_half ctg_2_half ordering 
    order_list = scipy.stats.bernoulli.rvs(0.5, size=total_n_pairs*5, random_state=seed).tolist()

    return sampled_inner_chrom, sampled_inter_homo, sampled_inter_nonhomo, order_list


def output_fasta(out_file, contig_dict, sampled_inner_chrom, sampled_inter_homo, sampled_inter_nonhomo, order_list):

    com_dict = {
            'a': 't', 'A': 'T', 'c': 'g', 'C': 'G',
            't': 'a', 'T': 'A', 'g': 'c', 'G': 'C',
            'n': 'n', 'N': 'N'
            }

    def revcom(seq):
        return ''.join([com_dict[base] for base in seq[::-1]])

    def write_lines(sampled_list, chimeric_type, n):
        
        for ctg_1, ctg_2 in sampled_list:
            ctg_1_seq = contig_dict[ctg_1]
            ctg_1_len = len(ctg_1_seq)
            ctg_2_seq = contig_dict[ctg_2]
            ctg_2_len = len(ctg_2_seq)
            # split contigs
            ctg_1_left = ctg_1_seq[:ctg_1_len//2]
            ctg_1_right = ctg_1_seq[ctg_1_len//2:]
            ctg_2_left = ctg_2_seq[:ctg_2_len//2]
            ctg_2_right = ctg_2_seq[ctg_2_len//2:]
            # change the order of randomly
            if order_list[5*n]:
                ctg_1_chimeric = ctg_1_left
                ctg_1_frag = ctg_1_right
            else:
                ctg_1_chimeric = ctg_1_right
                ctg_1_frag = ctg_1_left

            if order_list[5*n+1]:
                ctg_2_chimeric = ctg_2_left
                ctg_2_frag = ctg_2_right
            else:
                ctg_2_chimeric = ctg_2_right
                ctg_2_frag = ctg_2_left
            
            if order_list[5*n+2]:
                ctg_1_chimeric = revcom(ctg_1_chimeric)
            if order_list[5*n+3]:
                ctg_2_chimeric = revcom(ctg_2_chimeric)
           
            if order_list[5*n+4]:
                ctg_1_chimeric, ctg_2_chimeric = ctg_2_chimeric, ctg_1_chimeric
            
            fout.write('>{}_{}_{}_{}\n{}{}\n'.format(
                ctg_1, ctg_2, chimeric_type, ''.join([str(num) for num in order_list[5*n:5*n+5]]), ctg_1_chimeric, ctg_2_chimeric))

            fout.write('>{}_frag\n{}\n'.format(ctg_1, ctg_1_frag))
            fout.write('>{}_frag\n{}\n'.format(ctg_2, ctg_2_frag))
            
            chimeric_ctgs.add(ctg_1)
            chimeric_ctgs.add(ctg_2)
            n += 1

        return n


    chimeric_ctgs = set()

    with open(out_file, 'w') as fout:
        
        n = 0
        n = write_lines(sampled_inner_chrom, 'chimeric_inner_chrom', n)
        n = write_lines(sampled_inter_homo, 'chimeric_inter_homo', n)
        n = write_lines(sampled_inter_nonhomo, 'chimeric_inter_nonhomo', n)

        for ctg, seq in contig_dict.items():
            if ctg not in chimeric_ctgs:
                fout.write('>{}\n{}\n'.format(ctg, seq))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file of contigs')
    parser.add_argument(
            '--chimeric_ratio', type=float, default=0.05,
            help='final chimeric contigs simulated / total contigs, default: %(default)s')
    parser.add_argument(
            '--inner_chrom_weight', type=float, default=0.1, 
            help='weight of simulating inner-chromosome chimeric contigs, default: %(default)s')
    parser.add_argument(
            '--inter_homo_weight', type=float, default=0.7, 
            help='weight of simulating inter-homologous-chromosome chimeric contigs, default: %(default)s')
    parser.add_argument(
            '--inter_nonhomo_weight', type=float, default=0.2, 
            help='weight of simulating inter-nonhomologous-chromosome chimeric contigs, default: %(default)s')
    parser.add_argument(
            '--seed', type=int, default=12345, 
            help='seed for contig random combination, default: %(default)s')
   
    args = parser.parse_args()

    contig_dict = parse_fasta(args.fasta)


    sampled_inner_chrom, sampled_inter_homo, sampled_inter_nonhomo, order_list = sim_chimeric_ctgs(
            contig_dict, args.chimeric_ratio, args.inner_chrom_weight, 
            args.inter_homo_weight, args.inter_nonhomo_weight, args.seed)


    out_file = '{}_chimeric_{}_{}_{}_{}.fa'.format(
            args.fasta.rsplit('.', 1)[0], args.chimeric_ratio,
            args.inner_chrom_weight, args.inter_homo_weight, args.inter_nonhomo_weight
            )

    output_fasta(out_file, contig_dict, sampled_inner_chrom, sampled_inter_homo, sampled_inter_nonhomo, order_list)


if __name__ == '__main__':
    main()

