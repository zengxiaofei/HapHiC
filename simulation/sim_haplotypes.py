#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-04-20 17:27

import argparse
import numpy as np
import os
import collections
import scipy.stats

def parse_fasta(fasta):
    
    seq_list = list()
    
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                seq_list.append([ID, []])
            else:
                seq_list[-1][-1].append(line.strip().upper())

    for n, (ID, _list) in enumerate(seq_list):
        seq_list[n][-1] = ''.join(_list)

    return seq_list


def get_random_lists(seq_list, div, ploidy, weight, ts_tv_ratio, distribution, CV, seed):
    
    seq_len_list = [len(seq.replace('N', '')) for ID, seq in seq_list]
    total_seq_len = sum(seq_len_list)

    elements = [0, 1]
    coord_lists = list()
    changed_n_lists = list()
    for n in range(ploidy):
        coord_lists.append([])
        # bin size = 100kb
        bin_size=100000
        bin_num = total_seq_len//bin_size+1 if total_seq_len % bin_size else total_seq_len//bin_size
        probabilities_list = list()
        for v in scipy.stats.norm.rvs(loc=div, scale=CV*div, size=bin_num, random_state=seed*(n+1)):
            if v > 0:
                probabilities_list.append(v)
            else:
                probabilities_list.append(0)
        for m in range(bin_num):
            probabilities = [1-probabilities_list[m], probabilities_list[m]]
            np.random.seed(seed*(m+1)*(n+1))
            if m == bin_num-1 and total_seq_len % bin_size:
                bin_size = total_seq_len % bin_size
            coord_lists[-1].extend(np.random.choice(elements, bin_size, p=probabilities).tolist())
        changed_n_lists.append(coord_lists[-1].count(1))
    # substitution: A->TCG T->ACG C->ATG G->ATC (0, 1, 2)
    # insertion: -->ATCG (3, 4, 5, 6)
    # deletion: ATCG->- (7)

    weight_list = [float(w) for w in weight.split(',')]
    subw, insw, delw = [w/sum(weight_list) for w in weight_list]
    elements = [0, 1, 2, 3, 4, 5, 6, 7]
    r = ts_tv_ratio
    probabilities = [subw*r/(1+r), subw/(2*r+2), subw/(2*r+2), insw/4, insw/4, insw/4, insw/4, delw]
    change_lists = list()
    
    for n in range(ploidy):
        np.random.seed(seed*(n+1))
        change_lists.append(np.random.choice(elements, changed_n_lists[n], p=probabilities).tolist())
    return coord_lists, change_lists


def sim_haps(seq_list, coord_lists, change_lists, ploidy, fasta):

    sub_dict = {'A': 'GCT', 'T': 'CAG', 'C': 'TAG', 'G': 'ATC'}
    ins_dict = {3: 'A', 4: 'T', 5: 'C', 6: 'G'}

    prefix = os.path.basename(fasta).rsplit('.', 1)[0]
    
    allele_dict = collections.defaultdict(list)

    for n in range(ploidy):
        with open('{}_hap{}.fa'.format(prefix, n+1), 'w') as fout:
            nN = 0
            nbase = 0
            nchanged = 0
            for ID, seq in seq_list:
                # write ID
                fout.write('>{}_{}\n'.format(ID.rsplit('_', 1)[0], n+1))
                # write seq
                for coord, base in enumerate(seq):
                    if base != 'N':
                        if coord_lists[n][nbase]:
                            code = change_lists[n][nchanged]
                            if code in {0, 1, 2}:
                                newbase = sub_dict[base][code]
                                fout.write(newbase)
                                allele_dict[(ID, coord, base)].append((n, newbase))
                            elif code in {3, 4, 5, 6}:
                                newbase = base + ins_dict[code]
                                fout.write(newbase)
                                allele_dict[(ID, coord, base)].append((n, newbase))
                            else:
                                allele_dict[(ID, coord, base)].append((n, '-'))
                            nchanged += 1
                        # no change
                        else:
                            fout.write(base)
                        nbase += 1
                    # base is N
                    else:
                        fout.write(base)
                        nN += 1
                # write additional \n
                fout.write('\n')

    with open('allele_info.txt', 'w') as fout:
        haps = '\t'.join(['hap_{}'.format(p+1) for p in range(ploidy)])
        fout.write('Number\tChrom\tRef_coord\tRef_base\t{}\n'.format(haps))
        num = 0
        for ID, seq in seq_list:
            for coord, base in enumerate(seq):
                if (ID, coord, base) in allele_dict:
                    num += 1
                    haps_allele = [base]*ploidy
                    for p, newbase in allele_dict[(ID, coord, base)]:
                        haps_allele[p] = newbase
                    haps_allele = '\t'.join(haps_allele)
                    fout.write('{}\t{}\t{}\t{}\t{}\n'.format(num, ID, coord+1, base, haps_allele))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'fasta',
            help='input template genome of FASTA format')
    parser.add_argument(
            'div', type=float, 
            help='sequence divergence between simulated haplotypes and template genome, interval (0, 1)')
    parser.add_argument(
            'ploidy', type=int,
            help='how many haplotypes to simulate')
    parser.add_argument(
            '--weight', type=str, default='10,0.5,0.5',
            help='weight of substitutions, insertions and deletions, separated with commas, default: %(default)s')
    parser.add_argument(
            '--ts_tv_ratio', type=float, default=2.0,
            help='transition / transversion ratio of substitutions, default: %(default)s')
    parser.add_argument(
            '--distribution', type=str, choices={'normal', 'uniform'}, default='normal',
            help='distribution of sequence divergence (100kb bin), default: %(default)s')
    parser.add_argument(
            '--CV', type=float, default=0.5,
            help='coeffcient of variation of sequence divergence (100kb bin), default; %(default)s')
    parser.add_argument(
            '--seed', type=int, default=12345,
            help='seed for random processes, default: %(default)s')

    args = parser.parse_args()

    seq_list = parse_fasta(args.fasta)

    coord_lists, change_lists = get_random_lists(seq_list, args.div, args.ploidy, args.weight, args.ts_tv_ratio, args.distribution, args.CV, args.seed)
    
    sim_haps(seq_list, coord_lists, change_lists, args.ploidy, args.fasta)


if __name__ == '__main__':
    main()

