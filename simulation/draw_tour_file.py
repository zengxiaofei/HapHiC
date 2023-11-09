#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-07-27 13:54


import os
import argparse
import collections
import numpy as np
from decimal import Decimal
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def parse_fasta(fasta_file):

    ref_len_dict = collections.defaultdict(int)
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

    with open(tour_file) as f:
        for line in f:
            pass

    return line.split()


def analyze_sorted_ctgs(sorted_ctgs, ref_len_dict):
    
    # get the dominant source chromosome of the group
    chr_len_dict = collections.defaultdict(int)
    
    for ctg in sorted_ctgs:
        elements = ctg.split('_')
        source_chr = elements[0]
        chr_len_dict[source_chr] += int(elements[-1][:-1])

    dominant_chr = sorted(chr_len_dict.items(), key=lambda x: x[1])[-1][0]

    # get x and y of each line
    line_list = list()
    # used for CCC calculation
    coord_list = [[], []]

    accumulated_length = 0
    for ctg in sorted_ctgs:
        elements = ctg.split('_')
        source_chr = elements[0]
        if source_chr == dominant_chr:
            # get strand
            if elements[-2] == elements[-1][-1]:
                strand = '+'
            else:
                strand = '-'
            # get length
            length = int(elements[-1][:-1])
            # get x and y
            if strand == '+':
                line_list.append([accumulated_length+1, accumulated_length+length])
                coord_list[0].extend(list(range(accumulated_length+1, accumulated_length+length+1)))
            else:
                line_list.append([accumulated_length+length, accumulated_length+1])
                coord_list[0].extend(list(range(accumulated_length+1, accumulated_length+length+1))[::-1])
            line_list.append([int(elements[2]), int(elements[3])])
            coord_list[1].extend(list(range(int(elements[2]), int(elements[3])+1)))
            accumulated_length += length

    return line_list, coord_list, dominant_chr


def run_derange(sorted_ctgs, look_ahead):
    # prepare order.txt
    
    with open('order1.txt', 'w') as fout:
        out_list = list()
        for ctg in sorted_ctgs:
            if ctg[-1] == ctg.split('_')[4]:
                out_list.append(ctg.split('_')[1])
            else:
                out_list.append('-'+ctg.split('_')[1])
        out_list.append('0')
        fout.write('{}\n'.format(' '.join(out_list)))
    
    nctgs = len(sorted_ctgs)
    with open('order2.txt', 'w') as fout:
        out_list.clear()
        for ctg in sorted_ctgs:
            if ctg[-1] != ctg.split('_')[4]:
                out_list.append(str(nctgs-int(ctg.split('_')[1])+1))
            else:
                out_list.append('-'+str(nctgs-int(ctg.split('_')[1])+1))
        out_list.append('0')
        fout.write('{}\n'.format(' '.join(out_list)))
    
    # run derange2
    os.system('echo "" > derange.out && derange2 -S -L order1.txt {0} 1 1 2 derange.out && derange2 -S -L order2.txt {0} 1 1 2 derange.out'.format(look_ahead))
    
    costs = list()
    # parse derange.out
    with open('derange.out') as f:
        for line in f:
            if line.startswith('Total cost'):
                costs.append(int(float(line.split()[-1])))

    return min(costs)


def calculate_CCC(coord_list):
    
    y_pred, y_true = coord_list[0], coord_list[1]
    
    # Pearson product-moment correlation coefficient
    cor = np.corrcoef(y_true, y_pred)[0][1]
    # mean
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    # variance
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)
    # standard deviation
    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)
    # calculate CCC
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred)**2
    ccc = round(numerator / denominator, 8)
    return ccc


def draw_lineplot(line_list, dominant_chr, tour_file, program, N50, ccc, cost):
    
    group_name = os.path.basename(tour_file).split('.tour')[0]
    
    to_Mb = lambda x: Decimal(x/1000000)
    
    plt.figure(figsize=(2.5, 2.5))
    plt.rc('font', family='Arial')
    plt.xlabel('{} (Mb)'.format(group_name), fontsize=9)
    plt.ylabel('{} (Mb)'.format(dominant_chr), fontsize=9)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams.update({'font.size': 9})

    for n in range(len(line_list)//2):
        
        x = line_list[n*2]
        y = line_list[n*2+1]
        
        if x[1] > x[0]:
            color = '#E64B35'
            # color = '#d55e00'
        else:
            color = '#3C91BF'
            # color = '#0173b2'
        
        plt.plot(list(map(to_Mb, x)), list(map(to_Mb, y)), color=color, linewidth=2)
    
    title_list = list()
    stat_list = list()
    if ccc is not None:
        title_list.append("Lin's CCC = {:.2f}".format(ccc))
        stat_list.append(str(ccc))
    if cost is not None:
        title_list.append("derange cost = {}".format(cost))
        stat_list.append(str(cost))

    title = '\n'.join(title_list)

    if title_list:
        plt.legend(loc='upper center', title=title, frameon=False)
    
    plt.savefig('{}_{}_{}.pdf'.format(group_name, program, N50), bbox_inches='tight')
    
    if stat_list:
        print('{}\t{}\t{}\t{}'.format(group_name, program, N50, '\t'.join(stat_list)))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('tour', help='input tour file')
    parser.add_argument('program', help='program for ordering and orientation')
    parser.add_argument('N50', help='contig / scaffold N50')
    parser.add_argument('--CCC', default=False, action='store_true', help="calculate Lin's Concordance Correlation Coefficient")
    parser.add_argument('--derange', default=False, action='store_true', help='use derange2 to calculate cost, including inversions (weight 1), transpositions (weight 1) and transversions (weight 2)')
    parser.add_argument('--look_ahead', default=3, type=int, choices={2, 3, 4}, help='a parameter used in derange2, default: %(default)s')

    args = parser.parse_args()
    
    ref_len_dict = parse_fasta(args.fasta)

    sorted_ctgs = parse_tour(args.tour)
    
    line_list, coord_list, dominant_chr = analyze_sorted_ctgs(sorted_ctgs, ref_len_dict)

    if args.derange:
        cost = run_derange(sorted_ctgs, args.look_ahead)
    else:
        cost = None

    if args.CCC:
        ccc = calculate_CCC(coord_list)
    else:
        ccc = None

    draw_lineplot(line_list, dominant_chr, args.tour, args.program, args.N50, ccc, cost)


if __name__ == '__main__':
    main()

