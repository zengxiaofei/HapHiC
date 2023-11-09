#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-11-29 15:33

import argparse
import collections

def parse_stat(stat):

    def intepret_ID(ID):
        IDsplit = ID.split('_')
        if IDsplit[-1][-1] == '0':
            return int(IDsplit[6])//2
        else:
            assert IDsplit[-1][-1] == '1'
            return int(IDsplit[13])//2
    
    # (1) inter_homo_broken, (2) inter_homo_nonbroken
    # (3) inter_nonhomo_broken, (4) inter_nonhomo_nonbroken
    # (5) intra_chrom_broken, (6) intra_chrom_nonbroken
    # (7) nonchimeric_broken, (8) nonchimeric_nonbroken
    summary_dict = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0, 0, 0])
    
    # absolute distance to the actual error point
    # (1) inter_homo_broken, (2) inter_nonhomo_broken, (3) intra_chrom_broken
    break_point_distance_dict = collections.defaultdict(lambda: [[], [], []])
    
    with open(stat) as f:
        for line in f:
            cols = line.strip().split('\t')
            prog_N50 = tuple(cols[:2])
            ID = cols[2]
            if cols[3] == 'Inter_homo':
                if int(cols[4]) > 0:
                    summary_dict[prog_N50][0] += 1
                    for p in cols[5].split():
                        break_point_distance_dict[prog_N50][0].append(
                                abs(int(p) - intepret_ID(ID)))
                else:
                    assert int(cols[4]) == 0
                    summary_dict[prog_N50][1] += 1
            elif cols[3] == 'Inter_nonhomo':
                if int(cols[4]) > 0:
                    summary_dict[prog_N50][2] += 1
                    for p in cols[5].split():
                        break_point_distance_dict[prog_N50][1].append(
                                abs(int(p) - intepret_ID(ID)))
                else:
                    assert int(cols[4]) == 0
                    summary_dict[prog_N50][3] += 1
            elif cols[3] == 'Intra_chrom':
                if int(cols[4]) > 0:
                    summary_dict[prog_N50][4] += 1
                    for p in cols[5].split():
                        break_point_distance_dict[prog_N50][2].append(
                                abs(int(p) - intepret_ID(ID)))
                else:
                    assert int(cols[4]) == 0
                    summary_dict[prog_N50][5] += 1
            else:
                assert cols[3] == 'Non_chimeric'
                if int(cols[4]) > 0:
                    summary_dict[prog_N50][6] += 1
                else:
                    assert int(cols[4]) == 0
                    summary_dict[prog_N50][7] += 1

    return summary_dict, break_point_distance_dict


def output_summary(summary_dict, break_point_distance_dict):

    def get_dist_interval(dist):
        if dist <= 500:
            return '[0, 500]'
        elif dist <= 1000:
            return '(500, 1000]'
        elif dist < 5000:
            return '(1000, 5000]'
        elif dist < 10000:
            return '(5000, 10000]'
        elif dist < 50000:
            return '(10000, 50000]'
        elif dist < 100000:
            return '(50000, 100000]'
        else:
            return '>100000'
    
    with open('summary.txt', 'w') as f:
        for prog_N50, nctgs_list in summary_dict.items():
            prog, N50 = prog_N50
            f.write('{}\t{}\tInter_homo\tBroken\t{}\n'.format(prog, N50, nctgs_list[0]))
            f.write('{}\t{}\tInter_homo\tUnbroken\t{}\n'.format(prog, N50, nctgs_list[1]))
            f.write('{}\t{}\tInter_nonhomo\tBroken\t{}\n'.format(prog, N50, nctgs_list[2]))
            f.write('{}\t{}\tInter_nonhomo\tUnbroken\t{}\n'.format(prog, N50, nctgs_list[3]))
            f.write('{}\t{}\tIntra_chrom\tBroken\t{}\n'.format(prog, N50, nctgs_list[4]))
            f.write('{}\t{}\tIntra_chrom\tUnbroken\t{}\n'.format(prog, N50, nctgs_list[5]))
            f.write('{}\t{}\tNon_chimeric\tBroken\t{}\n'.format(prog, N50, nctgs_list[6]))
            f.write('{}\t{}\tNon_chimeric\tUnbroken\t{}\n'.format(prog, N50, nctgs_list[7]))

    with open('break_point_distance.txt', 'w') as f:
        for prog_N50, distance_list in break_point_distance_dict.items():
            prog, N50 = prog_N50
            for dist in distance_list[0]:
                dist_interval = get_dist_interval(dist)
                f.write('{}\t{}\tInter_homo\t{}\n'.format(prog, N50, dist_interval))
            for dist in distance_list[1]:
                dist_interval = get_dist_interval(dist)
                f.write('{}\t{}\tInter_nonhomo\t{}\n'.format(prog, N50, dist_interval))
            for dist in distance_list[2]:
                dist_interval = get_dist_interval(dist)
                f.write('{}\t{}\tIntra_chrom\t{}\n'.format(prog, N50, dist_interval))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stat', help='input correction_stat.txt')
    args = parser.parse_args()

    summary_dict, break_point_distance_dict = parse_stat(args.stat)
    output_summary(summary_dict, break_point_distance_dict)


if __name__ == '__main__':
    main()
