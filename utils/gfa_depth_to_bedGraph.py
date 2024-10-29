#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-10-16 10:40

import argparse
import gzip
import re
import sys
import math

# `bedGraphToBigWig` in kentUtils (https://github.com/ucscGenomeBrowser/kent) 
# is necessary for the downstream format conversion:
# `$ bedGraphToBigWig in.bedGraph chrom.sizes out.bw`

def parse_gfa(gfas, depth_tag):

    depth_dict = dict()
    depth_pattern = re.compile(r'.+{}:[if]:([\d.]+)'.format(depth_tag))

    for gfa in gfas:

        if gfa.endswith('.gz'):
            fopen = gzip.open
        else:
            fopen = open

        with fopen(gfa, 'rt') as f:
            for line in f:
                if not line.startswith('S\t'):
                    continue
                segment = line.split()[1]
                depth_match = depth_pattern.match(line)
                if not depth_match:
                    raise Exception('Cannot find the read depth for segment {}'.format(segment))
                depth = int(depth_match.groups()[0]) + 1
                depth_dict[segment] = depth

    return depth_dict


def parse_agp(agp, scale):

    ctg_list = []
    total_len = 0
    with open(agp) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.split()
            if cols[4] == 'W':
                ctg, ctg_len = cols[5], int(cols[2]) - int(cols[1]) + 1
                total_len += ctg_len
                ctg_list.append((ctg, ctg_len))

    estimated_scale = 1 + total_len // 2100000000
    if scale is None:
        assembly_size = total_len // estimated_scale
        print('Scale is estimated to be {}, and the assembly size is {}'.format(
            estimated_scale, assembly_size), file=sys.stderr)
        return ctg_list, estimated_scale, assembly_size 

    if scale != estimated_scale:
        print('Warning: designated scale ({}) differs from the estimated one ({})'.format(
            scale, estimated_scale), file=sys.stderr)
    return ctg_list, scale, total_len // scale


def output_bedGraph(ctg_list, scale, assembly_size, depth_dict):

    BIN_SIZE = 100000

    with open('chrom.sizes', 'w') as f:
        print('Writing assembly size into `chrom.sizes`...', file=sys.stderr)
        f.write('assembly\t{}\n'.format(assembly_size))

    with open('in.bedGraph', 'w') as f:
        print('Writing depths of each contigs into `in.bedGraph`...', file=sys.stderr)
        f.write('track type=bedGraph\n')
        accumulated_start = 0
        for ctg, ctg_len in ctg_list:
            for n in range(math.ceil(ctg_len / BIN_SIZE)):
                start = accumulated_start + n * BIN_SIZE
                end = accumulated_start + (n + 1) * BIN_SIZE if (n + 1) * BIN_SIZE < ctg_len else accumulated_start + ctg_len
                f.write('assembly\t{}\t{}\t{}\n'.format(
                    start // scale,
                    end // scale,
                    depth_dict[ctg]))
            accumulated_start += ctg_len

    print('Finished. Please execute the command: `$ bedGraphToBigWig in.bedGraph chrom.sizes out.bw`', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='the AGP file output by scaffolders like HapHiC and YaHS')
    parser.add_argument('gfa', nargs='+', help='the GFA file(s) output by assemblers, gzipped files are acceptable')
    parser.add_argument('--depth_tag', default='rd', help='the tag name for read depth, default: %(default)s')
    parser.add_argument('--scale', default=None, type=int, help='assembly scale, default: %(default)s (estimated automatically)') 
    args = parser.parse_args()

    depth_dict = parse_gfa(args.gfa, args.depth_tag)
    ctg_list, scale, assembly_size = parse_agp(args.agp, args.scale)
    output_bedGraph(ctg_list, scale, assembly_size, depth_dict)


if __name__ == '__main__':
    main()

