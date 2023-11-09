#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-10-08 21:20


import argparse


def parse_assembly(assembly):
    
    n = 0
    frag_dict = dict()
    with open(assembly) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if line.startswith('>Chr'):
                frag = cols[0]
                number = cols[1]
                length = int(cols[2])
                frag_dict[number] = (frag, length)
            elif not line.startswith('>'):
                # one group should have at least two contigs
                if len(cols) < 2:
                    continue
                n += 1
                with open('group{}.txt'.format(n), 'w') as fout:
                    for number_ori in cols:
                        number = number_ori.strip('-')
                        if number in frag_dict:
                            frag, length = frag_dict[number]
                            fout.write('{}\tNA\t{}\n'.format(frag, length))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('assembly', help='input assembly file (output file "genome.final.assembly" of 3D-DNA)')
    args = parser.parse_args()

    parse_assembly(args.assembly)


if __name__ == '__main__':
    main()
