#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2023-02-13 10:24

import argparse
import os
import re


def interpret_time(time_string):
    
    time_list = time_string.split(':')
    assert len(time_list) == 3
    return int(time_list[0])*60*60 + int(time_list[1])*60 + int(time_list[2])


def statistics(jobids, ndays):
    
    peakmem = 0
    accumulated_cput = 0
    accumulated_wallt = 0

    jobids = [jobid.split('.')[0] for jobid in jobids]
    if len(jobids) != len(set(jobids)):
        raise RuntimeError('Find repeated JobID!')

    for jobid in jobids:
        with os.popen('tracejob {} -n {} 2>&1'.format(jobid, ndays)) as f:
            for line in f:
                match = re.match(r'.+Exit_status=(\d).+resources_used.cput=([\w:]+) resources_used.mem=(\d+)kb.+resources_used.walltime=([\d:]+)', line)
                if match:
                    status, cput, mem, wallt = match.groups()
                    if status != '0':
                        raise RuntimeError('Exit status != 0 ({}, Exit_status={})'.format(jobid, status))
                    mem = int(mem)
                    accumulated_cput += interpret_time(cput)
                    accumulated_wallt += interpret_time(wallt)
                    if mem > peakmem:
                        peakmem = mem
                elif "Couldn't find Job Id" in line:
                    raise RuntimeError(line)
    
    print('Wall time = {} min\nCPU time = {} min\nPeak memory = {} GiB'.format(
        round(accumulated_wallt/60, 2), round(accumulated_cput/60, 2), round(peakmem/(1024*1024), 2)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('jobids', nargs='+', help='the jobids that need performance statistics')
    parser.add_argument('--ndays', type=int, default=100, help='number of days in the past to look for job, default: %(default)s')
    args = parser.parse_args()

    statistics(args.jobids, args.ndays)

if __name__ == '__main__':
    main()
