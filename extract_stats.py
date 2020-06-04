#!/usr/bin/python

import sys
import os
import random
import argparse
from os import listdir
from os.path import isfile, join
import construct_vcf_from_multifasta
from find_singleton_mutations import is_snp_line, vcfline_to_id
import subprocess

report_dir = "/Bmo/dantipov/covid/run_01_06/"



def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Support script to grep read-assembly stats generation")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
#    required_args.add_argument('--ref', required = True, help='Path to reference (NC_045512.2 aka MN908947) file')
    required_args.add_argument('--metric', help='metric to grep')

    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--cutoff', help='cutoff more or equal than X')    
    optional_args.add_argument('--tech', help='tech to grep')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def run_request(args):
    grep_line = "grep " +args.metric + " " + report_dir + "/*/*.report "
    if args.tech:
        grep_tech = "grep " + args.tech + " " + report_dir + "/*/*.report "  
        reports = set()
        out_tech = subprocess.check_output(grep_tech , shell=True).split('\n')
        for  line in out_tech:
            reports.add(line.split('\t')[0].split(':')[0])
#        print reports
    output = subprocess.check_output(grep_line , shell=True).split('\n')
#    print len(output)
#    print output    
#    print output
    cutoff = -1
    if args.cutoff:
        cutoff = args.cutoff
    filtered = []
    for line in output:
        arr = line.split('\t')
        report = arr[0].split(':')[0]
        if args.tech and report not in reports:
            continue
#        print arr
        paf = arr[0].split(':')[0].replace('report', 'paf')
        if len(arr) > 1 and int(arr[1]) >= int(cutoff):
            paf_line = "wc -l " + paf
            paf_check =  subprocess.check_output(paf_line , shell=True)
            paf_lines = int(paf_check.split()[0])
            if paf_lines == 1:
                print line
#            print paf

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    run_request(args)
