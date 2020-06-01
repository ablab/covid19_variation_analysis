#!/usr/bin/python

import sys
import os
import random
import argparse
from os import listdir
from os.path import isfile, join
import construct_vcf_from_multifasta
from find_singleton_mutations import is_snp_line, vcfline_to_id, is_valid_indel


def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Script for describing indels")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--multivcf', required = True, help='Path to vcf with all snps')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def construct_indels(large_db):
    old_snps = {}
    for line in open (large_db, "r"):
        if is_snp_line(line):
            if is_valid_indel(line):
                id = vcfline_to_id(line)
                if not id in old_snps.keys():
                    old_snps[id] = []
                old_snps[id].append(line.strip())
    sorted = {}
    for id in old_snps.keys():    
        if len (old_snps[id]) > 1:
            arr = id.split('_')
            pos = arr[0]
            if (len(arr[1]) - len(arr[2])) % 3 == 0 :
                print (id + " " + str(len(old_snps[id])))
    


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    construct_indels(args.multivcf)
