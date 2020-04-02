#!/usr/bin/python

import sys
import os
import random
import argparse
from os import listdir
from os.path import isfile, join
import construct_vcf_from_multifasta
from find_singleton_mutations import is_snp_line, vcfline_to_id


def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Script for comparison of large vcf file from GISAID SARS-CoV2 DB with given strain")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--ref', required = True, help='Path to reference (NC_045512.2 aka MN908947) file')
    required_args.add_argument('--genome', required = True, help='Path to your SARS-CoV2 strain')
    required_args.add_argument('--multivcf', required = True, help='Path to vcf with all snps')   
    required_args.add_argument('--out', required = True, help='Path to output vcf file')
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--minimap2', help='Path to minimap2 distribution')    
    optional_args.add_argument('--picard', help='Path to a folder with picard https://broadinstitute.github.io/picard/  script')
    optional_args.add_argument('--k8', help='Path to k8 binary, required to run paftools.js from minimap2 distribution')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def compare_vcfs(large_db, new_vcf):
    new_snps = {}
    snps_count = 0
    for line in open (new_vcf, "r"):
        if is_snp_line(line):
            new_snps[vcfline_to_id(line)] = line.strip()
            snps_count +=1
    print ("Total " + str(snps_count) + " variations")
    if snps_count > 100:
        print ("Seems that something went wrong - this is too much")
        exit()
    old_snps = {}
    for line in open (large_db, "r"):
        if is_snp_line(line):
            old_snps[vcfline_to_id(line)] = line.strip()   
    new_singletons = []
    for id in new_snps.keys():
        if not (id in old_snps.keys()):
            new_singletons.append(new_snps[id])
    print("There are " + str(len(new_singletons)) + " variations in your sample that are not present in DB")
    for line in new_singletons:
        print (line)

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    args.multifasta = args.genome
    construct_vcf_from_multifasta.construct_vcf(args)
    compare_vcfs(args.multivcf, args.out)
