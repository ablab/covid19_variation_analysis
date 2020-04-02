#!/usr/bin/python

import sys
import os
import random
import argparse
from os import listdir
from os.path import isfile, join

minimap2 = ""
picard = ""
k8 = ""
#also cleanups header from spaces and some other bad chars 
def split_fasta(filename, outputdir):
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
   
    outFile = None
    for line in open(filename):
        if line[0] == '>':
            if outFile:
                outFile.close()
            line = line.replace(" ","_").replace("|","_").replace("/","_")
            outFile = open(os.path.join(outputdir, line[1:].strip() + '.fa'), 'w')
        if outFile:
            outFile.write(line)
    if outFile: # if filename is empty
        outFile.close()

def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Script for construction of large vcf file from GISAID SARS-CoV2 DB. Output: all_sarscov2.vcf")
    parser._action_groups.pop()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--ref', required = True, help='Path to reference (NC_045512.2 aka MN908947) file')
    required_args.add_argument('--multifasta', required = True, help='Path to multifasta with SARS-CoV2 strains')
    optional_args = parser.add_argument_group('optional arguments')
    optional_args.add_argument('--minimap2', help='Path to minimap2 distribution')    
    optional_args.add_argument('--picard', help='Path to a folder with picard https://broadinstitute.github.io/picard/  script')
    optional_args.add_argument('--k8', help='Path to k8 binary, required to run paftools.js from minimap2 distribution')   
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()



def construct_vcf(ref, merged):
    splitted = "splitted"
    vcfs = "vcfs"
    if not os.path.isdir(splitted):
        os.mkdir(splitted)
    for f in listdir(splitted):
        os.remove(os.path.join(splitted, f))
    if not os.path.isdir(vcfs):
        os.mkdir(vcfs)
    for f in listdir(vcfs):
        os.remove(os.path.join(vcfs, f))

    split_fasta (merged, splitted)

    for f in listdir(splitted):
        l = (os.path.join(minimap2, "minimap2") + " -c --cs " + ref + " " + os.path.join(splitted, f) + " 2>/dev/null | sort -k6,6 -k8,8n | " + k8 + " " +  os.path.join(minimap2, "misc","paftools.js") + " call -L20000 -f " + ref + " -  > " + os.path.join(vcfs, os.path.basename(f) + ".vcf"))
        os.system (os.path.join(minimap2, "minimap2") + " -c --cs " + ref + " " + os.path.join(splitted, f) + " 2>/dev/null | sort -k6,6 -k8,8n | " + k8 + " " +  os.path.join(minimap2, "misc","paftools.js") + " call -L20000 -f " + ref + " -  > " + os.path.join(vcfs, os.path.basename(f) + ".vcf"))
    vcflist = open("vcfs.list", "w")
    for f in listdir(vcfs):
        if f.find("_bat_")== -1 and f.find("pangolin") == -1:
            vcflist.write(os.path.join(vcfs, f) + "\n")
    vcflist.close()
    os.system("java -jar " + os.path.join(picard, "picard.jar") + " MergeVcfs I=vcfs.list O=all_sarscov2.vcf")

#if len(sys.argv) != 3:
#    print "Usage: " + sys.argv[0] + " ref multifasta "
#    exit()
args = parse_args(sys.argv[1:])
if args.minimap2:
    minimap2= args.minimap2
if args.k8:
    k8=args.k8
if args.picard:
    picard = args.picard
construct_vcf(args.ref, args.multifasta)
