
#!/usr/bin/python


import sys
import os
import random
from os import listdir
from os.path import isfile, join

from joblib import Parallel, delayed
import csv
from StringIO import StringIO
import subprocess


#constants
low_coverage = 10
samtools = "samtools"
bcftools = "~/other_tools/bcftools/bcftools"
minimap2 = "~/other_tools/minimap2/minimap2"
ref = "/home/dantipov/scripts/covid/MN908947.3.fasta"
data_pref = "/Bmo/dantipov/data/sars_qc/"
#work_pref = "/Bmo/dantipov/covid/run_old/"
work_pref =  "/Bmo/dantipov/covid/run_01_06/"
sra_tools =  "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"
fixed_paftools = "/home/dantipov/scripts/covid/paftools_N.js"
k8 = "/home/dantipov/scripts/covid/k8-0.2.4/k8-Linux"
big_fasta ="/home/dantipov/scripts/covid/sequences_2020-05-25_16-36.fasta"
#big_list ="/home/dantipov/scripts/covid/sra_nextstrain_old.txt"
big_list = "/home/dantipov/scripts/covid/SARS_RunTable_30.05.txt"
work_tmp = "/Bmo/dantipov/tmp/"
nextstrain_good = "/Bmo/dantipov/tmp/good.ids"
nextstrain_bad = "/Bmo/dantipov/tmp/bad.ids"

def download_sample(id, outdir):
    pr_line = sra_tools + "prefetch --max-size 40000000 " + id
    print pr_line
    os.system(pr_line)
    fq_dump_line = sra_tools + "fastq-dump --gzip --split-files " + id + " -O " + outdir
    print fq_dump_line
    os.system(fq_dump_line)

#returns number of sequences in resulting fasta
def extract_reference(srr, gisaid, workdir):
    outfile = join(workdir, srr +".fasta")
    line = samtools + " faidx "  + big_fasta + " \""+ gisaid +"\" > " + outfile
    res = os.system(line) >> 8
    print line
    if res != 0:
        return 0
    grep_line = ("grep \">\" " + outfile + " | wc -l ")
    output = subprocess.check_output(grep_line , shell=True)
    print int(output)
    return int(output)

def extract_nextstrain_id(srr, gisaid, workdir):
    end_str = "/20"
    pos = gisaid.find(end_str)
    if pos == -1:
        return -1
#/2020 length
    pos += 5
    cropped = gisaid[:pos]
    arr = cropped.split("/")
    if len(arr) < 3:
        return -1
    if arr[-3] == "DK":
        arr[-3] = "Denmark"
        arr[-2] = arr[-2].replace("ALAB-HH-", "ALAB-HH")
    new_str = arr[-3] + "/" + arr[-2] + "/" + arr[-1]
    print srr + " " + new_str
    if extract_reference(srr, new_str, workdir)  != 1:
        return -1
    return 0

def extract_all_nextstrain(ids):
    good = open(nextstrain_good, "w")
    bad = open(nextstrain_bad, "w")
    count = 0
    for str in ids:
        if extract_nextstrain_id(str) !=0:
            bad.write(str[1]+ "\n")
        else:
            good.write(str[1]+ "\n")
        count +=1
#        if count == 100:
#            break
def run_sample(sample_descr):
    srr = sample_descr[0]
    gisaid = sample_descr[-1]
    workdir = join(work_pref, srr)
    if not os.path.isdir (workdir):
        os.mkdir(workdir)
    res = extract_nextstrain_id(srr, gisaid, workdir)
    print res
    if res != 0:
        return
    process_sample(sample_descr, workdir)

def process_list(inputlist, outdir):
    random.seed(239)
    ids = []
    for name in inputlist:
        rand = random.randint(0, 2)
#        if rand != 0:
#            continue
        id = name[0]
        ids.append(id)    
    Parallel(n_jobs=10)(delayed(download_sample)(id, outdir)
    for id  in ids)
#        exit()



def get_minimap_illumina_paired(srr_id, workdir):
    read1 = join(data_pref, srr_id+ "_1.fastq.gz")
    read2 = join(data_pref, srr_id+ "_2.fastq.gz")
    res = minimap2 + " -ax sr " + ref + " " + read1 +" " + read2  + " -t 10 > " + join(workdir, srr_id+".sam")
    return res

def get_minimap_illumina_single(srr_id, workdir):
    read1 = join(data_pref, srr_id+ "_1.fastq.gz")
    res = minimap2 + " -ax sr " + ref + " " + read1  + " -t 10 > " + join(workdir, srr_id+".sam")
    return res

def get_minimap_nanopore(srr_id, workdir):
    read1 = join(data_pref, srr_id+ "_1.fastq.gz")
    res = minimap2 + " -ax map-ont  " + ref + " " + read1  + " -t 10 > " + join(workdir, srr_id+".sam")
    return res

def run_ion_torrent():
    return
#TODO
def parse_paf ( paf_file):
    f = open (paf_file, 'r')
    res = {}
    for line in f:
        arr = line.split()
#    123 678
        res['qlen'] = int(arr[1])
        res['qstart'] = int(arr[2])
        res['qend'] = int(arr[3])
        res['rlen'] = int(arr[6])
        res['rstart'] = int(arr[7])
        res['rend'] = int(arr[8])
        res['name'] = arr[0].split('|')[0]
        break

    return res

def parse_vcf (vcf_file):
    f = open (vcf_file, 'r')
    res = []
    for line in f:
        arr = line.split()
        if len(arr) > 5 and arr[1].isdigit():
            res.append([int(arr[1]), arr[4], arr[3]])
    return res

def parse_depth (depth_file):
    res = {}
    pos = []
    for line in open (depth_file,'r'):
        arr = line.split()
        if len(arr) > 2:
            pos.append(int(arr[2]))
    res['trailing'] = 0
    k = 0
    while pos[k] == 0:
        k += 1
    res['leading'] = k
    k = len(pos) - 1
    while pos[k] == 0:
        k -= 1
    res['trailing'] = len(pos) - 1 - k
    res['zero_cov'] = set()
    res['low_cov'] = set()
    for k in range(res['leading'], len(pos) - res['trailing']):
        if pos[k] < low_coverage:
            res['low_cov'].add(k + 1)
        if pos[k] == 0:
            res['zero_cov'].add(k + 1)
    print res
    return res

def out_array(arr):
    res = str(len(arr)) + '\t'
    if len(arr) == 0:
        return res
    for i in range(0, len(arr) -1):
        res += str(arr[i])+ ','
    res += str(arr[len(arr) - 1])
    return res

def analyze_sample(srr_id, workdir, tech):
    paf_file = os.path.join(workdir, srr_id +".paf")
    vcf_file = os.path.join(workdir, srr_id +".vcf")
    depth_file = os.path.join(workdir, srr_id +".depth")
    report_file =  os.path.join(workdir, srr_id +".report")
    depth_dict = parse_depth(depth_file)
    snps = parse_vcf(vcf_file)
    paf = parse_paf(paf_file)
    extra_leading = depth_dict['leading'] - paf['rstart']    
    if extra_leading < 0:
        extra_leading = 0
    extra_trailing = depth_dict['trailing'] - paf['rlen'] + paf['rend']    
    if extra_trailing < 0:
        extra_trailing = 0
    zero_covered_snps = []
    low_covered_snps = []
    supported_0 = set()
    for snp in snps:
        if snp[1] != "N" and snp[1] != "n":
            if snp[0] in depth_dict['zero_cov']:
                zero_covered_snps.append(snp[0])
            elif snp[0] in depth_dict['low_cov']:
                low_covered_snps.append(snp[0])
        else:
            supported_0.add(snp[0])
    deleted = set()
    for snp in snps:
        if len(snp[2]) != 1:
            for i in range(1, len (snp[2])):
                print "excluding " + str(i + snp[0])
                deleted.add(snp[0] + i)
    unsupported_0 =  depth_dict['zero_cov'] - supported_0 - deleted
#    print unsupported_0
    report_f = open(report_file, "w")
    report_f.write("name\t" + paf['name'] + "\n")
    report_f.write("tech\t" + tech + "\n")
    report_f.write("unreported_leading_N\t" + str(extra_leading) + "\n")
    report_f.write("unreported_trailing_N\t" + str(extra_trailing) + "\n")
    report_f.write("ZERO_covered_snps\t" + out_array(zero_covered_snps)  + "\n")
    report_f.write("low_covered_snps\t" + out_array(low_covered_snps) + "\n")
    report_f.write("missing_N\t" + out_array(sorted(unsupported_0))+ "\n")
#    if paf['qstart'] != 0 or paf['qend'] != paf['qlen']:
    report_f.write("unmapped_leading "+ str(paf['qstart']) + "\n")
    report_f.write("unmapped_trailing "+ str(paf['qlen'] - paf['qend']) + "\n")
   
def run_minimap_reference (ref, query,srr_id, workdir):
    paf_file = os.path.join(workdir, srr_id +".paf")
    vcf_file = os.path.join(workdir, srr_id +".vcf")
    minimap2_str = minimap2 + " -c --cs  " + ref + " " + query +" > " +paf_file  + " 2>/dev/null"
    print minimap2_str   
    paftools_str = "sort -k6,6 -k8,8n " + paf_file +" | " + k8 + " " +  fixed_paftools + " call -L20000 -f " + ref + " -  > " + vcf_file
    print paftools_str
    os.system(minimap2_str)
    os.system(paftools_str)

def map_sample(sample_descr, workdir):
    srr_id = sample_descr[0]
    if  not os.path.isdir (workdir):
        os.mkdir(workdir)
    if sample_descr[-2] == "OXFORD_NANOPORE":
#TODO elif ionTORRENT    
        minimap_str = get_minimap_nanopore(srr_id, workdir)
    elif sample_descr[2] == "PAIRED":
        minimap_str = get_minimap_illumina_paired(srr_id, workdir)
    else:
        minimap_str = get_minimap_illumina_single(srr_id, workdir)
    print minimap_str
    sorted_bam = join (workdir, srr_id + ".sorted.bam")
    sam_f =  join(workdir, srr_id + ".sam")
    sort_str = samtools + " sort " + sam_f + " > " + sorted_bam
    depth_file = join (workdir, srr_id + ".depth")
    depth_str = samtools + " depth -a " + join(workdir, srr_id + ".sorted.bam") + " > " + depth_file
    mpileup_file = join(workdir, srr_id + ".vcf")
    mpileup_str = samtools + " mpileup -uf  "  + ref + " " + sorted_bam + " | " + bcftools + " call -mv > " + mpileup_file
    os.system(minimap_str)
    os.system(sort_str)
    os.remove(sam_f)
    os.system(depth_str)
    os.system(mpileup_str)
# ~/other_tools/minimap2/minimap2 -ax sr MN908947.3.fasta /Bmo/dantipov/data/sars_qc/SRR11494470_1.fastq.gz /Bmo/dantipov/data/sars_qc/SRR11494470_2.fastq.gz  -t 20 > aln.sam
# samtools sort SRR11494470.sam > SRR11494470.sorted.bam
# samtools depth -a SRR11494470.sorted.bam SRR11494470.sorted.depth
# samtools mpileup -uf MN908947.3.fasta SRR11494470.sorted.sam | ~/other_tools/bcftools/bcftools call -mv > SRR11494470.vcf

def process_sample(sample_descr, workdir):
    srr_id = sample_descr[0]
    if not os.path.isfile(join(data_pref, srr_id+ "_1.fastq.gz")):
        return
    print sample_descr
    map_sample(sample_descr, workdir)
    query = join(workdir, srr_id + ".fasta")
    run_minimap_reference(ref, query, srr_id, workdir)
    analyze_sample(srr_id, workdir, sample_descr[-2])

def savelist(names, outfile):
    outf = open (outfile, "w")
    for l in names:
        outf.write(l[0] + " " + l[1] + "\n")


def create_download_list(inputfile):
    names = []
    for line in open(inputfile, 'r'):
        arr = line.split(',')
        if len(arr) < 2:
            break
        sra_id = arr[1][2:-1]
        nextstrain_id = arr[2][10:-1]
        names.append([sra_id, nextstrain_id])
    return names 


if __name__ == "__main__": 
#data = StringIO(line) 
#reader = csv.reader(data, delimiter=',')     
    to_download = []
    nextstrain_ids = []
    data_descr = []
    for name in open (big_list, "r"):
        arr = list(csv.reader(StringIO(name)))[0]
#SRR MYSEQ PAIRED ILLUMINA
#        print arr[0] + " " + arr[7] + " " + arr[8] + " " + arr[12]
        tmp = [arr[0], arr[7], arr[8], arr[12]]
        for i in range (0, len(arr)):
            if arr[i].find("/2020") != -1 or arr[i].find("/2019") != -1:
#                print (str(i) +" " + arr[i])
                tmp.append(arr[i])
                data_descr.append(tmp)
                break
#                nextstrain_ids.append([arr[0],arr[i]])
#                to_download.append(arr[0])                               
#                break
#    exit()
#    extract_all_nextstrain(nextstrain_ids)
#                download_sample(arr[0], data_pref)
    
#    Parallel(n_jobs= 10) (delayed(download_sample)(srr, data_pref)
#    for srr in to_download)
#    names = create_download_list (big_list)       
    Parallel(n_jobs=10)(delayed(run_sample)(sample_descr)
    for sample_descr  in data_descr)



#    for name in names:
#        if i == 20:
#            run_sample(name)
#            exit(0)
#        i+= 1
#    process_list(names, data_pref)
#    savelist(names, sys.argv[2])
#    process_sample(sys.argv[1], sys.argv[2])
#process_list(sys.argv[1], sys.argv[2])
