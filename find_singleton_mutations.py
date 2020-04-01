import sys
import os
import random
import subprocess
import shutil
import re
from datetime import date
from genericpath import isdir, exists
from os.path import join
from joblib import Parallel, delayed
#import check_mash

#MN908947.3      29706   .       G       T       60      .       QNAME=hCoV-19_Iceland_232_2020_EPI_ISL_417560_2020-03-17;QSTART=29706;QSTRAND=+ GT      1/1
ref_id = "MN908947.3"


def get_date(str_date):
    arr = str_date.split('-')
#some have only year, that's bad
    year = int(arr[0])
    if len(arr) > 1:
        month = int(arr[1])
    else:
        month = 1
    if len(arr) > 2:
        day = int(arr[2])
    else:
        day = 1
    return date(year, month, day)

#for pos sorting
def compare(a, b):
    pos_a = int(a.split()[1])
    pos_b = int(b.split()[1])
    if pos_a < pos_b:
        return -1
    elif pos_a > pos_b:
        return 1
    else:
        if a < b:
            return -1
        elif a> b:
            return 1
        else:
             return 0

def parse(infile, before):
    table = {}
    limit = get_date(before)
    for line in open (infile, 'r'):
        arr = line.split()
        if arr[0] == ref_id:
            id = arr[1] + "_" + arr[3] + "_" + arr[4]
            if  not (id in table.keys()):
                table[id] = []
            table[id].append(line.strip())
    unique = []
    for i in table.keys():
        if len(table[i]) == 1:
            unique.append (table[i][0])
    old = []
    for u in unique:
        info = u.split()[7]
        current = info.split(";")[0].split('_')[-1]
        current_date = get_date(current)
        if current_date < limit:
           old.append(u)
    from functools import cmp_to_key
    old_sorted = sorted(old, key=cmp_to_key(compare))
    for u in old_sorted:
        print (u)



if len(sys.argv) != 3:
    print ("Usage: "+ sys.argv[0] + " snps.vcf date \n, date should be in YYYY-MM-DD format. Returns all unique snps happened before given date")
infile = sys.argv[1]
before = sys.argv[2]
time = 0
parse(infile, before)
