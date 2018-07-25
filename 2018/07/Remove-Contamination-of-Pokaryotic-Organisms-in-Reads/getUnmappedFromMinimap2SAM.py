#!/usr/bin/evn python

'''
Purpose: I've mapped PacBio reads to a contamination+insects library using minimap2, and I want to extarct those which aligned to insects or unmapped.

Usage: python xxx.py

'''


import sys
import pysam

# read the insects' genomes
seq_lib = {}
with open('insect_genomes.fa', 'r') as fin:
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            seq_id = line.split()[0].replace('>', '')
            seq_lib[seq_id] = ''


# read the minimap2.sam
clean_ids = {}
with open('minimap2.sam', 'r') as fin:
    for line in fin:
        line = line.strip()
        if line.startswith('@'):
            continue
        qname, flag, rname = line.split('\t')[:3]
        if flag == '4':
            clean_ids[qname] = ''
        else:
            if rname in seq_lib:
                clean_ids[qname] = ''

# print the number of reads left
print len(clean_ids)

# read all raw PacBio data
output = open('minimap2.clean.fa', 'w')
output_bad = open('minimap2.bad.id', 'w')
with open('third_all.fasta', 'r') as fin:
    tmp = 0
    count = 0
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            count += 1
            seq_id = line.split()[0].replace('>', '')
            if seq_id in clean_ids:
                output.write(line + '\n')
                tmp = 1
            else:
                output_bad.write(seq_id + '\n')
                tmp = 0
        else:
            if tmp == 1:
                output.write(line + '\n')
            else:
                continue

print count
