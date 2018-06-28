#!/usr/bin/env python

'''
purpose: mask the contaminants with N based on results of Kraken

usage: python xxx.py kraken_res input.fasta > masked.fasta

'''

import sys

# read the result of kraken
kraken_res = {}
contaminant_ratio = {}
with open('bettle.kraken', 'r') as fin:
    for line in fin:
        line = line.strip().split('\t')
        code, seqID, taxon, seq_len, classifies = line[:5]
        contaminant_ratio[seqID] = [int(seq_len)]
        if code == "U":
            kraken_res[seqID] = 0
            contaminant_ratio[seqID].append(0)
            contaminant_ratio[seqID].append(0)
        else:
            tmp_list = []
            gross_kmer_num = 0
            for classify in classifies.split(' '):
                # get all corrdinates (1-based) of contaminants
                taxon, current_kmer_num = classify.split(':')
                current_kmer_num = int(current_kmer_num)
                if taxon == "0":
                    pass
                else:
                    left = gross_kmer_num + 1
                    right = gross_kmer_num + current_kmer_num + 30
                    tmp_tuple = (left, right)
                    tmp_list.append(tmp_tuple)
                gross_kmer_num += current_kmer_num
            # merge overlapped regions
            if len(tmp_list) == 1:
                kraken_res[seqID] = tmp_list
                total_con_len = tmp_list[0][1] - tmp_list[0][0] + 1
                contaminant_ratio[seqID].append(total_con_len)
                contaminant_ratio[seqID].append(float(contaminant_ratio[seqID][1])/contaminant_ratio[seqID][0])
            else:
                contaminants_loci = []
                tmp_min, tmp_max = tmp_list[0]
                i = 1
                while (i < len(tmp_list)):
                    if tmp_list[i][0] < tmp_max:
                        tmp_max = tmp_list[i][1]
                    else:
                        contaminants_loci.append((tmp_min, tmp_max))
                        tmp_min, tmp_max = tmp_list[i]
                    i += 1
                else:
                    contaminants_loci.append((tmp_min, tmp_max))
                kraken_res[seqID] = contaminants_loci
                total_con_len = 0
                for x in contaminants_loci:
                    total_con_len += (x[1] -x[0] + 1)
                contaminant_ratio[seqID].append(total_con_len)
                contaminant_ratio[seqID].append(float(contaminant_ratio[seqID][1])/contaminant_ratio[seqID][0])

out_contaminant_ratio = open('contaminant_ratio.txt', 'w')

for seqID in contaminant_ratio:
   out_contaminant_ratio.write(seqID + '\t' + '\t'.join([str(i) for i in contaminant_ratio[seqID]]) + '\n')


def multi_sub(string, p, c):
    new = {}
    for i, s in enumerate(string):
        new[i] = s
    for i in p:
        new[i] = c
    return ''.join(new.values())


out_masked = open('contaminant_masked.fa', 'w')

# read the fasta of assembly
with open(sys.argv[1], 'r') as fin:	# the fasta
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            if 'seq_seq' in vars():
                if kraken_res[seqID] == 0:
                    out_masked.write('>' + seqID + '\n')
                    out_masked.write(seq_seq + '\n')
                else:
                    # if 'N' > 80%, not print
                    tmp = 0
                    for i in kraken_res[seqID]:
                        tmp += i[1]-i[0]+1
                    # if (tmp/float(len(seq_seq))) > 0.8:
                    #	continue
                    for contaminants_locus in kraken_res[seqID]:
                        seq_seq = multi_sub(seq_seq, range(
                            contaminants_locus[0]-1, contaminants_locus[1]), "N")
                    out_masked.write('>' + seqID + '\n')
                    out_masked.write(seq_seq + '\n')
            else:
                pass
            seqID = line.split()[0].replace('>', '')
            seq_seq = ""
        else:
            seq_seq += line
