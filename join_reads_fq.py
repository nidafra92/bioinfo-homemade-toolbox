#!/usr/bin/env python3

# written by Nicolás D. Franco-Sierra on 22/nov/2022
# a simple script for joining two readsets from two independent FASTQ files
# usage: join_reads_fq.py <reads_A.fastq> <reads_B.fastq>
# tested on python version 3.7.4
# dependencies: Biopython version 1.74
# outputs: 1_union.fast, 2_intersection.fastq, 3_uniqA.fastq, 4_uniqB.fastq
# stats.txt, containing set stats

import sys
from Bio import SeqIO

filename_A = sys.argv[1]
filename_B = sys.argv[2]

readsA = {}
readsB = {}

with open(filename_A) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        readsA[record.id] = record
        
with open(filename_B) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        readsB[record.id] = record
        

A = set(readsA.keys())
B = set(readsB.keys())


union = readsB
for read in A - B:
    union[read] = readsA.get(read)
    
uniq_A = {read: readsA[read] for read in A - B}
    
uniq_B = {read: readsB[read] for read in B - A}
    
intersect = {read: readsA[read] for read in A & B}
    
SeqIO.write(union.values(), "1_union.fastq", "fastq")
SeqIO.write(intersect.values(), "2_intersection.fastq", "fastq")
SeqIO.write(uniq_A.values(), "3_uniqA.fastq", "fastq")
SeqIO.write(uniq_B.values(), "4_uniqB.fastq", "fastq")

with open("stats.txt", "w") as stats:
    stats.write("{}\t{}\tunion\tintersect\tuniqA\tuniqB\n".format(filename_A, filename_B))
    stats.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(len(A), len(B), len(A | B), len(A & B), len(A - B), len(B - A)))
