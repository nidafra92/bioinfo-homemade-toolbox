#!/usr/bin/env python3
# 2nd/apr/2018
# by Nicolás D. Franco-Sierra
# get_fq_greater_than.py
# Usage: get_fq_greater_than.py <FASTQ file> <lower read size> <output FASTQ>

from Bio import SeqIO
import sys

input_fastq_file = sys.argv[1]
value = int(sys.argv[2])
output_fastq_file = sys.argv[3]

with open(input_fastq_file) as input_handle, open(output_fastq_file, "a") as output_handle:
    desired_reads = []
    for record in SeqIO.parse(input_handle, "fastq"):
        if len(record.seq) >= value:
            desired_reads.append(record)
    SeqIO.write(desired_reads, output_handle, "fastq")
