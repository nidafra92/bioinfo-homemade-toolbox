#!/usr/bin/env python3
# 2nd/apr/2018
# by Nicolás D. Franco-Sierra
# get_fast5_from_fastq.py
# Usage: get_fast5_from_fastq.py <FASTQ file> <input fast5 folder> <desired destination>

from Bio import SeqIO
import sys
import os
import shutil

input_fastq_file = sys.argv[1]
input_fast5_folder = sys.argv[2]
destination = sys.argv[3]

def get_data_from_fastq(fastq_filename):
    read_channel_set = set([])
    with open(fastq_filename) as fastq_file:
        for line in fastq_file.read().split("\n")[::4]:
            cur_line = line.split()
            print(cur_line)
            if len(cur_line) != 0:
                channel = cur_line[-2][3:]
                read = cur_line[-3][5:]
                search_string_pairs = (read, channel)
                read_channel_set.add(search_string_pairs)
    return read_channel_set

filenames_to_copy = get_data_from_fastq(input_fastq_file)

for root, _, filenames in os.walk(input_fast5_folder):
    for filename in filenames:
        read = filename.partition("_read_")[2].partition("_ch_")[0]
        channel = filename.partition("_ch_")[2].partition("_strand")[0]
        search_tuple = (read, channel)
        if search_tuple in filenames_to_copy:
            shutil.copy(os.path.join(root, filename), destination)
