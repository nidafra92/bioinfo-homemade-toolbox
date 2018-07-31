#!/usr/bin/env python3
# 27th/jun/2018
# by Nicolás D. Franco-Sierra
# ONT_read-timesorter.py
# Usage: get_fast5_from_fasta.py <fasta file unsorted> <fasta file sorted> destination>


import sys, os
from optparse import OptionParser
from itertools import groupby
import dateutil.parser

import subprocess
import time
from datetime import timedelta
from multiprocessing.dummy import Pool as ThreadPool
import math

start_time = time.monotonic() # I want to record execution time.

parser = OptionParser(version="%prog {}".format("1"),
            usage = "not written yet")
parser.add_option("-i", "--input_fasta", dest="input_fasta",
                    metavar="FILE",
                    type="string",
                    help="Name and path of unsorted ONT reads in fasta format. ")
parser.add_option("-o", "--output_fasta", dest="output_fasta",
                    metavar="FILE",
                    type="string",
                    default="sorted.fasta",
                    help="Name and path of where the sorted fasta should \
                    end up.")
parser.add_option("-r", "--result_file", dest="result_file",
                    metavar="FILE",
                    type="string",
                    default="sorting_stats.tab",
                    help="Name and path of where the stats file should end up. ")
parser.add_option("-f", "--filter_list", dest="filter_list",
                    metavar="FILE",
                    type="string",
                    default="no_file",
                    help="Blast 6 fmt of barcode positives.")
# parser.add_option("-t", "--temp_directory", dest="temp_directory",
#                     metavar="PATH",
#                     type="string",
#                     default="tmp",
#                     help="Name and path of temporary directory where \
#                         calculations should take place. ")
parser.add_option("-j", "--jobs", dest="jobs",
                    default=1,
                    metavar="jobs",
                    type="int",
                    help="Specifies the number of jobs (operations) to run \
                    in parallel.")

(options, remaining_args) = parser.parse_args()

jobs = options.jobs

if not options.input_fasta:
    parser.error("\n\n\tMissing parameter --input_fasta FILE\n\n")


#get argvs from flags
original_fasta = options.input_fasta
sorted_fasta_file = options.output_fasta
stats_result_file = options.result_file
filter_list = options.filter_list

def fasta_reader(fasta_filename):
    with open(fasta_filename) as fasta_handle:
        fasta_iter = (x[1] for x in groupby(fasta_handle, lambda line: line[0] == ">"))
        for header in fasta_iter:
            headerStr = header.__next__()[1:].strip()
            seq="".join(s.strip() for s in fasta_iter.__next__())
            yield (headerStr, seq)

import dateutil.parser

def parse_header(header_string):
    header_list = header_string.split()
    read_id = header_list[0]
    run_id = header_list[1].split("=")[-1]
    read_number = header_list[2].split("=")[-1]
    channel_number = header_list[3].split("=")[-1]
    start_time = dateutil.parser.parse(header_list[4].split("=")[-1])
    return (read_id, run_id, read_number, channel_number, start_time)


class ONTread:
    """docstring forONTreadf __init__(self, arg):
        superONTread__init__()
        self.arg = arg
    """
    def __init__(self, read_tuple):
        header_tuple = parse_header(read_tuple[0])
        self.readname = header_tuple[0]
        self.runid = header_tuple[1]
        self.readnum = header_tuple[2]
        self.chnum = header_tuple[3]
        self.time = header_tuple[4]
        self.seq = read_tuple[1]
        #self.keep = True
    def __str__(self):
        return "{}\t{}\t{}\t{} bp".format(self.readname, self.runid, self.time, len(self.seq))
    #def activate(self):
        #self.keep = False

unsorted_reads = []

for read_tuple in fasta_reader(original_fasta):
    unsorted_reads.append(ONTread(read_tuple))

print("Reading reads from {}".format(original_fasta))
print("Total of {} ONT reads were parsed!".format(len(unsorted_reads)))

# print("\n\nUNSORTED\n\n")
# for formatted_read in unsorted_reads:
#     print(formatted_read)
print("Sorting reads by start time.")

sorted_reads = sorted(unsorted_reads, key = lambda read: read.time)

print("Completed!")

if filter_list != "no_file":
    print("Filter list with read IDs detected!")
    filter_set = set([])
    reads_to_keep = []
    with open(filter_list) as ID_list:
        for line in ID_list:
            read_id = line.strip().split("\t")[0]
            filter_set.add(read_id)
    for read in sorted_reads:
        if read.readname in filter_set:
            reads_to_keep.append(read)
    print("Total of {} reads were kept.".format(len(reads_to_keep)))
    sorted_reads = reads_to_keep



with open(stats_result_file,"a") as stats:
#print("\n\nSORTED\n\n")
    stats.write("read_pos\tread_id\trunid\tstart_time\tread_length\n")
    read_count = 0
    for formatted_read in sorted_reads:
        read_count += 1
        stats.write("{}\t{}\n".format(read_count,formatted_read))

print("A total of {} reads were processed and sorted.".format(read_count))

print("Writing stats file in: {}".format(stats_result_file))

with open(sorted_fasta_file, "a") as out:
    for read in sorted_reads:
        out.write(">{}\n{}\n".format(read.readname,read.seq))

print("Writing sorted FASTA file in: {}".format(sorted_fasta_file))
