#!/usr/bin/env python3
# 31th/ene/2018
# by Nicolás D. Franco-Sierra
# single_fasta_from_multifasta.py
# Usage: single_fasta_from_multifasta.py <FASTA with multiple records> <Output file>

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

multifasta_file = sys.argv[1]
concatfasta_file = sys.argv[2]

with open(multifasta_file) as input_handle, open(concatfasta_file, "a") as output_handle:
    con_seq = ""
    for record in SeqIO.parse(input_handle, "fasta"):
        con_seq += record.seq
    new_record = SeqRecord(con_seq,'concat_seq', '', '')
    SeqIO.write(new_record, output_handle, "fasta")
