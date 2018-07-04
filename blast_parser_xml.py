#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET
import csv

xml_filename = sys.argv[1]
csv_output = sys.argv[2]

tree = ET.parse(xml_filename)
root = tree.getroot()

with open(csv_output, "w") as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=",")
    header = ["Query", "Closest match", "Identity", "E-value", "Bit score"]
    csv_writer.writerow(header)
    for query_id, match_name, e_value, identity, align_len, bit_score in zip(root.iter('{http://www.ncbi.nlm.nih.gov}query-title'), root.iter('{http://www.ncbi.nlm.nih.gov}sciname'), root.iter('{http://www.ncbi.nlm.nih.gov}evalue'), root.iter('{http://www.ncbi.nlm.nih.gov}identity'), root.iter('{http://www.ncbi.nlm.nih.gov}align-len'), root.iter('{http://www.ncbi.nlm.nih.gov}bit-score')):
        identity_percent = int(identity.text)/int(align_len.text) * 100.0
        row = [query_id.text, match_name.text, identity_percent, e_value.text, bit_score.text]
        csv_writer.writerow(row)

# with open(csv_output, "w") as csv_file:
#     csv_writer = csv.writer(csv_file, delimiter=",")
#     header = ["Query", "Closest match", "Identity", "E-value", "Bit score"]
#     csv_writer.writerow(header)
#     for query_id, match_name, e_value, identity, align_len, bit_score in zip(root.iter('Iteration_query-def'), root.iter('Hit_def'), root.iter('Hsp_evalue'), root.iter('Hsp_identity'), root.iter('Hsp_align-len'), root.iter('Hsp_bit-score')):
#         identity_percent = int(identity.text)/int(align_len.text) * 100.0
#         row = [query_id.text, match_name.text, identity_percent, e_value.text, bit_score.text]
#         csv_writer.writerow(row)
