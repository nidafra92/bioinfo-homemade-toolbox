#!/usr/bin/env python3
# 3rd/sep/2018
# by Nicolás D. Franco-Sierra
# read_generator.py
# Usage: read_generator.py <fasta full of reads>

import sys
import arrow
import time
from itertools import groupby
import subprocess
import multiprocessing

class AlreadyBaitTested(Exception):
    pass

class NotBaitTested(Exception):
    pass

class NotmtDNA(Exception):
    pass

class AlreadyIdentified(Exception):
    pass

def run_cmd(cmd_str, stdin=None, work_path=None):
    """
    This helper function excecute a bash command as a separate process
    and returns and expection if the command fails.

    It returns the standart output and the standard error of the command as
    str objects.

    """
    #cmd_str = cmd_str.encode('utf-8')
    process = subprocess.Popen(cmd_str, shell = True, stdout = subprocess.PIPE,
            stdin = subprocess.PIPE, stderr = subprocess.STDOUT, cwd=work_path,
            encoding='utf8')
    stdout_str, stderr_str = process.communicate(stdin)
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str,
                            process.returncode))
    return stdout_str, stderr_str


def fasta_reader(fasta_filename):
    """
    Helper function: opens a FASTA file and yields 'read tuples' from FASTA
    entries contained in the original file.
    Generator fuction
    """
    with open(fasta_filename) as fasta_handle:
        fasta_iter = (x[1] for x in groupby(fasta_handle, lambda line: line[0] == ">"))
        for header in fasta_iter:
            headerStr = header.__next__()[1:].strip()
            seq="".join(s.strip() for s in fasta_iter.__next__())
            yield (headerStr, seq)

def wait_until(specified_dt: arrow.Arrow):
    """
    Stay in a loop until the specified date and time.
    """
    # Initially check every 0.1 seconds.
    refresh = 0.1
    current_dt = arrow.utcnow()

    while current_dt < specified_dt:
    # Check every millisecond if close to the specified time.
        current_dt = arrow.utcnow()
        if (specified_dt - current_dt).seconds < 0.11:
            refresh = .001

    time.sleep(refresh)

class ONTread:
    """
    This is a basic ONTread class. It receives a 'read tuple' ("seq_id", "ATCG")
    and creates a read object extracting relevant information from the
    description string (>fasta_header)
    """
    @staticmethod
    def parse_header(header_string):
        header_list = header_string.split()
        read_id = header_list[0]
        run_id = header_list[1].split("=")[-1]
        read_number = header_list[2].split("=")[-1]
        channel_number = header_list[3].split("=")[-1]
        start_time = arrow.get(header_list[4].split("=")[-1])
        return (read_id, run_id, read_number, channel_number, start_time)
    def __init__(self, read_tuple):
        #print(read_tuple)
        raw_header = read_tuple[0]
        raw_sequence = read_tuple[1]
        header_tuple = ONTread.parse_header(raw_header)
        self.raw_header = raw_header
        self.readname = header_tuple[0]
        self.runid = header_tuple[1]
        self.readnum = header_tuple[2]
        self.chnum = header_tuple[3]
        self.time = header_tuple[4]
        self.seq = raw_sequence
        self.is_mito = "not_tested"
        self.IDresults = "not_tested"
    def set_simtime(self, time_delta):
        self.simtime = self.time + time_delta
    def __len__(self):
        return len(self.seq)
    def __str__(self):
        return "{}\t{}\t{}\t{} bp".format(self.readname, self.runid, self.simtime, len(self.seq))
    def use_as_bait(self, blastjob_obj):
        if self.is_mito == "not_tested":
            blastable = ">{}\n{}".format(self.readname, self.seq)
            results, error = run_cmd(str(blastjob_obj), stdin=blastable)
            if len(results) != 0:
                self.is_mito = True
                print("The read is mitochondrial.")
            else:
                self.is_mito = False
                print("This read is NOT mitochondrial.")
        else:
            raise AlreadyBaitTested("This read has been already tested for mtDNA homology. Cannot be tested again.")
    def blast_for_ID(self, blastjob_obj):
        if self.IDresults == "not_tested":
            if self.is_mito == True:
                blastable = ">{}\n{}".format(self.readname, self.seq)
                results, error = run_cmd(str(blastjob_obj), stdin=blastable)
                self.IDresults = results
                return results
                print("I have to implement it...")
            elif self.is_mito == False:
                raise NotmtDNA("This read not mitochondrial, so it is not going to be blasted for identification...")
            else:
                raise NotBaitTested("This read needs to be tested first for mtDNA homology before to be identified.")
        else:
            raise AlreadyIdentified("This read has been already identified. Cannot be blasted again.")

class BlastJob:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
    def __str__(self):
        attr_dict = vars(self)
        command_line = attr_dict.get("task")
        for attr, value in attr_dict.items():
            if attr != "task":
                command_line += " -{} {}".format(attr, value)
        return command_line

def library_loader(sorted_readlist, time_delta, minion):
    for read in sorted_readlist:
        #print("loading {}".format(read.readname))
        read.set_simtime(time_delta)
        #print("new time set")
        minion.put(read)

def sequencing_task(minion, mt_bait_queue):
    while not minion.empty():
        sequenced_read = minion.get(True)
        print("waiting...")
        wait_until(sequenced_read.simtime)
        print("A new read has been sequenced at ".format(sequenced_read.simtime.format("HH:mm:ss ZZ")))
        mt_bait_queue.put(sequenced_read)
        print(sequenced_read)

def mtbait_tester(mt_bait_queue, spp_id_queue, mt_bait_job):
    while not mt_bait_queue.empty():
        sequenced_read = mt_bait_queue.get(True)
        sequenced_read.use_as_bait(mt_bait_job)
        if sequenced_read.is_mito:
            spp_id_queue.put(sequenced_read)

def spp_id_blaster(spp_id_queue, spp_id_job):
    while not spp_id_queue.empty():
        mtdna_read = spp_id_queue.get()
        mtdna_read.blast_for_ID(spp_id_job)
        print("Identified as: {}".format(mtdna_read.IDresults))


if __name__ == "__main__":
    original_fasta = sys.argv[1]
    #define blast jobs
    mt_bait_job = BlastJob(task="blastn", max_target_seqs="1", word_size="7", gapopen="2", gapextend="2", penalty="-3", reward="2", max_hsps="1", perc_identity="60", db="mtbait_db/mtDNA_bait_Metazoa.fasta", outfmt='6')
    #
    spp_id_job = BlastJob(task="blastn", max_target_seqs="1", word_size="11", gapopen="2", gapextend="2", penalty="-3", reward="2", max_hsps="1", perc_identity="85", db="mtbait_db/mtDNA_bait_Metazoa.fasta", outfmt='"6 qaccver saccver pident ppos length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore staxid sciname scomname sskingdom stitle"')
    unsorted_reads = []
    for read_tuple in fasta_reader(original_fasta):
        unsorted_reads.append(ONTread(read_tuple))
    print("Reading reads from {}".format(original_fasta))
    print("Total of {} ONT reads were parsed!".format(len(unsorted_reads)))

    print("Sorting reads by start time.")
    sorted_reads = sorted(unsorted_reads, key = lambda read: read.time)

    original_start_time = sorted_reads[0].time
    # defining start time for simulation
    sim_start_time = arrow.utcnow().shift(minutes=+0.1)
    time_delta = sim_start_time - original_start_time
    print("Initializing simulation start time!")
    print("Start time (UTC): {}".format(sim_start_time.format('YYYY-MM-DD HH:mm:ss ZZ')))
    print(sim_start_time.humanize())

    #virtual_minion = multiprocessing.Queue()
    sim_manager = multiprocessing.Manager()
    virtual_minion = sim_manager.Queue()
    mtDNA_test = sim_manager.Queue()
    spp_ID_test = sim_manager.Queue()

    load_library = multiprocessing.Process(target=library_loader, args=(sorted_reads, time_delta, virtual_minion))
    sequencing_run = multiprocessing.Process(target=sequencing_task, args=(virtual_minion, mtDNA_test))
    mt_bait_checker = multiprocessing.Process(target=mtbait_tester, args=(mtDNA_test, spp_ID_test, mt_bait_job))

    enhanced_mtBlaster = multiprocessing.Process(target=spp_id_blaster, args=(spp_ID_test, spp_id_job))

    processes = [load_library, sequencing_run, mt_bait_checker, enhanced_mtBlaster]

    print("Readjusting 'birth time' for each read")

    #for process in processes:
    #    process.start()


    load_library.start()
    #load_library.join()
    sequencing_run.start()
    mt_bait_checker.start()
    #virtual_minion.close()
    #sequencing_run.join()
    #enhanced_mtBlaster.join()
    #mt_bait_checker.join()
    #print("Library successfuly loaded.")
    #print("Starting sequencing")


#    for read in sorted_reads:
#     read.set_simtime(time_delta)
# print("Completed!")
#
#
# print("Testing if time works ok...")
#
# for read in sorted_reads:
#     wait_until(read.simtime)
#     print("A new read has been sequenced at ".format(read.simtime.format("HH:mm:ss ZZ")))
#     print(read)
#
#
# mt_bait_job = BlastJob(task="blastn", max_target_seqs="1", word_size="7", gapopen="2", gapextend="2", penalty="-3", reward="2", max_hsps="1", perc_identity="60", db="mtbait_db/mtDNA_bait_Metazoa.fasta", outfmt='6')
#
# spp_id_job = BlastJob(task="blastn", max_target_seqs="1", word_size="11", gapopen="2", gapextend="2", penalty="-3", reward="2", max_hsps="1", perc_identity="85", db="mtbait_db/mtDNA_bait_Metazoa.fasta", outfmt='"6 qaccver saccver pident ppos length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore staxid sciname scomname sskingdom stitle"')
#
# print(read.is_mito)
# read.use_as_bait(mt_bait_job)
# print(read.is_mito)
# read.blast_for_ID(mt_bait_job)
# print(read.IDresults)
# print(read.is_mito)
# read.use_as_bait(mt_bait_job)
# read.blast_for_ID(mt_bait_job)
