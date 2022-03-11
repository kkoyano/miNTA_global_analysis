#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 05/30/2018

Description: Compare the results of the bed file and the alignment file

Input: a bed file and NTE sam file

Output:
concordinaces of the two files

# results that are the same
# results that are different

Notes:


'''
import argparse
import os
import sys
from subprocess import Popen, PIPE
import re
from itertools import permutations, izip
import gzip
import time
from datetime import datetime

# f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')

class miRNAclass:
    """ A class that holds the information of the annotation, does not have to be a miRNA"""

    def __init__(self, chrom, tstart, tend, readID, Qsize, strand, Qstart, Qend):
        self.chrom = [chrom]
        self.strand = [strand]
        self.tstart = [tstart]
        self.tend = [tend]
        self.readID = readID
        self.Qsize = Qsize
        self.Qstart = [Qstart]
        self.Qend = [Qend]

def get_arguments():
    parser = argparse.ArgumentParser(description='ADD DESCRIPTION')

    parser.add_argument('-s ', dest='sam_file', type=str,
                        help='name of NTE sam file')
    parser.add_argument('-level ', dest='level', type=int,
                        help='level of stringency(1- least, 3-most):\nLevel 1-- remove read that does not match BLAT (keep other reads with modification that correctly mapped\n2-- remove reads that do not match with BLAT and remove modification position, but keep other miRNA reads that match with BLAT for expression 3-- if read was originally mapped with 0 mismatches from bowtie, but does not match BLAT then remove entire miRNA annotation from sample')


    args = parser.parse_args()
    # if len(sys.argv) < 2:
    # parser.print_help()
    # exit()
    fin_sam = args.sam_file
    level = args.level

    return fin_sam,  level


def check_file_exist(file_path):
    if os.path.exists(file_path) == False:
        print 'ERROR!!!'
        print '{} does not exist'.format(file_path)
    return

def compare_coords(blat_obj, chrom, strand, start, end):
    for i in range(len(blat_obj.tstart)):
        if strand == '-':
            strand = '+'
        if blat_obj.chrom[i] == chrom and blat_obj.strand[i] == strand and blat_obj.tstart[i] == start and blat_obj.tend[i] == end:
            return True
    else:
        return False

def level2_level3_mm1(nte_profile, remove_contents_set):
    '''
    determines what to do when there is a read with a modification and the original alignment is wrong.
    for level 2 and 3, we remove this read and any other modifications at this position

    '''
    if nte_profile != '0':
        nte_profile_list = nte_profile.split('|')
        for i in nte_profile_list:
            remove_contents_set.add(i)

    return remove_contents_set

def filter_modifications(nte_profile, remove_contents_set ,nte_nucleotide):
    nte_profile = nte_profile.split('|')

    for i in nte_profile:
        if i in remove_contents_set:
            nte_profile,nte_nucleotide  = '0', '0'
    nte_profile = '|'.join(nte_profile)
    return nte_profile,nte_nucleotide

def output_stringency_level_2(remove_contents_set, nte5, nte3, internal, ncRNA_ID, nte5_profile, nte3_profile, internal_profile):
    if ncRNA_ID in remove_contents_set:
        # if the mirna position was identified in  any other reads, then remove modification
        nte5_profile, nte5   = filter_modifications(nte5_profile, remove_contents_set, nte5)
        nte3_profile, nte3 = filter_modifications(nte3_profile, remove_contents_set, nte3)
        internal_profile, internal = filter_modifications(internal_profile, remove_contents_set, internal)

    return output_line


def apply_BLAT_stringency_filter(fin_sam, level):

    # progress_list = progress(fin_bed)
    # progress_count = 0
    fout= fin_sam+".lvl{}".format(level)
    results_header = ''


    # parse sam file for reads that match the stringency
    # contains either the miRNA name or the modification position
    remove_contents_set = set()

    fin_sam_obj = gzip.open(fin_sam, 'rb') if re.search(r'\.gz$', fin_sam) else open(fin_sam, 'r')


    # opening sam file with original bowtie alignments
    header = ''
    for i in fin_sam_obj:
        if i.startswith("#chr"): #this is the header line
            continue
        i = i.strip().split()
        chrom, strand, start, end, readID, nm_flag, ncRNA_ID, ncRNA_species, nte5_profile, nte3_profile, internal_profile, blat_result = i[0], i[1], i[2], i[3], i[4], i[7], i[10], i[14], i[15], i[16], i[17], i[19]

        # Step1: Go through alignment file and store all of hte miRNA annotations and modification positions that need to be removed from the output results file. This is based on the level of stringency called by the user

        # if the blat result of the alignment is not correct
        if blat_result != 'CORRECT':

            #if the alignment had 0 mismatches
            if nm_flag == "NM:i:0":
                if level == 1:
                    continue # will filter out this read later, because by default we will remove all reads where the alignment does not match the BLAT alignment. So remove all reads but those that are with the "CORRECT flag"
                elif level == 2:
                    continue
                elif level == 3:
                    remove_contents_set.add(ncRNA_ID)
                else: print 'UNKNOWN LEVEL!'

            # the alignment has a mismatch in the read, probably "NM:i:1"
            else:
                if level == 1:
                    continue # will filter out this read later, because by default we will remove all reads where the alignment does not match the BLAT alignment. So remove all reads but those that are with the "CORRECT flag"
                elif level == 2 or level == 3:
                    remove_contents_set = level2_level3_mm1(nte3_profile, remove_contents_set)
                    remove_contents_set = level2_level3_mm1(nte5_profile, remove_contents_set)
                    remove_contents_set = level2_level3_mm1(internal_profile, remove_contents_set)
                else: print 'UNKNOWN LEVEL!'

    fin_sam_obj.close()
    fin_sam_obj = gzip.open(fin_sam, 'rb') if re.search(r'\.gz$', fin_sam) else open(fin_sam, 'r')

    # step2: Go through the blat alignment results file and check if the miRNA or the modification position is in the remove list set. If it is in the set then you can remove the miRNA or the modification
    fout_obj = open(fout, 'w')

    output_file_lines = []
    for i in fin_sam_obj:

        if i.startswith("#chr"): #this is the header line
            fout_obj.write(i)
            continue
        i = i.strip().split()
        readID, nte5, nte3, internal, nm_flag, ncRNA_ID, ncRNA_species, nte5_profile, nte3_profile, internal_profile, blat_result = i[4], i[5], i[6], i[18], i[7], i[10], i[14], i[15], i[16], i[17], i[19]

        readID_count = int(readID.split('|')[1])

        # if the read does not match with the
        if blat_result != 'CORRECT':
            continue

        # if the read matches with blat (CORRECT), then perform output based on stringency requested
        else:
            if level == 2:

                # Check if the modification positions should be removed from analysis
                # if the mirna position was identified in  any other reads, then remove modification
                nte5_profile, nte5   = filter_modifications(nte5_profile, remove_contents_set, nte5)
                nte3_profile, nte3 = filter_modifications(nte3_profile, remove_contents_set, nte3)
                internal_profile, internal = filter_modifications(internal_profile, remove_contents_set, internal)

            if level == 3:
                if ncRNA_ID in remove_contents_set: continue
                nte5_profile, nte5   = filter_modifications(nte5_profile, remove_contents_set, nte5)
                nte3_profile, nte3 = filter_modifications(nte3_profile, remove_contents_set, nte3)
                internal_profile, internal = filter_modifications(internal_profile, remove_contents_set, internal)

            # output line
            i[4], i[5], i[6], i[18], i[7], i[10], i[14], i[15], i[16], i[17], i[19] = readID, nte5, nte3, internal, nm_flag, ncRNA_ID, ncRNA_species, nte5_profile, nte3_profile, internal_profile, blat_result
        output_line = '\t'.join(i)
        output_file_lines.append(output_line)

    all_output = '\n'.join(output_file_lines)
    fout_obj.write(all_output)

    fout_obj.close()

    # some stats we may want to know:
    # 1) how many miRNA that were originaly perfectly mapped, do we lose from this process
    # 2) how many modification of each type do we lose (3' end, 5' end, internal)

    return


def progress(read1_file):
    cmd = ['wc', '-l', read1_file]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    num_lines = int(stdout.strip().split()[0])

    percentages_time = []
    for i in range(0, 10):
        percentages_time.append(num_lines * i / 10)

    return percentages_time


def main():
    start_time = datetime.now()
    fin_sam, level= get_arguments()


    # TESTING #
    # fin_sam = '/Users/kiku/mount/nte_global_analysis_2017-09/unique_reads/data/u87_inhouse/canonical_reads/U87_sRNA_R1_S1_L002_R1_001.NTE.sam.anno_overlap.t0.5.canonical_end'
    # fin_bed = '/Users/kiku/mount/nte_global_analysis_2017-09/unique_reads/data/u87_inhouse/canonical_reads/blat/U87_sRNA_R1_S1_L002_R1_001.NTE.sam.anno_overlap.t0.5.canonical_end.mirna.ALL.psl.best_score.bed'
    # fout = '/Users/kiku/mount/test/sam2bed.txt'
    apply_BLAT_stringency_filter(fin_sam, level)

    #NEXT
    # MODIFY OUTPUT SO THAT IT PRINTS OUT THE QSTART, QEND, AND # MATCHES, # MISMATCHES S
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'
    return


if __name__ == "__main__":
    main()
