#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 08/30/2018

Description: pipeline takes in bowtie parameters and fasta or fastq files and performs initial mapping and 2 s of 3' trimming, and 2 s of 5' trimming

Usage:
python /u/nobackup/gxxiao3/kkoyano/scripts/trimming_pipeline/NTE_trimming_pipeline.V4.py -i ./test/test.fa -b mapping_opt_unique.txt -x /u/nobackup/gxxiao2/apps/genomes/hg19/bowtieindex/all.chr.fa --fasta

Input:

Output:

Notes:
dependencies:
(must be in current directory, or an alias linking to the directory it is in)
built on python 2.7.2
bowtie version 1.1.2
samtools faidx
helper_function.py in same directory

---------------------------------
Updates
---------------------------------
V4:
1. Fix bug where if a di-3'NTE is seen as NOT a 3'NTE, then the script will still look for modifications on the 5' end. If a 5' modification exists, then add the 3'NTE as a internal mismatch and record 5'NTE

 - case1:                                     G------------TT
 - after checking 3' NTE:                     G------------T-
        - KEEP because 5' NTE
 - case2:                                     ----X--------TT
 - after checking 3' NTE:                     ----X--------T-
        - REMOVE because 2 internal mismatches, too many


2. Add a field that records which step the read was aligned/processed
    - 3nte-trim1
    - 3nte-trim2
    - 5nte-trim1
    - 5nte-trim2
    - both-trim1 [ NEW!! ]*

3. both-trim1 [ NEW!! ]*
        - if modification does not pass the last 5' end trim, then try trimming 1 nucletoide on both ends, this will find cases where there is an internal modifciation and a 5' NTE and 3' NTE i.e:  G--------X---------T

4. could try enhancing the speed by implementing bedtools to look up the sequences instead of GenomeFetch.py
do entire lookups in batches rather than one at a time, so for each round instead of each alignment
    bedtools getfasta -fi /u/nobackup/gxxiao2/apps/genomes/hg19/all.chr.fa -bed test.bed




'''
import argparse
import os
import sys
from subprocess import Popen, PIPE, call
import re
from itertools import permutations
from helper_functions import *
import cPickle as pickle
import time
from datetime import datetime


MISMATCH_NUMBER = 1

def get_arguments():

    parser = argparse.ArgumentParser(description="Trimming pipeline, takes in either fasta, fastq, or a sam file of alignments, and outputs a modified sam file that contains information about 3' and 5' NTE, it also will output a summary file with the breakdown ")
    #manditory arguments
    parser.add_argument('-i', dest='input_file', required=True, type=str,
                help='name and path of input sequencing file, needs to be either a fasta, fastq, or sam format')
    parser.add_argument('-x', dest='idx_path', required=True, type=str,
                help='path of the index file and prefix for bowtie alignment\nex: /u/nobackup/gxxiao2/apps/genomes/hg19/bowtieindex/all.chr.fa')
    #optional arguments
    parser.add_argument('-o', '--output', dest='output_dir', type=str,
                help='name of output directory')
    parser.add_argument('-p', '--path',dest='path', type=str,
                help='path of directory that contains input file')

    #different input types for sequences
    parser.add_argument('-f', '--fasta', dest='fasta_file', action='store_true',
                help='name of fasta file')
    parser.add_argument('-q','--fastq', dest='fastq_file', action='store_true',
                help='name of fastq file')
    parser.add_argument('-s', '--sam', dest='sam', action='store_true',
                help='input file is in SAM format')

    #bowtie mapping options
    parser.add_argument('-b', '--bowtie', dest='bowtie', type=str,
                help='bowtie mapping arguments in a file; see bowtie manual for details')

    parser.add_argument('-T', '--trim3', dest='trim3', action='store_true',
                help='flag that they want 3 prime trimming ')
    parser.add_argument('-F', '--trim5', dest='trim5', action='store_true',
                help='flag that they want 5 prime trimming')
    parser.add_argument('-A', '--trim3and5', dest='trim3and5', action='store_true',
                help='flag that they want 3 prime and 5 prime trimming')
    # parser.add_argument('-e ', dest='exp', type=str,
                    # help='a directory containing the experimental files from my predictions ')

    args = parser.parse_args()

    fin = args.input_file
    fout = args.output_dir
    directory = args.path
    fa = args.fasta_file
    fq = args.fastq_file
    sam = args.sam
    trim3 = args.trim3
    trim5 = args.trim5
    trim3and5 = args.trim3and5
    mapping_file=args.bowtie
    input_err_msg = 'INPUT ERROR: Program Exiting\n\nMust choose one of the three options\n1)Fasta file\n2)fastq file\n3)SAM format alignment file'
    if sum([fa,fq,sam]) > 1:
        print input_err_msg
        sys.exit()
    if sum([fa,fq,sam]) == 0:
        fa = True
    if sum([trim3,trim5,trim3and5]) > 1:
        print input_err_msg
        sys.exit()
    if sum([trim3,trim5,trim3and5]) == 0:
        trim3and5 = True

    input_options= [fa,fq,sam]
    analysis_options= [trim3,trim5,trim3and5]
    input_dict = {0: '-f', 1: '-q', 2: 'sam'}
    analysis_dict = {0: 'trim3', 1: 'trim5', 2: 'trim3and5'}
    seq_file_type = False
    for i in range(len(input_options)):
        if input_options[i] == True:
            seq_file_type = input_dict[i]
            break
    for i in range(len(analysis_options)):
        if analysis_options[i] == True:
            analysis_type = analysis_dict[i]
            break

    if fout == None:
        fout = fin+'.sam'

    #seq_file_type will hold one of the values: '-f, -q, or sam'
    return fin, fout, directory, seq_file_type, analysis_type, args.idx_path, mapping_file


class NTE_class:
    """ A class that holds the information of the annotation, does not have to be a miRNA"""
    def __init__(self, sequence):
        self.original_RNA_sequence = sequence
        self.nte3_trim = ''
        self.nte5_trim = ''


def bowtie_mapping(sequence_file, fout, mapping_file, idx_path, seq_file_type, rnd=1, trim_description='init', trim_round_number = 0):
    # makes bowtie mapping commands
    # file line contains:
    # bowtie -v 1 -k 100 --best --strata -S
    # full command: bowtie -v 1 -k 100 --best --strata -S [--trim3 1/--trim5 1] --un unmapped.fa -f ~/index_path/prefix sample.fa alignment.sam
    # idx_path = /u/nobackup/gxxiao2/apps/genomes/hg19/bowtieindex/all.chr.fa
    seq_file_ext = ''
    if seq_file_type == "-f":
        seq_file_ext ='fa'
    if seq_file_type == "-q":
        seq_file_ext = 'fq'

    with open(mapping_file, 'r') as mapping:
        for line in mapping:
            if line.startswith('#'):
                continue
            else:
                cmd = line.strip().split(' ')
                if '-S' not in cmd:
                    cmd.append('-S')
                if '--un' not in cmd:
                    cmd.append('--un')
                    unmapped = '{}-unmapped-{}-{}.{}'.format(fout, trim_description, str(rnd), seq_file_ext)
                    cmd.append(unmapped)

                # adding trim3 or trim5 to bowtie command
                if trim_description == 'trim3':
                    cmd.append('--trim3')
                    cmd.append(str(trim_round_number))
                elif trim_description == 'trim5':
                    cmd.append('--trim5')
                    cmd.append(str(trim_round_number))
                elif trim_description == 'both':
                    cmd.append('--trim5')
                    cmd.append(str(trim_round_number))
                    cmd.append('--trim3')
                    cmd.append(str(trim_round_number))
                else:
                    pass

                unmapped_fn_idx = cmd.index('--un')
                unmapped_fn = cmd[unmapped_fn_idx+1]

                cmd.append(idx_path)
                cmd.append(seq_file_type)
                cmd.append(sequence_file)


                print 'running command...\t{}'.format(' '.join(cmd))

                mapping_results_list, stderr = run_command(cmd)
                mapping_results_list = mapping_results_list.split('\n')

    # removing reads that are unmapped
    cleaned_mapping = []
    for i in range(len(mapping_results_list)):
        line = mapping_results_list[i].strip('\n').split()
        if line == []: continue
        if line[1] != '4':
            trim_name_description = trim_description +'-{}'.format(trim_round_number)
            line.append(trim_name_description)
            line = '\t'.join(line)
            line = line +'\n'
            cleaned_mapping.append(line)

    return unmapped_fn, cleaned_mapping


def trimming_pipeline(input_file, fout, trim_dict, seq_file_type, trim_end = 'three', trim_round_number = 1 ):
    # treating trimming for FASTA formatted files
    # input_file = FILE; regular fasta file, contains full sequence that will be trimmed depending on the "trim end" flag
    # fout = STRING; the output prefix
    # trim_dict = DICT obj.
    # seq_file_type = STRING <-f, -q, sam>
    # trim_end = STRING <'three', 'five'>
    # rnd = INT round of trimming

    if seq_file_type == 'sam': #mapping file #only perform the 'sam' type the first time, then after, you will just create a fasta file to feed into this function
        trim_dict = {}
        trimmed_reads_list = []
        mapped_reads_list = []
        unmapped_fn = []

        with open(input_file, 'r') as fin:
            for line in fin:
                if line.startswith('@'):
                    continue
                else:
                    line = line.strip().split()


                    if line[1] == '4': #unmapped, needs to get trimmed
                        ID, seq = '>'+line[0], line[9] # readID, seq
                        unmapped_string =  ID + "\n" + seq + '\n'
                        unmapped_fn.append(unmapped_string) # original unmapped sequences
                        # trimming the 3'-end (end of sequence) one nucleotide.
                        # key = ID and value is trimmed nucleotide
                        if trim_end == 'three':
                            trim_dict[ID] = line[9][-1:]
                            fa_output = ID + "\n" + seq[:-1] + '\n'
                            trimmed_reads_list.append(fa_output) # list of all new trimmed reads
                        # trimming the 5'-end (beginning of sequence) one nucleotide.
                        if trim_end == 'five':
                            trim_dict[ID] = line[9][:1]
                            fa_output = ID + "\n" + seq[1:] + '\n'
                            trimmed_reads_list.append(fa_output)

                    else: #line contains read that was mapped
                        mapped_reads_list.append('\t'.join(line))


        with open('{}-trim{}.fa'.format(fout, trim_round_number), 'w') as trim_fa:
            trimmed_read_sequences = '\n'.join(trimmed_reads_list)
            trim_fa.write(trimmed_read_sequences)
            # for i in trimmed_reads_list:
                # trim_fa.write(i)

        with open('{}-unmapped1.fa'.format(fout), 'w') as un_fn:
            unmapped_fn = '\n'.join(unmapped_fn)
            un_fn.write(unmapped_fn)
            # for i in unmapped_fn:
                # un_fn.write(i)

        return trim_dict, mapped_reads_list

    # the input sequence file is in a FASTA format
    if seq_file_type == '-f':

        with open(input_file, 'r') as input_file_obj:
            for i in input_file_obj:
                i = i.strip('\n')
                if i.startswith('>'):
                    readID = i[1:]
                    sequence = input_file_obj.next()
                    sequence = sequence.strip()
                    if readID in trim_dict:
                        nte_class_obj = trim_dict[readID]
                        if trim_end == 'three':
                            nte_class_obj.nte3_trim = sequence[-trim_round_number:]
                        elif trim_end == 'five':
                            nte_class_obj.nte5_trim = sequence[:trim_round_number]
                        elif trim_end == 'both':
                            nte_class_obj.nte3_trim = sequence[-trim_round_number:]
                            nte_class_obj.nte5_trim = sequence[:trim_round_number]
                        else:
                            print 'unidentified commdand... error on line 297'
                        trim_dict[readID] = nte_class_obj

                    else:
                        nte_class_obj = NTE_class(sequence)
                        if trim_end == 'three':
                            nte_class_obj.nte3_trim = sequence[-trim_round_number:]
                        elif trim_end == 'five':
                            nte_class_obj.nte5_trim = sequence[:trim_round_number]
                        elif trim_end == 'both':
                            nte_class_obj.nte3_trim = sequence[-trim_round_number:]
                            nte_class_obj.nte5_trim = sequence[:trim_round_number]
                        else:
                            print 'unidentified commdand... error on line 310'
                        trim_dict[readID] = nte_class_obj

        return trim_dict


def check_sequence(chrom, start, seq, genome_fa ='/u/nobackup/gxxiao2/apps/genomes/hg19/all.chr.fa'):
    # def check_sequence(chrom, start, seq, genome_fa ='/u/home/k/kkoyano/nobackup-gxxiao3/apps_and_genome_indexes/genomes/rn5/rn5.all.chr.fa'):
    # Note: samtools faidx will ALWAYS output the result with respect to the postive strand
    end = int(start) + len(seq)- 1
    cmd = 'samtools faidx {} {}:{}-{}'.format(genome_fa, chrom, start, end)
    stdout, stderr = run_command(cmd)
    seq = stdout.split('\n')[1]
    seq = seq.upper()
    return seq

def check_sequence_genome_pkl(chrom, start, end, genome):
    '''
    checking the sequencing using stephens pickled genome
    '''
    start=int(start)
    ref_seq = genome[chrom][start-1:end-1].upper()

    if len(ref_seq) is '':
        print 'Not a valid alignment, {}:{}:{}'.format(chrom, start-1, end-1)
    return ref_seq


def check_cigar_3nte(cigar, chrom, start, strand, seq, trim_dict, readID, genome):
    '''
    note: input seq is with respect to the positive strand

    Function has 3 major tasks:

    1. checks cigar string for 3'NTE modifications from the alignment
    2. checks 3'NTE from trimmed dictionary, if so then add to nte3 found in cigar
    3. checks if the trimmed NTEs are true modifications ( i.e if there is more than 1 3'NTE, the mismatch must be consecutive)
    4. If 3'end modifications exists from the trimemd dictionary, it adds these modifications to the cigar string.
    '''

    nte3 = ''
    cigar_string = cigar.split(':')[2]  # ['MD', 'Z', '14C8'] --> '14C8'
    matches_count_list = re.split(r'[ACTG]', cigar_string)  # ['14','8'] or ['21','0', '1'], of no mismatches will be 1 number
    add_nte3 = ''
    original_sequence = seq

    # 1. check cigar for 3'NTE
    if len(matches_count_list) != 1:
        for i in range(len(matches_count_list)):
            # any 3'NTE for reads on the positive strand, will appear at the end of the cigar string, and appear as a 0 position
            # ex: ['22','0']

            # Note: This assumes that there is only one mismatch allowed in the current alignment and thus, there is only one possible 3'-end modification (this will be the only allowed mismatch in alignment)
            # print matches_count_list
            # print matches_count_list[-1 - i]

            if strand == '0':
                if matches_count_list[-1 - i] == '0':
                    nte3_base = seq[-1 - i]
                    nte3 = nte3_base + nte3
                else:
                    break
            # any 3'NTE for reads on the negative strand, will appear at the beginning of the cigar string
            else:
                if matches_count_list[i] == '0':
                    nte3_base = seq[0 + i]
                    nte3 = nte3 + nte3_base  # this is w.r.t the positive strand
                    nte3 = rev_complement(nte3)
                else:
                    break

    # 2. checks 3'NTE from trimmed dictionary
    if readID in trim_dict:
        nte_class_obj = trim_dict[readID]
        add_nte3 = nte_class_obj.nte3_trim
        original_sequence = nte_class_obj.original_RNA_sequence

    nte3 =  nte3 + add_nte3 # this is currently w.r.t positive strand
    nte_len = len(nte3)

    # print 'nte3 before checking- {}'.format(nte3)
    # initializing variables
    count = 0
    cigar_ext = ''
    new_nte = ''
    nte_bool = True
    new_cigar = ''
    updated_sequence = ''
    add_nucleotides_to_ref = False

    if nte_len is 0:
        return new_nte, cigar, start, seq

    # 3. checks if the trimmed NTEs are true modifications ( i.e if there is more than 1 3'NTE, the mismatch must be consecutive)
    if strand == '0':
        # checking to see of the sequence is correct
        end = int(start) + len(seq) + len(add_nte3)
        refseq = check_sequence_genome_pkl(chrom, start, end, genome)
        if refseq is '':
            print 'alignment is out of bounds, skipping'
            return 0, 0, 0, 0
        # print refseq
        refseq = refseq[-nte_len:]  # get the last letters of reference sequence that correspond to 3'NTE
        # 4. If 3'end modifications exists from the trimemd dictionary, it adds these modifications to the cigar string.
        # print 'ref:{} nte:{}'.format(refseq, nte3)
        for i in range(0, nte_len):
            # the NTE nucleotide matches the reference sequence, so it is no longer considered a NTE
            if refseq[-i - 1] == nte3[-i - 1]:
                count +=1
                nte_bool = False
            # ref and nte do not match
            else:
                cigar_ext = refseq[-i - 1] + str(count) + cigar_ext
                if nte_bool is True:
                    new_nte = nte3[-i - 1] + new_nte  # adding NTE from L <-- R i.e: T  --> AT --> AAT
                count = 0
        # rewritting cigar string, if 3'-end modification initially existed, remove and rewrite.
        if matches_count_list[-1] is '0' :
            cigar_string = cigar_string[:-2]
        new_cigar = cigar_string + cigar_ext
        updated_sequence = seq + add_nte3

    # negative strand, so the 3'-NTE is at the beginning of SAM sequence and cigar.
    else:
        start = int(start) - len(add_nte3)

        end = int(start) + len(seq) + len(add_nte3)
        refseq = check_sequence_genome_pkl(chrom, start, end, genome)
        if refseq is '':
            print 'alignment is out of bounds, skipping'
            return 0, 0, 0, 0


        refseq = rev_complement(refseq)
        refseq = refseq[-nte_len:]
        for i in range(nte_len):
            if refseq[-i-1] == nte3[-i-1]:
                count +=1
                nte_bool = False
                # new_nte = ''
            else:
                # write cigar string from R --> L because need cigar string to be w.r.t positive strand for proper SAM format
                cigar_ext =  cigar_ext + str(count) + rev_complement(refseq[-i-1])
                if nte_bool is True:
                    new_nte = nte3[-i - 1] + new_nte  # adding NTE from L <-- R i.e: T  --> AT --> AAT
                count = 0

        # already had a 3'-end modification, we are going to re-write over it anyways
        if matches_count_list[0] is '0':
            cigar_string = cigar_string[2:]
        new_cigar = cigar_ext + cigar_string
        updated_sequence = rev_complement(add_nte3) + seq
    # print 'readID-{}\tnte3-{}\tnew_cigar={}'.format(readID, new_nte, new_cigar)

    cigar_element_list = cigar.split(':')  # ['MD', 'Z', '14C8']
    cigar_element_list[2] = new_cigar
    new_cigar_string = ':'.join(cigar_element_list)

    if nte_bool is False:
        new_nte = ''

    return new_nte, new_cigar_string, start, updated_sequence


def check_cigar_5nte(cigar, chrom, start, strand, seq, trim_dict, readID, genome):
    '''
    checks cigar string for 5'NTE modification, also looks into the trim_dict to see if there are modifications from trimming steps
    '''
    nte5 = ''
    cigar_string = cigar.split(':')[2]  # ['MD', 'Z', '14C8'] --> '14C8'
    matches_count_list = re.split(r'[ACTG]', cigar_string)  # ['14','8']
    add_nte5 = ''

    if len(matches_count_list) != 1:
        for i in range(len(matches_count_list)):
            # any 5'NTE for reads on the positive strand, will appear at the beginning of the cigar string, and appear as a 0 position
            # ex: ['22','0']
            if strand == '0':
                if matches_count_list[i] == '0':
                    nte5_base = seq[0 + i]
                    nte5 = nte5 + nte5_base
                # cigar does not contain any more NTEs
                else:
                    break
            # any 5'NTE for reads on the negative strand, will appear at the end of the cigar string
            else:
                if matches_count_list[-1 - i] == '0':
                    nte5_base = seq[-1 - i]
                    nte5 = nte5_base + nte5 # this is w.r.t the positive strand. will take Reverse complement of this after finished
                    nte5 = rev_complement(nte5)
                # cigar does not contain any more NTEs
                else:
                    break

    # 2. checks 3'NTE from trimmed dictionary
    if readID in trim_dict:
        nte_class_obj = trim_dict[readID]
        add_nte5 = nte_class_obj.nte5_trim
        original_sequence = nte_class_obj.original_RNA_sequence

    nte5 = add_nte5 + nte5 # this is currently w.r.t positive strand
    nte_len = len(nte5)

    # initializing variables
    count = 0
    cigar_ext = ''
    new_nte = ''
    nte_bool = True
    new_cigar = ''
    updated_sequence = ''

    if nte_len is 0:
        return new_nte, cigar, start, seq

    # 3. checks if the trimmed NTEs are true modifications ( i.e if there is more than 1 5'NTE, the mismatch must be consecutive)
    if strand == '0':
        # checking to see of the sequence is correct
        start = int(start) - len(add_nte5)
        end = int(start) + len(seq) + len(add_nte5)
        refseq = check_sequence_genome_pkl(chrom, start, end, genome)
        if refseq is '':
            print 'alignment is out of bounds, skipping'
            return 0, 0, 0, 0
        refseq = refseq[:nte_len]  # get the first letters of reference sequence that correspond to 5'NTE

        # 4. If 5'end modifications exists from the trimemd dictionary, it adds these modifications to the cigar string.
        for i in range(nte_len):
            if refseq[i] == nte5[i]:
                count +=1
                nte_bool = False

            else:
                cigar_ext =  str(count) + refseq[i] + cigar_ext
                if nte_bool is True:
                    new_nte = new_nte + nte5[i]  # because 5'nte, add from L --> R
                count = 0
        # rewritting cigar string, if 5'-end modification initially existed, remove and rewrite.
        if matches_count_list[0] is '0' :
            cigar_string = cigar_string[2:]
        new_cigar = cigar_ext + cigar_string
        updated_sequence = add_nte5 + seq

    # negative strand, so the 5'-NTE is at the beginning of SAM sequence and cigar.
    else:
        end = int(start) + len(seq) + len(add_nte5)
        refseq = check_sequence_genome_pkl(chrom, start, end, genome)

        if refseq is '':
            print 'alignment is out of bounds, skipping'
            return 0, 0, 0, 0

        refseq = rev_complement(refseq)
        refseq = refseq[:nte_len]
        for i in range(nte_len):
            if refseq[i] == nte5[i]:
                count +=1
                nte_bool = False

            else:
                # write cigar string from R --> L because need cigar string to be w.r.t positive strand for proper SAM format
                cigar_ext =  rev_complement(refseq[i]) + str(count) + cigar_ext
                if nte_bool is True:
                    new_nte = new_nte + nte5[i]  # adding NTE from R --> L
                count = 0

        # already had a 3'-end modification, we are going to re-write over it anyways
        if matches_count_list[-1] is '0':
            cigar_string = cigar_string[:-2]
        new_cigar = cigar_string + cigar_ext
        updated_sequence =  seq + rev_complement(add_nte5)

    cigar_element_list = cigar.split(':') # ['MD', 'Z', '14C8']
    cigar_element_list[2] = new_cigar
    new_cigar_string = ':'.join(cigar_element_list)

    if nte_bool is False:
        new_nte = ''

    return new_nte, new_cigar_string, start, updated_sequence


def internal_modification_count(cigar, nte3, nte5, num_mismatches_allowed=2):
    '''
    function to check if the ntes from the trimmed dict are true NTEs (consecutive NTEs)
    returns 3' and 5' NTEs and new cigar string

    '''
    mismatch_treshold_bool = False

    cigar_string = cigar.split(':')[2]  # ['MD', 'Z', '14C8'] --> '14C8'
    matches_count_list = re.split(r'[ACTG]', cigar_string)  # ['14', '8']

    internal_match_count = 0
    start = True

    if len(matches_count_list) - len(nte5) - len(nte3) - 1 >= num_mismatches_allowed:
        return True
    else:
        return False


def NTE_in_SAM(master_mapping_results, output_file, trim_dict, seq_file_type, genome, trim_type='three'):
    '''
    Function objective: Identify NTE's from alignments or trimming dictionary and add to alignment file.
    1) checking the cigar string
    2) checking the trimming dictionary for NTE (5' and 3')
    3) look up reference sequence to confirm NTE, this is important for consecutive NTEs
    4) modify cigar string to include NTEs if they are identified
    5) if sequence has 2 or more internal modifications, then remove alignment (not trustable)

    # trim type is either "init", "three", "five", or "both"

    '''
    final_output_list = []

    print trim_type
    for alignment in master_mapping_results:
        #skip all instances of unmapped reads
        # print alignment
        if alignment == '':
            continue
        alignment = alignment.split()

        strand = alignment[1]
        if strand == '4':
            continue
        else:
            readID, chrom, start, seq = alignment[0], alignment[2], alignment[3], alignment[9]

            # initializing NTEs as NULL
            nte3 = ''
            nte5 = ''
            mismatch_flag = alignment[13]
            cigar = alignment[12]

            # fetching the 3' end modification
            nte3, new_cigar_string, start, updated_sequence = check_cigar_3nte(cigar, chrom, start, strand, seq, trim_dict, readID, genome)

            # the alignment was out of bounds, skip
            if nte3 is 0 and new_cigar_string is 0:
                continue

            # fetching the 5' end modification
            nte5, new_cigar_string, start, updated_sequence = check_cigar_5nte(new_cigar_string, chrom, start, strand, updated_sequence, trim_dict, readID, genome)

            # the alignment was out of bounds, skip
            if nte5 is 0 and new_cigar_string is 0:
                continue
            # if trimmed nucleotides are true (consecutive NTEs, then edit cigar string)
            # print '{}\tnte3-{}\tnte5-{}'.format(readID, nte3, nte5)

            # will return True if more than 2 mismatches in the alignment
            exceeding_mismatches_bool = internal_modification_count(new_cigar_string, nte3, nte5, 2)
            if nte3 is '':
                nte3 = '0'
            if nte5 is '':
                nte5 = '0'

            if exceeding_mismatches_bool is True:
                continue
            else:
                alignment[3], alignment[6], alignment[7], alignment[8], alignment[9], alignment[12] = str(start), alignment[14], nte5, nte3, updated_sequence, new_cigar_string
                alignment.pop()
                output_alignment = '\t'.join(alignment)
                final_output_list.append(output_alignment)

    # for i in final_output_list:
    #     print i

    return final_output_list


def run_pipeline(fin, fout, directory, seq_file_type, analysis_type, idx_path, mapping_file, genome):

    rnd = 2  # MODIFY: MAKE THIS A USER DEFINED PARAMETER
    trim_dict = {}
    master_mapping_results = []
    # if input file is SAM format (indicates already one round of mapping), proceed to the trimming step
    sam = 0
    if seq_file_type == 'sam':
        sam = 1
        trim_dict, master_mapping_results = trimming_pipeline(fin, fout, trim_dict, seq_file_type, 'three', round_number = sam )
        seq_file_type = '-f'
        fin = '{}-trim0.fa'.format(fout)


    # first round of mapping is the initial mapping with no trimming
    trim_description = 'init-mapping'
    round_number = 0
    unmapped_fn, mapping_results_list = bowtie_mapping(fin, fout, mapping_file, idx_path, seq_file_type, round_number, trim_description = trim_description, trim_round_number = round_number )
    master_mapping_results.extend(mapping_results_list)

    # 3' trimming and mapping
    if seq_file_type == '-f': #sequence file is a fasta file format
        for round_number in range(1+sam,rnd+1): #if round=2, then this will loop through 1,2, and 3 rounds of mapping/trimming
            trim_description = 'trim3'
            trim_dict = trimming_pipeline(unmapped_fn, fout, trim_dict, seq_file_type, 'three', round_number)

            # trim_description = choice of "trim3", "trim5", or "both"
            unmapped_fn, mapping_results_list = bowtie_mapping(unmapped_fn, fout, mapping_file, idx_path, seq_file_type, round_number, trim_description = trim_description, trim_round_number = round_number )
            master_mapping_results.extend(mapping_results_list)
            print 'fin: {}; unmapped: {}'.format(fin,unmapped_fn)

        # for i in trim_dict:
        #     print i

        # removing the trimmed 3'end from dict from reads that did not map after 2 rounds of 3' end trimming.
        # we will try to map the original sequence but trimming from the 5'end
        NTE3_final_output = NTE_in_SAM(master_mapping_results, 'outputfile.sam', trim_dict, seq_file_type, genome, trim_type='three')
        print 'done with parsing NTE3 in SAM'


        fin = unmapped_fn
        trim_dict = {}
        master_mapping_results = []

        for round_number in range(1,rnd+1):
            trim_dict = trimming_pipeline(unmapped_fn, fout,trim_dict, seq_file_type, 'five', round_number)
            # trim_description = choice of "trim3", "trim5", or "both"
            trim_description = 'trim5'
            unmapped_fn, mapping_results_list = bowtie_mapping(unmapped_fn, fout, mapping_file, idx_path, seq_file_type, round_number, trim_description= trim_description, trim_round_number = round_number)
            master_mapping_results.extend(mapping_results_list)
            print 'fin: {}; unmapped: {}'.format(fin,unmapped_fn)
            #first  of 5' trimming
        # for i in master_mapping_results:
        #     print i
        NTE5_final_output = NTE_in_SAM(master_mapping_results, 'outputfile.sam', trim_dict, seq_file_type, genome, trim_type='five')
        print 'done with parsing NTE5 in SAM'
        NTE3_final_output.extend(NTE5_final_output)
        NTE5_final_output = []



        trim_dict = {}
        # last round of trimming, trim 1 nt on both the 5' and the 3' end
        master_mapping_results = []
        # trim_description = choice of "trim3", "trim5", or "both"
        trim_dict = trimming_pipeline(unmapped_fn, fout, trim_dict, seq_file_type, 'both', 1)
        unmapped_fn, mapping_results_list = bowtie_mapping(unmapped_fn, fout, mapping_file, idx_path, seq_file_type, 1 , trim_description = 'both', trim_round_number = 1)
        master_mapping_results.extend(mapping_results_list)

        trim_both_final_output = NTE_in_SAM(master_mapping_results, 'outputfile.sam', trim_dict, seq_file_type, genome, trim_type='both')

        NTE3_final_output.extend(trim_both_final_output)
        trim_both_final_output = []


    print 'printing output to: {}.NTE.sam'.format(fout)
    output= open('{}.NTE.sam'.format(fout), 'w')
    output.write('readID\tstrand\tchr\tstart\tX\tX\ttrimming_step\t5NTE\t3NTE\tbwt-seq\tquality\tsam_flags\n')
    print_lines = '\n'.join(NTE3_final_output)

    output.write(print_lines)

    return

def main():
    start_time = datetime.now()

    rnd = 2
    fin, fout, directory, seq_file_type, analysis_type, idx_path, mapping_file = get_arguments()
    if check_input_type(fin):
        print "input file is valid, opening {}".format(fin)
    start_time1 = datetime.now()
    filehandler = open("/u/home/k/kkoyano/nobackup-gxxiao3/scripts/genome.pickle", "r")
    genome = pickle.load(filehandler)
    filehandler.close()
    end_time1 = datetime.now()
    print "loading pickle took: {} ".format(end_time1 - start_time1)


    run_pipeline(fin, fout, directory, seq_file_type, analysis_type, idx_path, mapping_file, genome)


    os.system("rm {}-unmapped-*.fa".format(fout))
    # os.system('{}-NTE3_final_output.pkl'.format(fout))
    end_time = datetime.now()
    print 'Duration: {}'.format(end_time-start_time)
    print 'DONE!'



    return

if __name__ == "__main__":
	main()
