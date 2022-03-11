"""
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 08/30/2018

Description: Helper functions fro the main script: trimming_pipeline.py

trimming_pipeline.py takes in bowtie parameters and fasta or fastq files and performs initial mapping and 2 rounds of 3' trimming, and 2 rounds of 5' trimming

Input:

Output:

Notes:

"""
import argparse, os, sys
from subprocess import Popen, PIPE
import re
from itertools import permutations
import time
from datetime import datetime

def check_fasta_format(line2, line3):
    line2 = line2.strip()
    line3 = line3.strip()
    line2 = set(line2)
    nucleotides = set('ATCGN')
    if line2.difference(nucleotides):
        print 'this is not a fasta file, letters in first line are not all nucleotides'
        sys.exit()
    if line3.startswith('>'):
        return True
    print 'this is not a fasta file'
    sys.exit()


def check_fastq_format(line2, line3):
    """
    @label
    sequence
    +other
    quality_scores

    @SRR073731.13 PATHBIO-SOLEXA1:7:1:25:1994 length=38
    CAAGGTGAGGAAAGGGAAACAAGAATTCTTTTT
    +SRR073731.13 PATHBIO-SOLEXA1:7:1:25:1994 length=38
    BB@CC;B@BC=:4>=A39BABC>9A9?4:BBB>

    """
    line2 = line2.strip()
    line3 = line3.strip()
    line2 = set(line2)
    nucleotides = set('ATCGN')
    if line2.difference(nucleotides):
        print 'this is not a fastq file, letters  in first line are not all nucleotides'
        sys.exit()
    if line3.startswith('+'):
        return True
    print 'this is not a fastq file'
    sys.exit()


def check_sam_format(line2, line3):
    line2 = line2.strip()
    line3 = line3.strip()
    line2 = set(line2)
    nucleotides = set('ATCGN')
    if line2.difference(nucleotides):
        print 'this is not a fasta file, letters are not all nucleotides'
        sys.exit()
    if line3.startswith('>'):
        return True
    print 'this is not a fasta file'
    sys.exit()


def check_sam_flags(cigar, mismatch_flag, alignment):
    """
    Checking to see if the cigar is in the right location?
    """
    cigar_idx = 12
    if not cigar.startswith('MD:Z:') or mismatch_flag.startswith('NM:i:'):
        cigar_idx = [ i for i, item in enumerate(alignment) if re.search('MD:Z:', item) ]
        mismatch_idx = [ i for i, item in enumerate(alignment) if re.search('NM:i:', item) ]
        if len(cigar_idx) != 1:
            print alignment
            print cigar
            print mismatch_flag
            print cigar_idx
            print 'not correct sam format, need cigar index of --> NM:i:X\n please put in correct sam format '
            sys.exit(1)
        else:
            mismatch_flag = alignment[mismatch_idx[0]]
            cigar = alignment[cigar_idx[0]]
    return (
     cigar, cigar_idx[0], mismatch_flag)


def rev_complement(sequence):
    sequence = sequence.upper()
    complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    rev_seq = sequence[::-1]
    rev_compl = ''
    for i in rev_seq:
        rev_compl += complement[i]

    return rev_compl


def complement(sequence):
    sequence = sequence.upper()
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    complement = ''
    for i in sequence:
        complement += complement_dict[i]

    return complement


def check_3NTE(chrom, start, seq, position):
    p = position - len(seq)
    alt_nucleotide = seq[p:]
    ref_nucleotide = check_sequence(chrom, start, seq)[p:]
    print ('alt= {}, ref = {}').format(alt_nucleotide, ref_nucleotide)
    for i in range(len(ref_nucleotide)):
        if ref_nucleotide[i] == alt_nucleotide[i]:
            return False
            continue

    NTE3 = alt_nucleotide
    return NTE3


def check_5NTE(chrom, start, seq, position):
    p = position + 1
    alt_nucleotide = seq[:p]
    ref_nucleotide = check_sequence(chrom, start, seq)[:p]
    print ('alt= {}, ref = {}').format(alt_nucleotide, ref_nucleotide)
    for i in range(len(ref_nucleotide)):
        if ref_nucleotide[i] == alt_nucleotide[i]:
            return False
            continue

    NTE5 = alt_nucleotide
    return NTE5


def check_input_type(input_file):
    file_exists(input_file)
    with open(input_file, 'r') as (fin):
        line = fin.readline()
        line2, line3 = fin.readline(), fin.readline()
        if line.startswith('@'):
            check_fastq_format(line2, line3)
        if line.startswith('>'):
            line2, line3 = fin.readline(), fin.readline()
            check_fasta_format(line2, line3)
    return True


def run_command(cmd, stdin=0):
    if not isinstance(cmd, list):
        cmd = cmd.split(' ')
    if stdin == 0:
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
    else:
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
        stdout, stderr = p.communicate(stdin)
    if p.returncode != 0:
        sys.stderr.write('command failed')
        sys.stderr.write(stdout)
        sys.stderr.write(stderr)
        sys.exit(1)
    return (stdout, stderr)


def file_exists(fname):
    if os.path.isfile(fname):
        return
    print ('ERROR!! \n\nfilename {} does not exist').format(fname)
    sys.exit(1)

