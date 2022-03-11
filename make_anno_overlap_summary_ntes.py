#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date:


!!!!!!!!!! NOTE: THIS SCRIPT ASSUMES HG19 !!!!!!!!!!!!!!

Description: Condense *.anno_overlap files so that each line is an annotation that contains all the NTEs
 Usage: python make_anno_overlap_summary_ntes.py -anno ./test/test.2mirnas.anno_overlap -summary ./test/test.2mirnas.summary -o . -b /u/nobackup/gxxiao3/kkoyano/annotations/mirbase/hg19.gff3.mature.mirs.with_seq.bed

or
# don't have to specify output file, will just output in the same directory as the anno_overlao
python make_anno_overlap_summary_ntes.py -anno ./test/test.2mirnas.anno_overlap -summary ./test/test.2mirnas.summary -b /u/nobackup/gxxiao3/kkoyano/annotations/mirbase/hg19.gff3.mature.mirs.with_seq.bed
Input:
1) input is a sam file annotation overlap file
2) annotation file containing coordinates and the sequence
3)

Output:
1) file with all annotations and mismatches counted
2)

Dependencies:
- python /u/home/k/kkoyano/nobackup-gxxiao3/scripts/GenomeFetch.py
- assumes hg19 is the build, to modify change function "get_annotation_sequence()" in helper_functions_mismatch_profile.py


Notes:

1) Decide if want to keep
    get_annotation_sequence_from_file() : uses a file to get the annotation sequence (user input as bed file) OR
    get_annotation_sequence() : uses genome fetch to get sequence

    LINES: 271 (as of 04/04..
    function --> annotation_overlap_mismatch_summary()
                ....
                 anno_seq = get_annotation_sequence_from_file(annotation_name)


2)


'''
import argparse
import os
import sys
import subprocess
# from subprocess import Popen, PIPE
import re
from itertools import permutations, izip
import gzip
import time
from datetime import datetime
# from helper_functions_mismatch_profile import get_annotation_sequence, get_total_reads


# BOUNDARY = 5  # this is the boundary extension to add SNPs to set -5 and +5 outside of annotation start and end coordinates, respectively
GENOME_FETCH_PATH = "/u/home/k/kkoyano/nobackup-gxxiao3/scripts/GenomeFetch.py"


def get_arguments():

    parser = argparse.ArgumentParser(description='ADD DESCRIPTION')

    parser.add_argument('-anno', dest='annotation_overlap', type=str,
                        help='name of annotation overlap')
    parser.add_argument('-summary', dest='summary_file', type=str,
                        help='name annotation summary file')
    parser.add_argument('-o', dest='output_dir', type=str,
                        help='name of output directory, will default to the current directory')
    # parser.add_argument('-b', dest='bed_file', type=str,
    #                     help='name of bed file that contains the reference sequence for all annotations (sequence must be in the 6th column)')
    args = parser.parse_args()

    annotation_overlap_file = args.annotation_overlap
    summary_file = args.summary_file
    output_dir = args.output_dir
    # bed_file = args.bed_file

    return annotation_overlap_file, summary_file, output_dir
    # return annotation_overlap_file, summary_file, output_dir, bed_file


class miRNA_annotation:
    """ A class that holds the information of the annotation, does not have to be a miRNA"""

    def __init__(self, chromosome, strand, start, end, species, name, seq, read_count, FIRST_X_READ_POSITIONS, LAST_X_READ_POSITIONS):
        self.chromosome = chromosome
        self.strand = strand
        self.start = int(start)
        self.end = int(end)
        self.species = species
        self.name = name
        self.perfect_match_reads_count = 0
        self.mismatch_reads_count = 0
        self.seq = seq  # actual RNA sequence, this is strand specific
        self.read_count = int(read_count)
        self.nte3 = {}  # key = {'genomic_position:REF>NTE':count, 'genomic_position:REF>NTE':count ...}

        self.nte5 = {}  # key = {'genomic_position:REF>NTE':count, 'genomic_position:REF>NTE':count ...}

        self.internal_mismatch = {}  # key = {'genomic_position:REF>NTE':count, 'genomic_position:REF>NTE':count ...}


def get_total_reads(summary_file):
    library_depth = ''
    with open(summary_file, 'r') as fin:
        for i in fin:
            if i.startswith("total number of reads:"):
                i = i.strip('\n').split('\t')
                library_depth = float(i[1])
                break
# total_read_count = get_total_reads(summary_file)
    return library_depth


def progress(progress_file):
    cmd = ['wc', '-l', progress_file]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    num_lines = int(stdout.strip().split()[0])

    percentages_time = []
    for i in range(0, 10):
        percentages_time.append(num_lines * i / 10)
    # print percentages_time
    return percentages_time


def reverse_complement(sequence):
    complement_dict = {'A': 'T',
                       'T': 'A',
                       'C': 'G',
                       'G': 'C',
                       'N': 'N'}
    sequence = sequence.upper()
    rev_seq = sequence[::-1]
    rev_compl = ''
    for base in rev_seq:
        rev_compl += complement_dict[base]
    return rev_compl


def get_annotation_sequence(annotation_name, chromosome, strand, start, end):
    # The function gets the sequence of interest based on it's coordinates and strand, calling sequencines i 1 based

    cmd = "python {} -o hg19 -c {} -s {} -f {} -t {}".format(GENOME_FETCH_PATH, chromosome, strand, start, end)
    cmd = cmd.split()
    bash_result = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    annotation_sequence = bash_result.communicate()[0]
    annotation_sequence = annotation_sequence.strip('\n')
    # print 'annotation seq for {}: {}'.format(annotation_name, annotation_sequence)
    #'AAAAGCTGGGTTGAGAGGGCGA'
    return annotation_sequence


def get_annotation_sequence_from_file(annotation_name, bed_file):
    # bed_file = '/Users/kiku/mount/annotations/mirbase/hg19.gff3.mature.mirs.with_seq.bed'
    cmd = 'grep -m 1 {} {}'.format(annotation_name, bed_file)
    cmd = cmd.split(' ')
    bash_result = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    annotation_sequence = bash_result.communicate()[0]
    if annotation_sequence == '':
        print 'did not have annotation sequence'
        return
    annotation_sequence = annotation_sequence.strip('\n').split('\t')[6]

    # print annotation_name, annotation_sequence

    return annotation_sequence


def get_nte3_position(ref, read_nte, read_start, read_end, anno_strand):

    if anno_strand == "+":
        if len(read_nte) > 1:  # if the nte is more than one nucleotide then need to have position cover 2 genomic positions
            shift_start_position = len(read_nte) - 1  # ex: shift 1 position to the left if there is 2 NTEs
            start_position = int(read_end) - shift_start_position
            position = str(start_position) + '-' + read_end

        else:
            position = read_end

        # reference = ref[::len(read_nte)]
        reference = ref[-len(read_nte):]  # get the last 2 nucleotides in the reference string, which are the NTEs

    else:  # negative strand
        if len(read_nte) > 1:
            shift_start_position = len(read_nte) - 1  # ex: shift 1 position to the left if there is 2 NTEs
            start_position = int(read_start) + shift_start_position
            position = str(start_position) + '-' + read_start
        else:
            position = read_start
        reference = ref[:len(read_nte)]  # get the first 2 nucleotides in the reference string, which are the 3NTEs, but since neg strand, get the reverse complement
        reference = reverse_complement(reference)
    return reference, position


def update_nte_dict(nte_dict, reference, position, read_nte, read_count):
    """
    update nte_dictionary with key ('position:ref>nte') and value (read_count) so we can give it back to the annotation object

    """
    position = str(position)
    join_string = ''.join([reference, '>', read_nte])
    key = ':'.join([position, join_string])  # '22102508:TT>AA'
    # print key
    if key in nte_dict:
        nte_dict[key] = nte_dict[key] + read_count
    else:
        nte_dict[key] = read_count

    return nte_dict


def update_nte_info(annotation_obj, anno_strand, read_start, read_end, read_nte5, read_nte3, cigar, read_seq, read_count):
    '''
    update annotation class object
    cigar = MD:Z:0A0C20

    '''

    # annotation_obj.read_count = annotation_obj.read_count + read_count
    ### add to nte3, nte5 and internal #
    # key = 'position:nte'

    mismatch_count = annotation_obj.mismatch_reads_count
    mismatch_count = mismatch_count + read_count
    annotation_obj.mismatch_reads_count = mismatch_count

    cigar = cigar.split(':')[2]  # MD:Z:0A0C20 --> '0A0C20'
    x = re.split("[ACTG]", cigar)  # ex: ['0', '20']  #2 5' ex: ['0', '0', '20']
    ref = re.sub("[^a-zA-Z]", "", cigar)  # 'AC'
    original_ref = ref

    reference, position = '', 0
    if read_nte3 != "0":
        # if the read has a 3'NTE, then find out the 1) positions of the mismatches 2) NTE nucleotides 3) the reference nucleotides and
        reference, position = get_nte3_position(ref, read_nte3, read_start, read_end, anno_strand)

        # updating nte dictionary counts
        nte3_dict = update_nte_dict(annotation_obj.nte3, reference, position, read_nte3, read_count)
        #(nte_dict, reference, position, read_nte, read_count):

        annotation_obj.nte3 = nte3_dict  # updated

    if read_nte5 != "0":
        # key = 'position:nte', value = count
        # if the read has a 5'NTE, then find out the 1) positions of the mismatches 2) NTE nucleotides 3) the reference nucleotides and
        if anno_strand == "+":
            position = read_start
            if len(read_nte5) > 1:  # if the NTE is more than 1 nucleotide, record the position as a window (i.e NTE_start-NTE_end)

                position = str(read_start) + '-' + str(int(read_start) + len(read_nte5) - 1)

            reference = ref[:len(read_nte5)]  # get the first N nucleotides of the reference string
        else:  # strand is negative
            position = read_end
            if len(read_nte5) > 1:
                start_position = int(position) - len(read_nte5) + 1
                position = position + '-' + str(start_position)
            reference = ref[-len(read_nte5):]               # get the last 2 nucleotides of hte string
            reference = reverse_complement(reference)  # the reference nucleotides are from the cigar string, the output of the cigar string is always with respect to the positive strand, so we need to get the reverse complement to get the reference of the actual RNA

        # updating nte dictionary
        nte5_dict = update_nte_dict(annotation_obj.nte5, reference, position, read_nte5, read_count)
        annotation_obj.nte5 = nte5_dict  # updated dictionary for the annotation object

    # if the read has no NTE then set string to 0, this is because we will use the length of the NTE string to count how many NTEs there are, currently the NTE is set to 0 if no NTE.
    if read_nte3 == '0':
        read_nte3 = ''
    if read_nte5 == '0':
        read_nte5 = ''

    x_temp = x                          # example: ['0', '0', '4', '20', '0'], a list of the positions in the cigar string
    if anno_strand == '+':
        shift_start = len(read_nte5)  # used to remove the beginning 0s in the cigar string
        shift_end = len(read_nte3)  # used to remove the ending 0s in the cigar string
    else:  # negative strand
        shift_start = len(read_nte3)  # used to remove the beginning 0s in the cigar string, because RNA is on other strand, the NTE should be flipped.
        # positive strand -->        [0, 0 , 4, 20]  ; has 2 5' NTEs
        # negative strand cigar --> [0, 0 , 4, 20]  ; has 2 3' NTEs
        shift_end = len(read_nte5)  # used to remove the ending 0s in the cigar string

    ref = list(ref)  # ['A', 'A', 'C', 'T']
    for i in range(shift_start):
        x_temp.pop(0)  # ['0', '0', '4', '20', '0'] -> #['0', '4', '20', '0'] -> ['4', '20', '0']
        ref.pop(0)  # ['A', 'A', 'C', 'T'] -> ['A', 'C', 'T'] ->['C', 'T']
    for i in range(shift_end):
        x_temp.pop()  # ['4', '20', '0'] --> #['4', '20']
        ref.pop()  # ['C', 'T'] -> ['C']

    if len(x_temp) > 1:  # if value > 1,then there is some internal mismatch # "2G17" -> [2, 17]
        # print x_temp
        internal_mismatch = annotation_obj.internal_mismatch
        index_counter = 0

        # print x
        for m in range(len(x_temp) - 1):  # positive strand = [0,0,4,12]
            minor_allele = ''
            # print read_seq, int(x_temp[m]), shift_start, index_counter
            position = int(read_start) + int(x_temp[m]) + shift_start
            minor_allele = read_seq[int(x_temp[m]) + shift_start + index_counter]  # find the minor allele
            reference = ref[m]
            if anno_strand == "-":  # strand is negative need to get compl.
                reference = reverse_complement(reference)
                minor_allele = reverse_complement(minor_allele)

            # updating dictionary
            internal_mismatch_dict = update_nte_dict(annotation_obj.internal_mismatch, reference, position, minor_allele, read_count)
            annotation_obj.internal_mismatch = internal_mismatch_dict
            index_counter += 1
    return annotation_obj


def dict_to_col(nte_dict):
    # key = 'position:NTE', value = count
    # return 2 values
    #   1. col_line - string of  'position1:ref>nte:count|position2:ref>nte:count|....|position3:ref>nte:count'
    #   2. total NTE counts for the

    if len(nte_dict) == 0:  # if there are no mismatches, return 0 to fill column position
        col_line = '0'
        return col_line, 0
    else:
        total_nte_count = 0
        col_list = ['{}:{}'.format(key, nte_dict[key]) for key in nte_dict]
        for key in nte_dict:
            total_nte_count += int(nte_dict[key])
        col_line = '|'.join(col_list)
        return col_line, total_nte_count


def write_output(annotation_overlap_file, annotation_obj_dict, output_dir, library_depth):
    """

    printing output file.
        def __init__(self, chromosome, strand, start, end, species, name, seq, read_count, FIRST_X_READ_POSITIONS, LAST_X_READ_POSITIONS):
        self.chromosome = chromosome
        self.strand = strand
        self.start = int(start)
        self.end = int(end)
        self.species = species
        self.name = name
        self.perfect_match_reads_count = 0
        self.mismatch_reads_count = 0
        self.seq = seq  # actual RNA sequence, this is strand specific
        self.read_count = int(read_count)
        self.nte3 = {}  # key = {'genomic_position:REF>NTE':count, 'genomic_position:REF>NTE':count ...}

        self.nte5 = {}  # key = {'genomic_position:REF>NTE':count, 'genomic_position:REF>NTE':count ...}

        self.internal_mismatch = {}  # key = {'genomic_position:REF>NTE':count, 'genomic_position:REF>NTE':count ...}


    =====================
    output columns
    =====================

    1. annotation_name
    1. species
    2. chromosome
    3. strand
    4. start
    5. end
    6. sequnece
    7. total_read_count
    8. perfect_match_count
    9. read_mismatch_count
    10. RPM
    11. nte5
    12. nte3
    13. internal modifications

    """

    # outputfn = ''

    header = ['annotation_name', 'species', 'chromosome', 'strand', 'start', 'end', 'anno_sequence',
              'total_read_count', 'perfect_match_count', 'anno_mismatch_count', 'RPM',
              'total_nte5_count', 'total_nte3_count', 'total_internal_count',
              'nte5', 'nte3', 'internal_mismatch\n']

    header = '\t'.join(header)

    output_fn = ''
    if output_dir == '':
        output_fn = annotation_overlap_file + '.anno_summary'
    else:
        if output_dir == '.':
            cwd = os.getcwd()
            output_fn = cwd + '/' + annotation_overlap_file.split('/')[-1] + '.anno_summary'
        else:
            if not output_dir.endswith('/'):
                output_fn = output_dir + '/' + annotation_overlap_file + '.anno_summary'
            else:
                output_fn = output_dir + annotation_overlap_file + '.anno_summary'
    with open(output_fn, 'w') as out_fn:
        out_fn.write(header)
        for key in annotation_obj_dict:  # key = miRNA name, value = mirna object
            annotation_obj = annotation_obj_dict[key]
            nte3_dict = annotation_obj.nte3
            internal = annotation_obj.internal_mismatch
            nte5_dict = annotation_obj.nte5

            chromosome = annotation_obj.chromosome
            strand = annotation_obj.strand
            start = annotation_obj.start
            end = annotation_obj.end
            species = annotation_obj.species
            name = annotation_obj.name
            perfect_match_reads_count = annotation_obj.perfect_match_reads_count
            mismatch_reads_count = annotation_obj.mismatch_reads_count
            sequence = annotation_obj.seq
            read_count = annotation_obj.read_count
            rpm = int(read_count) * 1000000.0 / int(library_depth)
            rpm = round(rpm, 1)

            nte5_col, total_nte5_count = dict_to_col(nte5_dict)  # returns 2 strings, (1) long string of all 'pos:ref>nte:count' and are tab delimited, and the (2) is the total count for 5' NTEs
            nte3_col, total_nte3_count = dict_to_col(nte3_dict)
            internal_col, total_internal_count = dict_to_col(internal)

            annotation_output_line = [key, species, chromosome, strand, start, end, sequence,
                                      read_count, perfect_match_reads_count, mismatch_reads_count, rpm,
                                      total_nte5_count, total_nte3_count, total_internal_count,
                                      nte5_col, nte3_col, internal_col]
            annotation_output_line = [str(i) for i in annotation_output_line]
            annotation_output_line = '\t'.join(annotation_output_line)
            out_fn.write(annotation_output_line)
            out_fn.write('\n')
    return


def annotation_overlap_mismatch_summary(annotation_overlap_file, output_dir, library_depth, FIRST_X_READ_POSITIONS=10, LAST_X_READ_POSITIONS=-10):
    # def annotation_overlap_mismatch_summary(annotation_overlap_file, output_dir, library_depth, bed_file, FIRST_X_READ_POSITIONS=10, LAST_X_READ_POSITIONS=-10): # removing bed file dir
    """
    This function goes through the *.anno_overlap file and collects the annotation information and the
    1) 3NTE, 5NTE and internal modification positions, types and counts,
    2) the mismatch position counts (10 positions after start site and 10 positions before end site) relative to the read
    3) snp summary stats. i.e. how many 3'NTEs were lost due to SNPs and what was the mismatch nucleotide,

    input:
    chr8    -   22102487    22102508    251397|108  0   0   NM:i:0  MD:Z:22 TTCGCCCTCTCAACCCAGCTTT  hsa-miR-320a    chr8    22102488    22102509mirna
    chr8    -   22102485    22102508    251398|3    0   TT  NM:i:1  MD:Z:0T0T22 AATTCGCCCTCTCAACCCAGCTTT    hsa-miR-320a    chr8    2210248822102509    mirna
    chr8    -   22102486    22102507    251399|10   0   C   NM:i:1  MD:Z:0T21   GTTCGCCCTCTCAACCCAGCTT  hsa-miR-320a    chr8    22102488    22102509    mirna
    chr8    -   22102486    22102508    251400|11   T   T   NM:i:1  MD:Z:0T21T0 ATTCGCCCTCTCAACCCAGCTTA hsa-miR-320a    chr8    22102488    22102509    mirna


    for each annotation species (ex: hsa-miR-1, hsa-miR-2...etc.)
    There will be a dictionary for all the annotation
        key=annotation name
        value =annotation class object that contains all the information needed to print out information
    *See class miRNA_annotation* -althoough could be any annotation, i just named it this for clarification

    """
    annotation_dict = {}

    # dictionary where the key is the annotation NAME and the value is the annotation class object
    # that carries all the other information (i.e. chrom, annotation start/end, strand, read_count, nte5, nte3 and internal mismatches)
    progress_bar = progress(annotation_overlap_file)
    with open(annotation_overlap_file, 'r') as fin:

        progress_count = 0
        for i in fin:
            if i.startswith('#'):
                continue
            i = i.strip().split()
            annotation_name, anno_chrom, anno_start, anno_end, anno_strand, species = i[10], i[11], i[12], i[13], i[1], i[14]
            read_start, read_end, read_nte5, read_nte3, mismatch_number, cigar, read_seq = i[2], i[3], i[5], i[6], i[7], i[8], i[9]

            read_count = int(i[4].split('|')[1])
            mismatch_number = mismatch_number.split(':')[2]
            # Adding annotation to dictionary
            # annotation_dict = {'hsa-miR-320a': <miRNA_annotation object>, 'hsa-miR-379-5p':<miRNA_annotation object>}
            if annotation_name not in annotation_dict:
                anno_seq = get_annotation_sequence(annotation_name, anno_chrom, anno_strand, anno_start, anno_end)
                # anno_seq = get_annotation_sequence_from_file(annotation_name, bed_file)
                annotation_obj = miRNA_annotation(anno_chrom, anno_strand, anno_start, anno_end, species, annotation_name, anno_seq, read_count, FIRST_X_READ_POSITIONS, LAST_X_READ_POSITIONS)  # chrom, strand, start, end, species, anno_name, seq
                annotation_dict[annotation_name] = annotation_obj

            else:
                annotation_obj = annotation_dict[annotation_name]
                add_count = annotation_obj.read_count
                add_count = add_count + read_count
                annotation_obj.read_count = add_count

            ###################
            # Counting NTES/Modifications
            ###################
            # if there are 15 fields, then look at the last field to see if it has internal modifications, this field is the same as nte3 and nte5, where it displays the minor nucleotide, with respect to the **RNA** not the positive strand. Also if there is no internal modification, the field will be '0'
            # if len(x) == 15:
            #    mismatch_number == i[14]
            if mismatch_number == '0' and read_nte5 == '0' and read_nte3 == '0':  # mismatch_number  = NM:i:0, no modifications or NTE to add, just add read count
                add_count = annotation_obj.perfect_match_reads_count
                add_count = add_count + read_count
                annotation_obj.perfect_match_reads_count = add_count
                continue

            else:
                annotation_obj = update_nte_info(annotation_obj, anno_strand, read_start, read_end, read_nte5, read_nte3, cigar, read_seq, read_count)
                annotation_dict[annotation_name] = annotation_obj
                ### add to nte3, nte5 and internal #

            # progress markers
            progress_count += 1
            if progress_count in progress_bar:
                x = progress_bar.index(progress_count)
                percentage_val = (x + 1) * 10
                print '{}% Complete'.format(percentage_val)
            # progress markers

    ## print statements for results #
    # for key in annotation_dict:  # mirna
    #     annotation_obj = annotation_dict[key]
        # print 'name = {}; count = {}'.format(annotation_obj.name, annotation_obj.read_count)

        # internal = annotation_obj.internal_mismatch
        # print internal

    write_output(annotation_overlap_file, annotation_dict, output_dir, library_depth)

    return


def main():
    """
    FINISHED CORRECT OUTPUT, NEED TO ADJUST INPUT AND OUTPUT PARAMETRES, BUT OTHERWISE READY TO GO!
    """
    #main run section, comment out for now so i can test the script #
    print 'Start script....'
    print 'This script assumes HG19 and also has a dependency on GenomeFetch.py'
    start_time = datetime.now()
    # annotation_overlap_file, summary_file, output_dir, bed_file = get_arguments()
    annotation_overlap_file, summary_file, output_dir = get_arguments()

    if output_dir == None:
        output_dir = ''

    library_depth = get_total_reads(summary_file)
    # print 'total reads:\t{}'.format(library_depth)
    FIRST_X_READ_POSITIONS = 10
    LAST_X_READ_POSITIONS = 10
    annotation_overlap_mismatch_summary(annotation_overlap_file, output_dir, library_depth, FIRST_X_READ_POSITIONS, LAST_X_READ_POSITIONS)
    # annotation_overlap_mismatch_summary(annotation_overlap_file, output_dir, library_depth, bed_file, FIRST_X_READ_POSITIONS, LAST_X_READ_POSITIONS) #removing bed file
    end_time = datetime.now()

    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'

    # testing section
    # annotation_overlap_file = "./test/test.2mirnas.anno_overlap"
    # output_dir = ''
    # library_depth = 1200000
    # FIRST_X_READ_POSITIONS = 10
    # LAST_X_READ_POSITIONS = 10
    # bed_file = "/Users/kiku/mount/annotations/mirbase/hg19.gff3.mature.mirs.with_seq.bed"
    # annotation_overlap_mismatch_summary(annotation_overlap_file, output_dir, library_depth, bed_file, FIRST_X_READ_POSITIONS, LAST_X_READ_POSITIONS)

    return


if __name__ == "__main__":
    main()
