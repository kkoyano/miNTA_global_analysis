#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 06/28/2018

Description: script to run for each sample, filter modificaitons based on the phred score


Input:

Output:

Notes:
07/23/2018

CANNOT JUST GREP MIRNA MODIFICATIONS, POSSIBLE FOR MODIFICATIN TO OCCUR ON BOTH 5' AND INTERNAL OR 3' AND INTERNAL, THIS RESULTS IN DOUBLE COUNTING IN THE PIPELIN


'''
import argparse
import os
import sys
import subprocess
import re
from itertools import permutations, izip
import gzip
import time
from datetime import datetime

# f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')


class modification_class:
    """ A class that holds the information of the modifications, does not have to be a miRNA"""

    def __init__(self, miRNA_name, modification_profile, RNA_modification_index):
        self.miRNA_name = miRNA_name
        self.modification_profile = modification_profile #'chr:strand:position:ref>alt'
        self.RNA_modification_index = RNA_modification_index #[20] or [19,20] #if there are more than 2 numbers, this means the length of the modificaion/NTE
        self.modification_count = 0 # modification count


def get_arguments():

    parser = argparse.ArgumentParser(description='Filter out modifications based on their phred score that is in the fastq file. default is phred 33')

    parser.add_argument('-a ', dest='anno_summary', type=str,
                        help='path of annotation summary file (where the rows are genes/annotations) ')
    parser.add_argument('-q ', dest='fastq_file', type=str,
                        help='path of fastq file')
    parser.add_argument('-s ', dest='sam_file', type=str,
                        help='path of *.NTE.sam file, where the rows are condensed reads (readID|abundance))')
    parser.add_argument('-phred ', dest='phred_score', type=int,
                        help='phred_score threshold, filter modifications that have a phred less than this value (number 0-42), based on phred +33 ')
    parser.add_argument('--phred64', action='store_true', default=False,
                        help='flag that indicates if the FQ file is based on phred-64 score ')
    parser.add_argument('--test', action='store_true', default=False,
                        help='flag that indicates if I am running a test ')
    args = parser.parse_args()

### Testing Section

    test = args.test
    phred_score = 30
    # test = True
    phred64 = args.phred64
    if test is True:
        anno_summary= "/Users/kiku/test_files/99_58_SER.NTE.sam.anno_overlap.t0.5.canonical_end.nte_fields.mirna.ALL.psl.best_score.bed.results.lvl1.anno_summary.nosnp.credible.rmsk.snv-m45M55.credible-r2-N5"
        fastq_file="/Users/kiku/test_files/99_58_SER.noadapters.fastq.gz"
        sam_file= "/Users/kiku/test_files/99_58_SER.NTE.sam.anno_overlap.t0.5.canonical_end.nte_fields.mirna.ALL.psl.best_score.bed.results.lvl1"
        print "This is a test"

        print "phred score threshold, remove if phred is less than {}".format(phred_score)
        return anno_summary, fastq_file, sam_file, phred_score, phred64

    anno_summary = args.anno_summary
    fastq_file = args.fastq_file
    sam_file = args.sam_file
    phred_score = args.phred_score
    phred64 = args.phred64
    check_file_exist(anno_summary)
    check_file_exist(fastq_file)
    check_file_exist(sam_file)

    print "phred score threshold, remove if phred is less than {}".format(phred_score)

    if 0 < phred_score > 42:
        sys.exit("phred score is not valid, musct pick value between 0-42")

    if phred64 is True:
        print "TODO: Calculating Phred score based on phred 64"
        sys.exit("HEY KIKU! You Need to implement this option! Exiting...")
    else:
        print "Calculating Phred score based on phred 33"

    return anno_summary, fastq_file, sam_file, phred_score, phred64


def check_file_exist(file_path):
    if os.path.exists(file_path) is False:
        print 'ERROR!!!'
        output_err_msg = '{} does not exist'.format(file_path)
        sys.exit(output_err_msg)
    return

# function modfied from add_modification_fields_to_SAM.py

def reverse_complement(sequence):
    '''
    get the reverse complement of given sequence
    '''
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


# end

def add_modification_set(modification_set, nte_profile, chrom, strand):
    '''
    modification set is a set() data type
    '''

    if nte_profile is '0':
        return modification_set
    else:
        nte_profile_list = nte_profile.split('|')
        for i in nte_profile_list: # ex: 50004088:C>A:1
            mod_elements = i.split(':')
            mod_elements = [chrom, strand] + mod_elements
            mod_elements.pop() # removing the count
            x = ':'.join(mod_elements)
            # modification_set.set(x) #x = 'chr:strand:position:ref>alt'
            modification_set.append(x) #x = 'chr:strand:position:ref>alt'

    return modification_set


def get_read_index_position(modification_profile, read_start, read_end, strand):
    '''
    gets index position of the modification in the read, this is obtained by using the read's mapped genomic coordinates relative to the position of the modification.

    The index position returned is with respect to the RNA!!!!

    This means that if the read is mapped to the negative strand, the index is the RNA, not the reverse complement sequence that's observed in the SAM file.

    will return list of modifications

    '''
    # modification_profile = 'chr1:+:41220050:A>C'

    position = modification_profile.split(':')[2]

    position_list= position.split('-')

    # handling if the NTE covers more than one position (in this case the RNA modification index position will be a list of all indexes with modifications)
    if len(position_list) > 1:
        if strand is "+" or strand is "0":
            positions = range(int(position_list[0]), int(position_list[1]) + 1)
            RNA_modification_index_position = [i - int(read_start) for i in positions]

        else:
            positions = range(int(position_list[1]), int(position_list[0]) + 1)
            RNA_modification_index_position = [int(read_end) - i for i in positions]
            # do I need to reverse list ? x[::-1]

        return RNA_modification_index_position
    else:
        if strand is "+" or strand is "0":
            RNA_modification_index_position = [int(position) - int(read_start)]

        else:  # strand is negative/16
            RNA_modification_index_position = [int(read_end) - int(position)]

        return RNA_modification_index_position


def filter_quality_score2(fastq_file, check_fastq_seqs, phred_score_cutoff):
    '''
    check_fastq_seqs = dictionary, key = RNA_sequence, value = list of modificaiton class objects. each object contains the miRNA name, modification profile ('chr:strand:position:ref>minor'), RNA modification index, read_count

    updated version where read_count is now filtered for only modifications > phred30 (or whatever user specifies)

    next: go through check_fastq_seq and modify

    d = {'modificationID':updated_read_count, ....}
    then we can use this to go through the annotation summary file, if the modification comes up in this, then we will modify the read count

    '''
    progress_list = progress(fastq_file)

    fin_obj = gzip.open(fastq_file, 'rb') if re.search(r'\.gz$', fastq_file) else open(fastq_file, 'r')

    if check_fastq_seqs is '0':
        return '0'
    count = 0
    print 'going through fastq file'
    for i in fin_obj:
        # first line is the readID
        # line 2 is the read sequence

        # @SRR5507957.1 HWI-ST388-W7D:192:C5JBDACXX:2:1101:1394:2208 length=40  ---> i
        # TCAAACTCTGAATGCCTAAAGCCAACGTA                                         ---> seq *
        # +SRR5507957.1 HWI-ST388-W7D:192:C5JBDACXX:2:1101:1394:2208 length=40  ---> skip
        # ?=?DDDDFHDFHFCEHG@CABIHG?AEAA                                         ---> quality_seq *

        seq = next(fin_obj)
        seq = seq.strip('\n')

        if seq in check_fastq_seqs:
            skip = next(fin_obj) # this is some other line we can ignore
            quality_seq = next(fin_obj) # quality score line

            # print quality_seq
            # go through the each index and calculate qscore
            modification_obj_list =  check_fastq_seqs[seq]

            # go through all modifications associated with this sequence and add 1 if the phred is greater than 30
            for modification_obj in modification_obj_list:
                index_list = modification_obj.RNA_modification_index
                # print index_list
                q_score_count = 0
                for index in index_list:
                    ascii_char = quality_seq[index]
                    q_score_num = ord(ascii_char)-33
                    if q_score_num < phred_score_cutoff:
                        # print 'failed index-{} is phred score:{} ascii:{}'.format(index, q_score_num, ascii_char)
                        break
                    else:
                        q_score_count += 1

                    if q_score_count == len(index_list):
                        modification_obj.modification_count += 1

                # print '{} : {} '.format(seq, check_fastq_seqs[seq])
                # skip2 = next(fin_obj) # sequence ID name of next sequence, read this line to queue up next actual read
                count += 4
        else:
            next(fin_obj)
            next(fin_obj)
            count += 4
        # tracking progress, this should be the time consuming part.

        if count in progress_list:
            percent_index = progress_list.index(count)
            percentage_val = (percent_index + 1) * 10
            print '{}% completed'.format(percentage_val)


    fin_obj.close()
    # this is a dictionary with all of the read_seuqences as the key and the value is another dictionary where the key = index, and the value is a list of modification objects, all with the updated version of now many reads
    return check_fastq_seqs




def get_read_sequences_from_annosummary(mir_mod_list, sam_file, fastq_file, phred_score_cutoff, nte_type= "nte3"):
    # nte_type= is either a choice of "nte3", "nte5", or "internal", this will decide
    check_fastq_seqs = {}

    if len(mir_mod_list) is 0:
        return '0'

    updated_nte_count = 0

    for i in mir_mod_list:
        ## UPDATE, NEED TO PUT BACK TOGETHER THE MODIFICATION AFTER FILTERING OUT FOR LOW QUALITY MODIFICATIONS
        # get all lines in SAM File that contain the modification profile of interest, only look grp from the specific column the NTE is found in, sometimes the NTE will be both an internal modification and also a 3' or 5' end modificaiton, otherwise we will double count.

        if nte_type is "nte5":
            cmd = "cut -f 1-16 {}".format(sam_file)
        elif nte_type is "nte3":
            cmd = "cut -f 1-15,17 {}".format(sam_file)
        elif nte_type is "internal":
            cmd = "cut -f 1-15,18 {}".format(sam_file)
        else:
            print "unknown nte type, check function get_read_sequences_from_annosummary"

        cmd = cmd.split(' ')

        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['grep', '-w', i ], stdin=p1.stdout, stdout=subprocess.PIPE)

        # cmd = "grep -w {} {}".format(i, sam_file)
        # cmd = cmd.split(' ')
        # p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = p2.communicate()[0]
        sam_lines = stdout.split('\n')
        sam_lines = filter(None, sam_lines) #removes list elements that are empty (i.e: '')

        new_modification_count = 0

        # Step 3:
        # go through each SAM line with the modifcation, grab the sequence, check the FQ file, and report back with new modificaiton count
        for sam_line in sam_lines:
            sam_line = sam_line.strip().split()
            anno_start, anno_end, anno_strand, species = sam_line[12], sam_line[13], sam_line[1], sam_line[14]
            read_start, read_end, read_nte5, read_nte3, mismatch_number, cigar, read_seq, miRNA_name, read_ID = sam_line[2], sam_line[3], sam_line[5], sam_line[6], sam_line[7], sam_line[8], sam_line[9], sam_line[10], sam_line[4]

            # modification index
            RNA_modification_index = get_read_index_position(i, read_start, read_end, anno_strand)

            if anno_strand is '-' or anno_strand is '16':
                read_seq = reverse_complement(read_seq)

            # add read_sequence to dictionary, key = RNA sequence, value = modification class object
            # dictionary for each read_sequence, this will be used to look up from sam files.
            # each list item is a different sam file read that had the modification in it.
            if read_seq not in check_fastq_seqs:
                modification_obj = modification_class(miRNA_name, i, RNA_modification_index) # miRNA name , modification profile ('chr:strand:position:ref>min'), RNA modification index
                check_fastq_seqs[read_seq] = [modification_obj]
            else:
                modification_obj = modification_class(miRNA_name, i, RNA_modification_index)
                check_fastq_seqs[read_seq].append(modification_obj)



    print 'done adding sequences to get qscore seq\nsize of check_fastq_seq dictionary = {}'.format(sys.getsizeof(check_fastq_seqs))

    return check_fastq_seqs

def fastq_phred_results_to_mirna(check_fastq_seqs):
    print 'condensing phred results --> miRNA results'
    '''
    check_fastq_seqs is a dictionary with the key = read_sequence, value = modification class object
    '''
    if check_fastq_seqs == '0':
        return '0'

    mirna_modification_dict = {}
    for i in check_fastq_seqs:
        modification_obj_list = check_fastq_seqs[i]
        for modification_obj in modification_obj_list:
            if modification_obj.modification_count is 0:
                print 'dropped modification {}'.format(modification_obj.modification_profile)
                continue
            if modification_obj.miRNA_name in mirna_modification_dict:
                mod_profile = modification_obj.modification_profile
                mod_profile = mod_profile.split(':')[2:]  # 'chr:strand:position:ref>alt'
                mod_profile.append(str(modification_obj.modification_count))
                updated_profile = ':'.join(mod_profile)  # 'position:ref>alt:count'
                mirna_modification_dict[modification_obj.miRNA_name].append(updated_profile)
            else:
                mod_profile = modification_obj.modification_profile
                mod_profile = mod_profile.split(':')[2:] # 'chr:strand:position:ref>alt' --> 'position:ref>alt:count'
                mod_profile.append(str(modification_obj.modification_count))
                updated_profile = ':'.join(mod_profile)
                mirna_modification_dict[modification_obj.miRNA_name] = [updated_profile]

    return mirna_modification_dict

def turn_mirna_dict_to_output(annotation_name, nte_profile_dict):


    if nte_profile_dict is '0':
        return '0', '0'
    if annotation_name not in nte_profile_dict:
        return '0', '0'

    else:
        unique_modification_dict = {} # key = modification_profile ('position:ref>alt', value = count

        phred_filter_nte_list = nte_profile_dict[annotation_name]
        # phred_filter_nte_string = '|'.join(phred_filter_nte_list)
        modification_type_count = 0
        # print phred_filter_nte_list
        # print phred_filter_nte_list
        for i in phred_filter_nte_list: # i = 'position:ref>alt:count'
            i =  i.split(':')  # [ position, ref>alt, count]
            profile = ':'.join(i[:2])  # 'position:ref>alt'

            count = int(i[2])
            # print 'i={}, profile = {} '.format(i, profile)
            if profile in unique_modification_dict:
                unique_modification_dict[profile] += count
            else:
                unique_modification_dict[profile] = count
            modification_type_count += count

        phred_filter_nte_list_string = []
        for profile in unique_modification_dict:
            value = ':'.join([ profile, str(unique_modification_dict[profile])])
            phred_filter_nte_list_string.append(value)

        phred_filter_nte_string = '|'.join(phred_filter_nte_list_string)

        return phred_filter_nte_string, str(modification_type_count)


def filter_modifications_from_phred(anno_summary, fastq_file, sam_file, phred_score_cutoff, phred64):

    # Step 1: Get modifications from the anno_summary file
    check_modifications_set = set()
    output_lines = []

    progress_list = progress(anno_summary)

    line_count = 0
    internal_modifcation_list = []
    nte3_modifcation_list = []
    nte5_modifcation_list= []
    with open(anno_summary, 'r') as anno_summary_obj:
        for i in anno_summary_obj:
            # skip header line
            if i.startswith('#') or i.startswith('annotation_name'):
                header = i.strip()
                output_lines.append(header)
                continue

            i = i.strip().split()
            annotation_name, species, chrom, strand, start, end, anno_sequence, total_read_count, perfect_match_count, anno_mismatch_count, rpm, total_nte5_count, total_nte3_count, total_internal_count, nte5, nte3, internal_mismatch = i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10], i[11], i[12], i[13], i[14], i[15], i[16]

            # if there are no annotated modifications in summary file, then add miRNA line to output and continue to next miRNA species
            if nte5 is '0' and nte3 is '0' and internal_mismatch is '0':
                continue

            # get lists of modifications specific for each modification type
            internal_modifcation_list = add_modification_set(internal_modifcation_list, internal_mismatch, chrom, strand)
            nte3_modifcation_list = add_modification_set(nte3_modifcation_list, nte3, chrom, strand)
            nte5_modifcation_list = add_modification_set(nte5_modifcation_list, nte5, chrom, strand)

    # print all_mir_modification_list
    print 'done collecting all miRNA modifications'
    # print '{} modifications collected'.format(len(all_mir_modification_list))

    # from the SAM file and list of modification profiles, get the RNA read sequences with the modifications in list. 'all_mir_modification_list'
    print 'starting internal modifications'
    # check_fastq_seqs will be created, this is a dictionary of the read sequences as keys, value = modification_class(miRNA_name, modification profile, RNA_modification_index). If a read shows up multiple times SAM file as a modificaiton, then it will add a new modification class object.
    check_fastq_seqs = get_read_sequences_from_annosummary(internal_modifcation_list, sam_file, fastq_file, phred_score_cutoff, "internal")
    # add counts of all modifications that pass
    check_fastq_seqs = filter_quality_score2(fastq_file, check_fastq_seqs, phred_score_cutoff)
    mirna_modification_dict_internal = fastq_phred_results_to_mirna(check_fastq_seqs)


    print 'starting nte3 modifications'
    check_fastq_seqs = get_read_sequences_from_annosummary(nte3_modifcation_list, sam_file, fastq_file, phred_score_cutoff, "nte3")
    check_fastq_seqs = filter_quality_score2(fastq_file, check_fastq_seqs, phred_score_cutoff)
    mirna_modification_dict_nte3 = fastq_phred_results_to_mirna(check_fastq_seqs)

    print 'starting nte5 modifications'
    check_fastq_seqs = get_read_sequences_from_annosummary(nte5_modifcation_list, sam_file, fastq_file, phred_score_cutoff, "nte5")
    check_fastq_seqs = filter_quality_score2(fastq_file, check_fastq_seqs, phred_score_cutoff)

    #dictionary has miRNA name as key and the specific type of modification as the value (it is in a list)
    mirna_modification_dict_nte5 = fastq_phred_results_to_mirna(check_fastq_seqs)

    output_lines = []
    with open(anno_summary, 'r') as anno_summary_obj:
        for i in anno_summary_obj:
            if i.startswith('#') or i.startswith('annotation_name'):
                header = i.strip()
                output_lines.append(header)
                continue

            i = i.strip().split()
            miRNA_name = i[0]
            i[16], i[13] = turn_mirna_dict_to_output(miRNA_name, mirna_modification_dict_internal)
            i[14], i[11]  = turn_mirna_dict_to_output(miRNA_name, mirna_modification_dict_nte5)
            i[15], i[12]  = turn_mirna_dict_to_output(miRNA_name, mirna_modification_dict_nte3)

            output_line = '\t'.join(i)
            output_lines.append(output_line)

    output_string = '\n'.join(output_lines)
    output_fn = anno_summary+'.phred{}'.format(phred_score_cutoff)
    with open(output_fn, 'w') as out_obj:
        out_obj.write(output_string)
        print 'output file: {}'.format(output_fn)

    return


def progress(read1_file):
    cmd = ['wc', '-l', read1_file]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    num_lines = int(stdout.strip().split()[0])

    percentages_time = []
    for i in range(0, 10):
        percentages_time.append(num_lines * i / 10)

    return percentages_time


def main():
    start_time = datetime.now()
    anno_summary, fastq_file, sam_file, phred_score_cutoff, phred64 = get_arguments()
    filter_modifications_from_phred(anno_summary, fastq_file, sam_file, phred_score_cutoff, phred64)

    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'
    return


if __name__ == "__main__":
    main()

# progress_list = progress(fin)
# progress_count = 0

# fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')
# for i in fin_obj:

#     i = i.strip().split()

#     # updating progress
#     if progress_count in progress_list:
#         percent_index = progress_list.index(progress_count)
#         percentage_val = (percent_index + 1) * 10
#         print '{}% completed'.format(percentage_val)
#     progress_count += 1
