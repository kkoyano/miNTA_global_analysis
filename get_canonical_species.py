#!/bin/python

'''

Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 04/20/2018
Usage:
python get_canonical_species.py -i ./test/get_canonical_species_test/input.anno_overlap -t .5
correct output: ./test/get_canonical_species_test/input.anno_overlap.canonical.summary




Description:
- get_canonical_species.py
    -goal: get the canonical <miRNA/piRNA/snoRNA> species. For a single miRNA species,  will be defined as miRNA with the same 3' end position, can have different internal sequence. Also, does not count NTEs
    NOTE!!: -t 0.5 THRESHOLD= the most abundant end position needs to contribute to at least 50% of all reads. Otherwise we will not consider this miRNA species.

        ex:
        miR-1 reads             3'end result
        -----------------       no
        -------------------T    no
        ---------------T        no
        ----------------TT      yes
        ----------------        y
        --------X-------A       y
          --------------T       y
       -----------------        y
    --------------------        y
    ------X-------------        y

So we would just keep all the reads with the check marks because that is the most abundant ending position.

Input:
- <sample_name>.anno_overlap

Output:
- it would be convenient to pipe the results right into a file that is compatible with the count_mismatch_positions.py
-- should modify this script to take in both types of intput (NTE SAM FILE)  and also *.anno_overlap
----output: <sample_name>.anno_overlap
test.2mirnas.anno_overlap
chr8    -   22102487    22102508    251397|108  0   0   NM:i:0  MD:Z:22 TTCGCCCTCTCAACCCAGCTTT  hsa-miR-320a    chr8    22102488    22102509    mirna
chr8    -   22102485    22102508    251398|3    0   TT  NM:i:1  MD:Z:0T0T22 AATTCGCCCTCTCAACCCAGCTTT    hsa-miR-320a    chr8    22102488    22102509    mirna
chr8    -   22102486    22102507    251399|10   0   C   NM:i:1  MD:Z:0T21   GTTCGCCCTCTCAACCCAGCTT  hsa-miR-320a    chr8    22102488    22102509    mirna

Notes:
# also


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
# from helper_functions_mismatch_profile import reverse_complement

# f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')


class position_class:
    """ A class that holds the read coutn and the IDs"""

    def __init__(self, count, IDs):
        self.count = count
        self.IDs = [IDs]


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


class miRNA_annotation:
    """ A class that holds the information of the annotation, does not have to be a miRNA"""

    def __init__(self, chromosome, strand, total_reads):
        self.chromosome = chromosome
        self.strand = strand
        self.total_reads = total_reads

        self.end_position_dict = {}  # key = 'chr:strand:genomic_position1' , value = position class, .count, .IDs
        #{'chr:strand:genomic_position1':count, 'chr:strand:genomic_position2':count ...}
        self.nte3_count = 0
        self.nte5_count = 0
        self.canonical_end_count = 0
        self.number_end_positions = 0


def get_arguments():

    parser = argparse.ArgumentParser(description='ADD DESCRIPTION')

    parser.add_argument('-i ', dest='input_file', type=str,
                        help='name of input file <sample_name>.NTE.sam.anno_overlap')
    parser.add_argument('-s ', dest='summary_file', type=str,
                        help='name of *.NTE.sam.summary')
    parser.add_argument('-o ', dest='output_file', type=str,
                        help='name of output file')
    parser.add_argument('-t ', dest='threshold', type=float,
                        help='Number (float) between 0-1, threshold of reads to have canonical read end position, default= .5')

    # parser.add_argument('-e ', dest='exp', type=str,
    # help='a directory containing the experimental files from my predictions ')

    args = parser.parse_args()

    fin = args.input_file
    fout = args.output_file
    threshold = args.threshold
    summary_file = args.summary_file

    if fin is None or summary_file is None:
        parser.print_help()
        exit()
    if threshold is None:
        threshold = .5
    if fout is None:
        fout = ''
    return fin, fout, summary_file, threshold


def get_total_reads(summary_file):
    '''
    Summary file to get the number of total mapped reads

    '''
    library_depth = 0
    with open(summary_file, 'r') as summary_fn:
        for i in summary_fn:
            if i.startswith("total number of reads:"):
                i = i.strip().split('\t')
                library_depth = i[-1]
                return library_depth
                # break


def check_file_exist(file_path):
    if os.path.exists(file_path) is False:
        print 'ERROR!!!'
        print '{} does not exist'.format(file_path)
    return


def find_end_position(anno_chrom, anno_strand, read_start, read_end, read_nte3):
    # position = 'chrom:strand:genomic_position'
    end_position = ''

    if read_nte3 == '0':
        read_nte3 = ''

    # if the read has a 3' NTE and is on the positive strand, shift the end position back 1nt
    if anno_strand == "+":
        if len(read_nte3) > 0:
            end_position = int(read_end) - len(read_nte3)
        else:
            end_position = int(read_end)

    # if the read has a 3' NTE and is on the negative strand, shift the end position up 1nt from the start site (on the negative strand, start site of RNA is really the 3' end position)
    else:
        if len(read_nte3) > 0:
            end_position = int(read_start) + len(read_nte3)
        else:
            end_position = int(read_start)
    position_key = ':'.join([anno_chrom, anno_strand, str(end_position)])
    return position_key


def get_canonical_species(fin, fout, threshold, library_depth):

    progress_list = progress(fin)
    progress_count = 0
    annotation_dict = {}
    fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')
    for i in fin_obj:
        if i.startswith('#'):
            continue
        i = i.strip().split()

        # get column values
        anno_name, anno_chrom, anno_start, anno_end, anno_strand, species = i[10], i[11], i[12], i[13], i[1], i[14]
        read_strand, read_start, read_end, readID, read_nte5, read_nte3, mismatch_number, cigar, read_seq = i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]
        read_count = int(readID.split('|')[1])

        position = find_end_position(anno_chrom, anno_strand, read_start, read_end, read_nte3)
        # print '{} \t {} '.format(readID, position)
        # position = 'chrom:strand:genomic_position'

        # annotation_dict = {'mirna_name': position_dict, 'mirna2_name': position_dict }
        # position_dict = {'chr:strand:position': <position_class_obj> }
        # <position_class_obj> =
        #       position_class_obj.count #int
        #       position_class_obj.IDs   #list

        if anno_name in annotation_dict:
            anno_name_class = annotation_dict[anno_name]
            anno_name_class.total_reads += read_count
            # print anno_name_class.total_reads
            position_dict = anno_name_class.end_position_dict
            if position in position_dict:  # if position is already in the dictionary, add to the canonical read count
                position_dict[position].count += read_count
                position_dict[position].IDs.append(readID)

            else:
                position_class_object = position_class(read_count, readID)
                position_dict[position] = position_class_object
            annotation_dict[anno_name].end_position_dict = position_dict

        else:
            annotation_dict[anno_name] = miRNA_annotation(anno_chrom, anno_strand, read_count)
            position_class_object = position_class(read_count, readID)
            annotation_dict[anno_name].end_position_dict = {position: position_class_object}

        if progress_count in progress_list:
            percent_index = progress_list.index(progress_count)
            percentage_val = (percent_index + 1) * 10
            # print '{}% completed'.format(percentage_val)
        progress_count += 1
    fin_obj.close()

    # pick the canonical sequence
    # canonical_annotation_dict = {}  # this is a dictionary with all the canonical end position for each miRNA, key = miRNA_ID, value = set of IDs
    canonical_end_IDs_set = set()

#################################
# UPDATE: add user defined threshold of % of reads needing to have the canonical end position (default = 50% total reads)
# record:
    # 1. how many reads are lost
    # 2. How many miRNAs are not included
    # 3. How many 3'NTEs are lost in total
#
##################################
    output_summary_fn = fin + '.t' + str(threshold) + '.canonical_end.summary'

    if fout is not '' or fout is not '.':  # fout is a directory that directs where the output file should be written out to
        if fout.endswith('/'):
            output_summary_fn = fout + fin + '.t' + str(threshold) + '.canonical_end.summary'
        else:
            output_summary_fn = fout + '/' + fin + '.t' + str(threshold) + '.canonical_end.summary'

    output_summary_obj = open(output_summary_fn, 'w')
    header_line = '\t'.join(['annotation_ID', 'total_read_count', 'canonical_end_count', 'percent_canonical_end_count', 'library_depth', 'RPM_total_reads', 'RPM_canonical', 'number_end_positions', 'canonical_end_position', 'decision'])
    output_summary_obj.write(header_line + '\n')  # write header

    for annotation_key in annotation_dict:
        anno_class_obj = annotation_dict[annotation_key]
        position_dict = anno_class_obj.end_position_dict

        number_end_locations = len(position_dict)
        key_max = max(position_dict.keys(), key=(lambda k: position_dict[k].count))  # this gets the key with the highest
        max_position_obj = position_dict[key_max]
        anno_class_obj.canonical_end_count = max_position_obj.count
        percent_canonical_end_reads = max_position_obj.count / float(anno_class_obj.total_reads)
        decision = 'REMOVE'
        if percent_canonical_end_reads >= threshold:  # TRESHOLD = 0.5, if the canonical end is >= 50% then keep, otherwise don't count

            # canonical_annotation_dict[annotation_key] = canonical_end_pos
            canonical_end_IDs_set.update(max_position_obj.IDs)  # if it passes the threshold then add IDs to the canonical_ids set, this will be used to output only reads that are found in this set, (i.e annotation species that have at least 50% of their reads being the canonical)
            decision = 'KEEP'
        percent_canonical_end_reads = str(round(percent_canonical_end_reads, 2) * 100)
        output_line = '\t'.join([annotation_key, str(anno_class_obj.total_reads), str(max_position_obj.count),
                                 str(percent_canonical_end_reads), str(library_depth),
                                 str(round((anno_class_obj.total_reads * library_depth / 1000000), 2)),
                                 str(round((max_position_obj.count * library_depth / 1000000), 2)),
                                 str(number_end_locations), key_max, decision])
        output_summary_obj.write(output_line + '\n')
        # print output_line

    # output the lines that have the IDs with the canonical end positions,
    fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')

    output_filename = fin + '.t' + str(threshold) + '.canonical_end'

    if fout is not '' or fout is not '.':  # fout is a directory that directs where the output file should be written out to
        if fout.endswith('/'):
            output_filename = fout + fin + '.t' + str(threshold) + '.canonical_end'
        else:
            output_filename = fout + '/' + fin + '.t' + str(threshold) + '.canonical_end'

    output_obj = open(output_filename, 'w')
    for i in fin_obj:
        if i.startswith('#'):
            continue
        original_line = i
        i = i.strip().split()
        readID = i[4]
        if readID in canonical_end_IDs_set:
            output_obj.write(original_line)
    output_obj.close()

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
    fin, fout, summary_file, threshold = get_arguments()
    library_depth = get_total_reads(summary_file)
    print library_depth

    get_canonical_species(fin, fout, threshold, float(library_depth))
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'

    return


if __name__ == "__main__":
    main()
