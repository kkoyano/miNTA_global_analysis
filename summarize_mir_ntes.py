#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 04/25/2018

Description:
    summary script to get the %NTE

Input:

Output:

Notes:


'''
import argparse
import os
import sys
from subprocess import Popen, PIPE
import re
from itertools import permutations, izip, product
import gzip
import time
from datetime import datetime

# f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')


def get_arguments():

    parser = argparse.ArgumentParser(description='ADD DESCRIPTION')

    parser.add_argument('-i ', dest='input_file', type=str,
                        help='name of input <sample_name>.anno_summary file ')
    parser.add_argument('-o ', dest='output_directory', type=str,
                        help='name of output directory')

    args = parser.parse_args()

    fin = args.input_file
    fout = args.output_directory

    return fin, fout


def check_file_exist(file_path):
    if os.path.exists(file_path) is False:
        print 'ERROR!!!'
        print '{} does not exist'.format(file_path)
    return


def update_nte_dict(nte_string, nte_dict):
    # INPUT: nte_string = '101488429-101488430:CG>AA:45|101488429:C>A:14|101488429-101488430:CG>TT:6|101488429:C>T:5|101488429-101488430:CG>CA:4'
    # RETURN:  nte_dict = {key=NTE, value = COUNT}
    if nte_string is '0':
        return nte_dict

    nte_list = nte_string.split('|')  # ['101488429-101488430:CG>AA:45', '101488429:C>A:14', '101488429-101488430:CG>TT:6', '101488429:C>T:5', '101488429-101488430:CG>CA:4']
    for nte in nte_list:
        nte = nte.split(':')    # ['101488429', 'C>T', '5']
        count = int(nte[2])  # '5'
        nte_change = nte[1]  # 'C>T'
        nte_nucleotides = nte_change.split('>')[-1]  # 'C>T' --> ['C', 'T'] --> 'T'
        nte_nucleotides = re.sub('T', 'U', nte_nucleotides)
        if nte_nucleotides in nte_dict:
            nte_dict[nte_nucleotides] += count
        else:
            if 'N' in nte_nucleotides:
                nte_dict['Ns'] += count
            else:
                print 'NTE= {} not found in dict. check dictionary.'.format(nte_nucleotides)
    return nte_dict


def make_nte_dictionary(n=3, nucleotides='AUGC'):
    # the value of n is the length of the product iteration
    # the key is the NTE type (i.e 'AA', 'A', 'CC', 'GT'.... ) and the value is the count of the occurance
    # i.e {'A': 0, 'AA': 0, 'C': 0, 'GT': 0}
    d = {}
    # d['non_NTE_mm'] = 0

    for i in range(1, n + 1):
        x = product(nucleotides, repeat=i)
        for combination in x:
            combination = ''.join(combination)
            d[combination] = 0
    d['Ns'] = 0
    return d


def make_nte_dict_order(n=3, nucleotides='AUGC'):
    # the value of n is the length of the product iteration
    # the key is the order for output, and the value is the NTE type (i.e 'AA', 'A', 'CC', 'GT'.... )
    # i.e {1: 'A', 2: 'T', 3: 'G', 4: 'C', 5: 'N', 6: 'AA', 7: 'AT', 8: 'AG'}
    d = {}
    # d[0] = 'non_NTE_mm'
    count = 0
    for i in range(1, n + 1):
        x = product(nucleotides, repeat=i)
        for combination in x:
            combination = ''.join(combination)
            d[count] = combination
            count += 1
    d[count] = 'Ns'

    # d ={0: 'A', 1: 'T', 2: 'G', 3: 'C', 4: 'AA', 5: 'AT', 6: 'AG', ..., 84: 'Ns'}
    return d


def write_nte_output(nte_output_lines, output_file_obj, output_header_line):
    # nte_output_lines = <list>
    output_file_obj.write(output_header_line + '\n')
    for i in nte_output_lines:
        output_file_obj.write(i + '\n')

    output_file_obj.close()

    return


def count_ntes(fin, fout):
    '''
    Counts the number reads that are 3'NTE for a specific nucleotide combination
    fin = <sample_name>.NTE.sam.anno_overlap.anno_summary
    # annotation_name    species chromosome  strand  start   end anno_sequence   total_read_count    perfect_match_count anno_mismatch_count RPM nte5    nte3    internal_mismatch
    hsa-miR-379-5p  mirna   chr14   +   101488408   101488428   TGGTAGACTATGGAACGTAGG   140 11  129 9333.3  101488411:G>A:32|101488411:G>T:15|1014884010-1014884011:GT>AG:14    101488429-101488430:CG>AA:45|101488429:C>A:14|101488429-101488430:CG>TT:6|101488429:C>T:5|101488429-101488430:CG>CA:4   101488410:G>A:21|101488411:T>A:19

    '''
####### INITIALIZING STUFF ###########

    progress_list = progress(fin)
    progress_count = 0

    fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')  # Input files

    NUMBER_NTES = 3
    NUCLEOTIDES = 'AUGC'
    output_order = make_nte_dict_order(NUMBER_NTES, NUCLEOTIDES)  # making a dictionary of different nucleotide NTE combinations of up to 3 nucleotides, key = NTE, value = read count, initializes value to 0
    num_combination = len(output_order)  # NTE dictionary space

    fout_obj_nte3 = open('{}.nte3_counts'.format(fout), 'w')              # nte output files
    fout_obj_nte5 = open('{}.nte5_counts'.format(fout), 'w')              # nte output files
    fout_obj_internal = open('{}.internal_mm_counts'.format(fout), 'w')   # nte output files

    print '{}.nte3_counts'.format(fout)
    output_header_line = '\t'.join(["annotation_ID", "species", "total_read_count", "RPM", "total_nte5_count", "total_nte3_count", "total_internal_count"])   # HEADER

    for j in range(0, num_combination):  # d ={0: 'A', 1: 'T', 2: 'G', 3: 'C', 4: 'AA', 5: 'AT', 6: 'AG', ..., 84: 'Ns'}
        order_key = output_order[j]  # eg: when i = 0, order_key ='A'
        s = '\t' + str(order_key)
        output_header_line += s
    # print output_header_line
####### INITIALIZING STUFF ###########


#######    ACTUAL SCRIPT   ###########
    output_string3_list = []
    output_string5_list = []
    output_string_internal_list = []

    for i in fin_obj:
        if i.startswith('#') or i.startswith("annotation"):
            continue

        nte3_dict = make_nte_dictionary(NUMBER_NTES, NUCLEOTIDES)
        nte5_dict = make_nte_dictionary(NUMBER_NTES, NUCLEOTIDES)
        internal_modification_dict = make_nte_dictionary(NUMBER_NTES, NUCLEOTIDES)  # initialized dict, key = NTE, value = count. inital val = 0

        i = i.strip().split()
        # print i
        nte3_string = i[15]  # 101488429-101488430:CG>AA:45|101488429:C>A:14|101488429-101488430:CG>TT:6|101488429:C>T:5|101488429-101488430:CG>CA:4
        nte3_dict = update_nte_dict(nte3_string, nte3_dict)

        nte5_string = i[14]
        nte5_dict = update_nte_dict(nte5_string, nte5_dict)

        internal_modification_string = i[16]
        internal_modification_dict = update_nte_dict(internal_modification_string, internal_modification_dict)

        annotation_ID, species, total_read_count, RPM, total_nte5_count, total_nte3_count, total_internal_count = i[0], i[1], i[7], i[10], i[11], i[12], i[13]
        # print annotation_ID
        # print nte3_string
        # for key in nte3_dict:
        #     if nte3_dict[key] != 0:
        #         print '{}:{}'.format(key, nte3_dict[key])
        # writing output line
        beginning_line = '\t'.join([annotation_ID, species, total_read_count, RPM, total_nte5_count, total_nte3_count, total_internal_count])
        output_string3, output_string5, output_string_internal = beginning_line, beginning_line, beginning_line

        for j in range(0, num_combination):  # d ={0: 'A', 1: 'T', 2: 'G', 3: 'C', 4: 'AA', 5: 'AT', 6: 'AG', ..., 84: 'Ns'}
            order_key = output_order[j]  # eg: when i = 0, order_key ='A'
            output_string3 = output_string3 + '\t' + str(round(nte3_dict[order_key], 2))
            output_string5 = output_string5 + '\t' + str(round(nte5_dict[order_key], 2))
            output_string_internal = output_string_internal + '\t' + str(round(internal_modification_dict[order_key], 2))

        # print output_string3
        # sys.exit()
        # add values to their respective lists
        output_string3_list.append(output_string3)
        output_string5_list.append(output_string5)
        output_string_internal_list.append(output_string_internal)

        # GET THE OUTPUT TO PRINT OUT IN THE EXACT SAME FOR THE NTE DICTIONARY (BASE THIS OFF OF THE KEY VALUES)
        # updating progress
        if progress_count in progress_list:
            percent_index = progress_list.index(progress_count)
            percentage_val = (percent_index + 1) * 10
            print '{}% completed'.format(percentage_val)
        progress_count += 1

    write_nte_output(output_string3_list, fout_obj_nte3, output_header_line)
    write_nte_output(output_string5_list, fout_obj_nte5, output_header_line)
    write_nte_output(output_string_internal_list, fout_obj_internal, output_header_line)

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
    fin, fout = get_arguments()

    # make sure fout is a valid output directory

    if fout is None or fout is '.':
        fout = fin
    else:
        if not fout.endswith('/'):
            fout = fout + '/' + fin.split('/')[-1]
        else:
            fout = fout + fin.split('/')[-1]

    count_ntes(fin, fout)
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'
    return


if __name__ == "__main__":
    main()
