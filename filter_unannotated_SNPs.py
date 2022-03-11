#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 06/20/2018

Description: filter out all modifications that maybe SNPs
Homozygous SNPs: modification occurs in 100% of all the reads
Heterozygoud SNPs: modification occurs in 50% of all the reads (option to put buffer of 45-55%)

Input:
file to check for homozygous or heterozygous snps

Usage:
filter_unannotated_SNPs.py -i <sample_name>.anno_summary #note the filename extension does not have the be anno_summary but it need to be in anno_summary format (beinf that rows are miRNA species)


Output:

*.snp2
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


def get_arguments():

    parser = argparse.ArgumentParser(description='ADD DESCRIPTION')

    parser.add_argument('-i ', dest='input_file', type=str,
                        help='name of input file')
    parser.add_argument('-m ', dest='minimum_het', type=int, default=50,
                        help='lower threshold to consider modification a heterozygous SNV (Default= 50)')
    parser.add_argument('-M', dest='Maximum_het', type=int,default=50,
                        help='upper threshold to consider modification a heterozygous SNV (Default= 50)')


    # parser.add_argument('-e ', dest='exp', type=str,
    # help='a directory containing the experimental files from my predictions ')

    args = parser.parse_args()
    # if len(sys.argv) < 2:
    # parser.print_help()
    # exit()
    fin = args.input_file
    minimum_het = args.minimum_het
    Maximum_het = args.Maximum_het


    return fin, minimum_het, Maximum_het


def check_file_exist(file_path):
    if os.path.exists(file_path) == False:
        print 'ERROR!!!'
        print '{} does not exist'.format(file_path)
    return

def check_if_homo_het_snp(mir_name , mir_total_count, nte_count, nte_profile, homozygous_count, heterozygous_count, minimum_het, Maximum_het):
    # check if any of the modifications appear to be homozygous SNPs (in 100% of reads) or heterozygous SNPs (in 50% of total reads)

    if nte_profile == '0':
        return nte_count, '0', homozygous_count, heterozygous_count

    nte_profile_list = nte_profile.split('|')

    return_modifications = []
    for nte_modification in nte_profile_list:
        modification_count = int(nte_modification.split(':')[2])

        # check if modification is homozygous, ALL of hte miRNA species reads have this modification
        if modification_count == mir_total_count:
            homozygous_count += 1
            nte_count -= modification_count
            # print '{} - homozygous! '.format(mir_name)

        # check if modification is heterozygous, level of modification is between min and max
        elif  minimum_het <= (modification_count*100/mir_total_count) <= Maximum_het:
            heterozygous_count += 1
            nte_count -= modification_count
            # print '{} - heterozygous! '.format(mir_name)
        else:
            return_modifications.append(nte_modification)

    if len(return_modifications) is 0:
        return_modifications = '0'
    return_modifications = '|'.join(return_modifications)

    return nte_count, return_modifications, homozygous_count, heterozygous_count

def filter_unannotated_snps(fin, minimum_het, Maximum_het):

    progress_list = progress(fin)
    progress_count = 0

    header = ''
    fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')
    homozygous_count, heterozygous_count = 0, 0
    all_lines = []
    for i in fin_obj:
        if i.startswith('annotation_name'):
            header = i
            continue
        i = i.strip().split()
        mirna_species, total_reads, total_nte5_count, total_nte3_count, total_internal_count, nte5, nte3, internal_mismatch = i[0], int(i[7]), int(i[11]), int(i[12]), int(i[13]), i[14],  i[15], i[16]

        total_nte5_count, nte5, homozygous_count, heterozygous_count= check_if_homo_het_snp(mirna_species, total_reads,total_nte5_count, nte5, homozygous_count, heterozygous_count, minimum_het, Maximum_het)
        total_nte3_count, nte3, homozygous_count, heterozygous_count= check_if_homo_het_snp(mirna_species, total_reads,total_nte3_count, nte3, homozygous_count, heterozygous_count, minimum_het, Maximum_het)
        total_internal_count, internal_mismatch, homozygous_count, heterozygous_count= check_if_homo_het_snp(mirna_species, total_reads,total_internal_count, internal_mismatch, homozygous_count, heterozygous_count, minimum_het, Maximum_het)

        i[0], i[7], i[11], i[12], i[13], i[14],  i[15], i[16] = mirna_species, total_reads, total_nte5_count, total_nte3_count, total_internal_count, nte5, nte3, internal_mismatch

        i = map(str, i)
        line = '\t'.join(i)
        all_lines.append(line)


    output ='\n'.join(all_lines)
    file_ext = '.snv-m{}M{}'.format(minimum_het, Maximum_het)
    with open(fin+ file_ext, 'w') as out_obj:
        out_obj.write(header)
        out_obj.write(output)


        # updating progress
        # if progress_count in progress_list:
        #     percent_index = progress_list.index(progress_count)
        #     percentage_val = (percent_index + 1) * 10
        #     print '{}% completed'.format(percentage_val)
        # progress_count += 1

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
    fin, minimum_het, Maximum_het = get_arguments()
    filter_unannotated_snps(fin, minimum_het, Maximum_het)
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'
    return


if __name__ == "__main__":
    main()
