#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 04/16/2018

Description: script to filter out snps from NTE summary files


Input:
- NTE summary file
- SNP list or SNP directory

Output:
- updated NTE summary file

Notes:

Usage: python filter_nte_snps.py -i <*.NTE.sam> -snp_file /Users/kiku/mount/annotations/hg19/snps/hg19.comprehensive.7.14.2017.snp.mirna.txt -boundary 10

-input_file
-snp_file
/u/home/k/kkoyano/nobackup-gxxiao3/annotations/hg19/snps/hg19.comprehensive.7.14.2017.snp.mirna.txt
    chr1 25245882 C>A rs766263437|2016AD|dbSNP147
    chr1 156905955 T>C |gnomAD-exomes
    chr1 226109831 C>T rs751426926|2016AD|dbSNP147
    chr1 1104455 G>A rs368678282|2016AD|dbSNP147
    chr1 231155599 C>A rs750328882|2016AD|dbSNP147
    chr1 12639050 G>A rs182935655|2016AD|dbSNP147
    chr1 19209751 G>A rs773230705|2016AD|dbSNP147

-SNP_DIR ex:
/u/home/k/kkoyano/nobackup-gxxiao3/annotations/hg19/snps/hg19.comprehensive.7.14.2017/hg19.comprehensive.7.14.2017.snp.


pseudocode:

1. take in SNP file or directory and make snp set chr:position:ref>minor
2. open the <sample_name>.anno_overlap.summary file
    - take in the 3' end modification line
        - split all modifications up,
        - work from end of RNA to beginning, if SNP occurs, STOP, remove NTE.
        - modify NTE count, leave all others the same.
    - take in the 5' end modification line
        - split all modifications up,
        - work from beginning of RNA, if SNP occurs, STOP, remove NTE.
        - modify NTE count, leave all others the same.
    - take in the internal modification line
        - split all modifications up,
        - if SNP occurs, STOP, remove
        - modify NTE count, leave all others the same.


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
from helper_functions_mismatch_profile import reverse_complement, progress, parse_nte3_list, parse_nte5_list, check_snp_position, RC_ref_minor, check_file_exist

# f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')


def get_arguments():

    parser = argparse.ArgumentParser(description='Usage:\npython filter_nte_snps.py -i <*.NTE.sam> -snp_file /Users/kiku/mount/annotations/hg19/snps/hg19.comprehensive.7.14.2017.snp.mirna.txt -boundary 10')

    parser.add_argument('-anno_summary ', dest='anno_summary', type=str,
                        help='name of anno_summary file')
    parser.add_argument('-o ', dest='output_file', type=str,
                        help='name of output file')
    parser.add_argument('-snp_file ', dest='snp_file', type=str,
                        help='name of snp file containing all snps, 4cols: (1)chromosome (2)position (3)ref>alt (4)source')
    parser.add_argument('-SNP_DIR ', dest='SNP_DIR', type=str,
                        help='Prefix of the SNP directory that contains the files of all SNPs')
    # parser.add_argument('-e ', dest='exp', type=str,
    # help='a directory containing the experimental files from my predictions ')

    args = parser.parse_args()

    anno_summary = args.anno_summary
    output_file = args.output_file
    snp_file = args.snp_file
    SNP_DIR = args.SNP_DIR

    if anno_summary is None:
        print 'Need to add an anno_summary file with using the -anno_summary flag'
        sys.exit()
    # return anno_summary, output_file, snp_file, SNP_DIR
    return anno_summary, snp_file, SNP_DIR, output_file


def get_snp_positions(chromosome, snp_dir, is_file_bool):
    '''
    make the snp_set for given chromosome or file.
    keys in set  --> 'chr:genomic_position:ref>minor_allele'
    return --> a set of snps from given file or directory
    '''
    # snp_dir = '/u/nobackup/gxxiao3/kkoyano/annotations/hg19/snps/hg19.comprehensive.7.14.2017/hg19.comprehensive.7.14.2017.snp.'
    # file = "hg19.comprehensive.7.14.2017.snp.10.gz"

    snp_set = set()

    filename = ''

    # checking to see if snp_dir is a file or not
    if os.path.isfile(snp_dir) or is_file_bool:
        filename = snp_dir
    else:
        chrom = re.sub("chr", "", chromosome)
        filename = snp_dir + chrom + ".gz"
    print 'opening snp file {}'.format(filename)

    if not os.path.isfile(filename):
        print 'no snp file for {}'.format(filename)
        return snp_set

    f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')

    # look through all the SNPs in the file
    for j in f:
        j = j.strip().split()  # chr1 1104455 G>A rs368678282|2016AD|dbSNP147
        # j.pop()  # remove the source --> ["chr1", "1104455", "G>A "]
        key = ':'.join(j[:3])  # key = 'chr19:29194762:G>A' --> 'chromosome:position:ref>minor_allele'
        snp_set.add(key)
    print 'finished making SNP set for {}'.format(filename)
    return snp_set  # this is the snp set of this directory or the file


def filter_nte_snps(fin, snp_dir, snp_file_bool, output_file):
    """
    #annotation_name    species chromosome  strand  start   end anno_sequence   total_read_count    perfect_match_count anno_mismatch_count RPM nte5    nte3    internal_mismatch

    hsa-miR-379-5p  mirna   chr14   +   101488408   101488428   TGGTAGACTATGGAACGTAGG   140 11  129 9333.3  101488411:G>A:32|101488411:G>T:15|1014884010-1014884011:GT>AG:14    101488429-101488430:CG>AA:45|101488429:C>A:14|101488429-101488430:CG>TT:6|101488429:C>T:5|101488429-101488430:CG>CA:4   101488410:G>A:21|101488411:T>A:19
    hsa-miR-320a    mirna   chr8    -   22102488    22102509    AAAAGCTGGGTTGAGAGGGCGA  310 108 202 20666.7 22102508:A>T:11|22102508-22102509:AA>TT:7|22102508-22102509:AA>TG:9 22102486-22102485:AA>TT:3|22102487:A>T:23|22102486-22102485:AA>TC:40|22102486:A>T:11|22102486:A>C:10    22102492:G>T:30|22102490:G>C:108

    """

    progress_list = progress(fin)

    snp_set = set()  # chr:strand:position:minor

    # user given snp_file so just use the snps given in this specified file.
    if snp_file_bool:
        temp_chr = ''
        snp_set = get_snp_positions(temp_chr, snp_dir, snp_file_bool)

    header = ''
    output_lines = []
    # print snp_set
    with open(fin, 'r') as fin_obj:
        prev_chom = ''
        progress_count = 0

        for i in fin_obj:
            position_added = False  # used to only include reads that were counted in the dictionary (because we don't count the other reads that don't fall in the dictionaries )
            if i.startswith("annotation_name") or i.startswith("#"):
                header = i
                continue
            i = i.strip().split('\t')

            annotation_name, species, anno_chrom, anno_strand, total_read_count, perfect_match_count, anno_mismatch_count, RPM = i[0], i[1], i[2], i[3], i[7], i[8], i[9], i[10]
            nte5_count, nte3_count, internal_mm_count = i[11], i[12], i[13]
            read_nte5, read_nte3, internal_mm = i[14], i[15], i[16]

            read_nte3_nosnp, update_nte3_count = parse_nte3_list(read_nte3, snp_set, anno_chrom, anno_strand, nte3_count)
            read_nte5_nosnp, update_nte5_count = parse_nte5_list(read_nte5, snp_set, anno_chrom, anno_strand, nte5_count)
            read_internal_mm_nosnp, update_internal_mm_count = parse_nte3_list(internal_mm, snp_set, anno_chrom, anno_strand, internal_mm_count)

            # update annotation line so that there are no SNPs

            i[11], i[12], i[13], i[14], i[15], i[16] = update_nte5_count, update_nte3_count, update_internal_mm_count, read_nte5_nosnp, read_nte3_nosnp, read_internal_mm_nosnp

            output_lines.append(i)
            # if the chromosome does not have an associated SNP file, then do not upload SNP
            if prev_chom != anno_chrom and not snp_file_bool:  # if new chromosome, then need to get new snp_list for specific chromosome
                snp_set = get_snp_positions(anno_chrom, snp_dir, snp_file_bool)
                prev_chom = anno_chrom

            if progress_count in progress_list:
                percent_index = progress_list.index(progress_count)
                percentage_val = (percent_index + 1) * 10
                print '{}% completed'.format(percentage_val)

            progress_count += 1
            if position_added:  # initially false, if a position was added, then add read_counts to mismatches_read_count
                mismatch_read_count += read_count

    # output_fn = fin + '.nosnp'
    if output_file is None:
        output_file = fin + '.nosnp'
    outfn_obj = open(output_file, 'w')
    outfn_obj.write(header)
    for i in output_lines:
        i = '\t'.join(i)
        outfn_obj.write(i + '\n')

    outfn_obj.close()
    return


def main():
    # Testing Section

    # REAL
    print 'Starting script...'
    start_time = datetime.now()
    # anno_summary, output_file, snp_file, SNP_DIR = get_arguments()
    anno_summary, snp_file, SNP_DIR, output_file = get_arguments()
    snp_file_bool = True
    if snp_file is None:  # if snp_file is not provided
        snp_file_bool = False
        snp_file = SNP_DIR  # check to see if dirctory then

    else:
        check_file_exist(snp_file)  # check to see if file exists
        snp_file_bool = True

    check_file_exist(anno_summary)
    filter_nte_snps(anno_summary, snp_file, snp_file_bool, output_file)
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print '\nDONE!'

    return


if __name__ == "__main__":
    main()
