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
    parser.add_argument('-b ', dest='bed_file', type=str,
                        help='bed file that is 0-based and the best scores from psl file')


    args = parser.parse_args()
    # if len(sys.argv) < 2:
    # parser.print_help()
    # exit()
    fin_sam = args.sam_file

    fin_bed = args.bed_file

    return fin_sam, fin_bed


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

def compare_coordiantes(fin_sam, fin_bed):

    # progress_list = progress(fin_bed)
    # progress_count = 0
    fout_noblat= fin_bed+".results"
    fout_noblat_obj = open(fout_noblat, 'w')
    fin_obj = gzip.open(fin_bed, 'rb') if re.search(r'\.gz$', fin_bed) else open(fin_bed, 'r')
    # parse bed file for coordinates
    bed_file_dict = {}  #key = readID,  value = miRNA class obj


    # opening bed file of best BLAT alignments
    for i in fin_obj:
        if i.startswith("#"): continue

        i = i.strip().split()
        chrom, Tstart, Tend, readID, Qsize, strand, Qstart, Qend  = i[0], int(i[1]), int(i[2]), i[3], int(i[4]), i[7], int(i[8]), int(i[9])
        Tstart_updated = (Tstart - Qstart) + 1
        Tend_updated  = Tend + (Qsize - Qend)
        # add to dictionary
        if readID not in bed_file_dict:
            miRNA_obj = miRNAclass(chrom, str(Tstart_updated), str(Tend_updated), readID, str(Qsize), strand, str(Qstart), str(Qend))
            bed_file_dict[readID] = miRNA_obj
        else:
            miRNA_obj = bed_file_dict[readID]
            miRNA_obj.chrom.append(chrom)
            miRNA_obj.tstart.append(Tstart_updated)
            miRNA_obj.tend.append(Tend_updated)
            miRNA_obj.strand.append(strand)
            miRNA_obj.Qstart.append(Qstart)
            miRNA_obj.Qend.append(Qend)


    fin_obj.close()
    fin_sam_obj = gzip.open(fin_sam, 'rb') if re.search(r'\.gz$', fin_sam) else open(fin_sam, 'r')
    # print bed_file_dict
    correct_alignments =0
    correct_alignment_mmap = 0
    diff_alignment_and_mmap = 0
    diff_and_uniq = 0
    no_blat_result = 0
    perfect =0
    mismatch_count = 0

    correct_alignments_readcount =0
    correct_alignment_mmap_readcount = 0
    diff_alignment_and_mmap_readcount = 0
    diff_and_uniq_readcount = 0
    no_blat_result_readcount = 0
    perfect_readcount =0
    mismatch_count_readcount = 0

    # opening sam file with original bowtie alignments
    header = ''
    for i in fin_sam_obj:
        if i.startswith("#chr"): #this is the header line
            i = i.strip()
            fout_noblat_obj.write(i+"\tCONSENSUS\n")
            continue
        i = i.strip().split()
        chrom, strand, start, end, readID, nm_flag, ncRNA_species = i[0], i[1], i[2], i[3], i[4], i[7], i[14]
        readID_count = int(readID.split('|')[1])
        if ncRNA_species == 'mirna':
            if nm_flag == "NM:i:1":
                mismatch_count +=1
                mismatch_count_readcount += readID_count
            else:
                perfect += 1
                perfect_readcount += readID_count

            if readID in bed_file_dict:

                blat_obj = bed_file_dict[readID]

                result = compare_coords(blat_obj, chrom, strand, start, end) #will return true or false
                if result is True and len(blat_obj.tstart) > 1:
                    # print 'Same, but non-uniquely mapped- {}'.format(readID)
                    # print 'same-multimap'
                    i.append('SAME-MULTIMAPP')
                    correct_alignment_mmap += 1
                    correct_alignment_mmap_readcount += readID_count
                    line = '\t'.join(i)
                    fout_noblat_obj.write(line + '\n')
                elif result is False:
                    if len(blat_obj.tstart) > 1:
                        i.append('DIFF-MULTIMAPP')
                        diff_alignment_and_mmap += 1
                        diff_alignment_and_mmap_readcount += readID_count
                        # mmap_positions = []
                        # for j in range(len(blat_obj.tstart)):
                        #     alt_coord = str(blat_obj.tstart[j]) + '-' + str(blat_obj.tend[j])
                        #     alt_location = ':'.join([blat_obj.chrom[j],blat_obj.strand[j], alt_coord])
                        #     mmap_positions.append(alt_location)
                        # mmap_string = '|'.join(mmap_positions)
                        # i.append(mmap_string)
                        line = '\t'.join(i)
                        fout_noblat_obj.write(line + '\n')
                    else: # the result is uniquely mapped, but mapped to a different location
                        i.append('DIFF-UNIQ')
                        diff_and_uniq += 1
                        diff_and_uniq_readcount += readID_count
                        alt_coord = blat_obj.tstart[0] + '-' + blat_obj.tend[0]
                        alt_location = ':'.join([blat_obj.chrom[0],blat_obj.strand[0], alt_coord])
                        # print alt_location
                        i.append(alt_location)
                        # print i
                        line = '\t'.join(i)
                        fout_noblat_obj.write(line + '\n')
                else: # the result is true and it only has 1 mapping location (uniq)
                    i.append('CORRECT')
                    line = '\t'.join(i)
                    fout_noblat_obj.write(line + '\n')
                    correct_alignments += 1
                    correct_alignments_readcount += readID_count
            else:
                i.append('NO_BLAT_RESULT')
                line = '\t'.join(i)
                no_blat_result += 1
                no_blat_result_readcount += readID_count
                fout_noblat_obj.write(line + '\n')


                continue #meaning that the read is perfectly mapped
        else: continue
    fout_noblat_obj.close()
    fout_summary= fin_bed+".summary"
    with open(fout_summary, 'w') as fout_summary_obj:
        fout_summary_obj.write("Result\tmapping\tuniq_seq_count\tuniq_seq_count(%)\tread_count\tread_count(%)\tperfect_sequences\tperfect_readcount\tmismatch_sequences\tmismatch_readcount\n")

        fout_summary_obj.write('Correct\tuniq:\t{}\t{}%\t{}\t{}%\t{}\t{}\t{}\t{}\n'.format(correct_alignments, round(correct_alignments*100.0/mismatch_count, 2), correct_alignments_readcount, round(correct_alignments_readcount*100.0/mismatch_count_readcount, 2), perfect, perfect_readcount, mismatch_count, mismatch_count_readcount))
        fout_summary_obj.write('Correct\tmultimapped:\t{}\t{}%\t{}\t{}%\t{}\t{}\t{}\t{}\n'.format(correct_alignment_mmap,round(correct_alignment_mmap*100.0/mismatch_count, 2), correct_alignment_mmap_readcount, round(correct_alignment_mmap_readcount*100.0/mismatch_count_readcount, 2), perfect, perfect_readcount, mismatch_count, mismatch_count_readcount))
        fout_summary_obj.write('Different\tmultimapped:\t{}\t{}%\t{}\t{}%\t{}\t{}\t{}\t{}\n'.format(diff_alignment_and_mmap,round(diff_alignment_and_mmap*100.0/mismatch_count, 2), diff_alignment_and_mmap_readcount, round(diff_alignment_and_mmap_readcount*100.0/mismatch_count_readcount, 2), perfect, perfect_readcount, mismatch_count, mismatch_count_readcount))
        fout_summary_obj.write('Different\tuniq:\t{}\t{}%\t{}\t{}%\t{}\t{}\t{}\t{}\n'.format(diff_and_uniq,round(diff_and_uniq*100.0/mismatch_count, 2), diff_and_uniq_readcount, round(diff_and_uniq_readcount*100.0/mismatch_count_readcount, 2), perfect, perfect_readcount, mismatch_count, mismatch_count_readcount))
        fout_summary_obj.write('None\tNone:\t{}\t{}%\t{}\t{}%\t{}\t{}\t{}\t{}\n'.format(no_blat_result,round(no_blat_result*100.0/mismatch_count, 2), no_blat_result_readcount, round(no_blat_result_readcount*100.0/mismatch_count_readcount, 2), perfect, perfect_readcount, mismatch_count, mismatch_count_readcount))


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
    fin_sam, fin_bed= get_arguments()


    # TESTING #
    # fin_sam = '/Users/kiku/mount/nte_global_analysis_2017-09/unique_reads/data/u87_inhouse/canonical_reads/U87_sRNA_R1_S1_L002_R1_001.NTE.sam.anno_overlap.t0.5.canonical_end'
    # fin_bed = '/Users/kiku/mount/nte_global_analysis_2017-09/unique_reads/data/u87_inhouse/canonical_reads/blat/U87_sRNA_R1_S1_L002_R1_001.NTE.sam.anno_overlap.t0.5.canonical_end.mirna.ALL.psl.best_score.bed'
    # fout = '/Users/kiku/mount/test/sam2bed.txt'
    compare_coordiantes(fin_sam, fin_bed)

    #NEXT
    # MODIFY OUTPUT SO THAT IT PRINTS OUT THE QSTART, QEND, AND # MATCHES, # MISMATCHES S
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'
    return


if __name__ == "__main__":
    main()
