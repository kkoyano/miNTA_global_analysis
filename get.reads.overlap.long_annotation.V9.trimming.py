#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 03/20/2017

Version 1: Once a read alignment matches an annotation, then choose that annotation.
For multimapped reads this will assign all the reads to the first annotation it finds, no prioritization in the annotation species (i.e miRNA, piRNA, etc. )

Version 2: Prioritizes the annotation species from a given list. Must go through all read's alignments first to decide annotation. Then based off of the prioritization
given by a file after the flag -PriorityFile, once read is aligned to an annotation, it does not read subsequent alignments for mmap annotations,
-if it is mmap, it assigns all reads to first annotation that read is assigned.

Version 3:
- imposes a 1/n distribution of reads for mmap reads. Where each possible mapping location recieves 1/n reads, where 'n' is the number of possible mapped locations
- NOTE: PARSER STILL ASSUMES --BEST --STRATA. if using a diff alignment, then must pre-filter sam file.
- extra column of the read sequence is added, this sequence from the sequencing read, not the annotation.

Version 4: long anntoation (for annotations where the molecule is longer than the read length)
*NEW* - the read must be containg within the annotaiton
- imposes a 1/n distribution of reads for mmap reads. Where each possible mapping location recieves 1/n reads, where 'n' is the number of possible mapped locations
- NOTE: PARSER STILL ASSUMES --BEST --STRATA. if using a diff alignment, then must pre-filter sam file.
- extra column of the read sequence is added, this sequence from the sequencing read, not the annotation.
- no priority file, if an alignment matches to multiple places, then allocate reads in both annotations

Version 6: for trimming sam files, same as version 4, except it will output the NTE fields that were added in the trimming analysis.
**note : assumes the first line is the header so the program just reads in the first line of the sam file

Version 7: get.reads.overlap.long_annotation.V7.trimming.py
*NEW* priority for coordinates assigned to more than 2 annotations
-takes in both normal sam file and also sam file produced from trimming pipeline (fields included for 3NTE and 5NTE )
- annotations are within -3 of start coordinate and +3 of end coordinate

Version 8: get.reads.overlap.long_annotation.V8.trimming.py
-takes in both normal sam file and also sam file produced from trimming pipeline (fields included for 3NTE and 5NTE )

*NEW*
    -annotations of species that are greater than 50nt will be subjected to the constraint of having the start and end alignment
    being within +3,-3 of the start and end annotation coordinate
    - If annotation is larger than 50nt, I just allow the alignment location to fall within the annotation coordinate,
    REASONING: This is because for longer annotations (ex: a 300nt ncRNA molecule), we expect that the read would be a fragment of this
    However, for small RNA (like a 34nt piRNA), if this truely is a piRNA molecule, we expect it to be close to the reported full length.
    if there is a 15nt read that maps to a piRNA location, we don't expect this to actually be the annotated piRNA molecule because it isn't the full length.
    Usage :python get.reads.overlap.long_annotation.V8.trimming.py -s overlap.test.NTE.sam -a annotation_paths_mir_pir_trf_rfam.txt -PriorityAnno /u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/scripts/priority_file.txt --sam

Version 9: get.reads.overlap.long_annotation.V8.trimming.py
- will add 3 fields, 3'NTE, 5'NTE and internal, these fields contain the modification profile for the read's alignment
- profile will be 'chr:strand:position:ref>alt'


Defintions:
classes: different classes of RNA or noncoding species. i.e mirna, pirna, tRf, snoRNA
species: names various species within each class of RNA, i.e let-7a, piR-12352, sno-12461

Usage:
python get.overlap.long_annotation.V6.trimming.py -s sam_filenames.txt -a ann_filenames.txt -PriorityAnno priority_file.txt [--bowtie or --sam]

ann_filenames.txt:
./blat_results/tRf-1.coords
./blat_results/tRf-2.coords
./blat_results/tRf-3.coords
./blat_results/tRf-4.coords
./blat_results/piR-1.coords

sam_filenames.txt
./GCPool_sam_ctrl/GCLPool11-B6-2.sam
./GCPool_sam_ctrl/GCLPool11-B6-3.sam
./GCPool_sam_ctrl/GCLPool11-B6-4.sam
./GCPool_sam_ctrl/GCLPool12-B6-8.sam

priority_file.txt
mirna
tRF
pirna
rfam

bowtie input format:
2|1	+	chr18	21191746	TGGGTCAGGGTTTCGTACGTA	IIIIIIIIIIIIIIIIIIIII	2	18:A>G
2|1	-	chr12	127650511	TACGTACGAAACCCTGACCCA	IIIIIIIIIIIIIIIIIIIII	2	6:C>T
2|1	+	chrUn_gl000220	118338	TGGGTCAGGGTTTCGTACGTA	IIIIIIIIIIIIIIIIIIIII	2	6:G>A
13|13	+	chr1	229673604	GCATTGGTGGTTCAGTGGTAAAATTCTCA	IIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	20:G>A
13|13	-	chrX	64237838	TGAGAATTTTACCACTGAACCACCAATGC	IIIIIIIIIIIIIIIIIIIIIIIIIIIII	1	20:C>T
14|7	+	chr9	90606462	GGCCATACCACCCTGAACGCA	IIIIIIIIIIIIIIIIIIIII	20


Description:
1) Takes in coordinates of piRNA/tRFs and *.sam or bowtie alignment files.
2) If a read from the *.sam file overlaps one of the coordinates in the
piRNA/tRF coordinates, it will output the read and corresponding annotation will be writtten to an output file.
3) Each read is only counted once.
4) once a match occurs, then the matching ends and the next read is analyzed

Overlap Definition
Annotation start must be within -3 or +2 of read alignment start site.
Annotation end must be within +3 or -2 of read alignment end site.

Output:
output from bowtie:
#chr	strand	start	end	ID	mismatch_position	mismatch_type	name	chr	start	end	species
chr17	-	1617210	1617230	229|3	4	A>T	hsa-miR-22-3p	chr17	1617208	1617229	mirna
chr9	-	131154912	131154930	469|12	17	C>T	hsa-miR-219a-2-3p	chr9	131154911	131154932	mirna
chr8	-	22102487	22102510	546|15	8	C>A	hsa-miR-320a	chr8	22102488	22102509	mirna
chr16	+	69967006	69967028	719|5	21	G>A	hsa-miR-140-5p	chr16	69967006	69967027	mirna
chr3	+	183959193	183959213	771|4	4	G>A	hsa-miR-1224-5p	chr3	183959193	183959211	mirna
chr20	+	33578211	33578233	829|1	11	A>G	hsa-miR-499a-5p	chr20	33578211	33578231	mirna

output from sam file:
2970813|6	16	chr5	157403919	255	30M	*	0	0	GGAACAGATACTACACTTGATCTTAGCCAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:1	MD:Z:0A29	NM:i:1
2970820|11	16	chrX	37393165	255	18M	*	0	0	CATCCGATCAGATTAAAT	IIIIIIIIIIIIIIIII  XA:i:1	MD:Z:5T12	NM:i:1
2970820|11	16	chr4	132350522	255	18M	*	0	0	CATCCGATCAGATTAAAT	IIIIIIIIIIIIIIIII XA:i:1	MD:Z:5C12	NM:i:1
2970820|11	16	chr5	103932966	255	18M	*	0	0	CATCCGATCAGATTAAAT	IIIIIIIIIIIIIIIII XA:i:1	MD:Z:10C7	NM:i:1
2970821|1	0	chr6	72122472	255	20M	*	0	0	AGGAATGAGGACATGTTGCT	IIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:20	NM:i:0

summary file;
file:	./exRNA_plasma/Sample_1S3.v1.k100.beststrata.bwt
total number of reads:	21914717
breakdown of reads aligned to annotations
trf-1	6
trf-5	1961229
trf-3	3849
mirna	1979809
piRNA	18859


Notes:

annotation files have the format of
#chr	strand	start	end	species	name
chr1	-	100230891	100230921	piRNA	hsa_piR_016792
chrX	+	10040642	10040672	piRNA	hsa_piR_016792
chr9	+	100968492	100968522	piRNA	hsa_piR_016792


'''
import argparse
import os
import sys
import math
import re
import itertools
import operator
from subprocess import Popen, PIPE
import subprocess
from collections import defaultdict
import time
from datetime import datetime

MINLENGTH = 18



def __unicode__(self):
   return unicode(self.some_field) or u''

def get_arguments():

    parser = argparse.ArgumentParser(description='take in default bowtie alignments and output 3 files,\
    1) all reads with their assigned annotation\
    2) summary file of all annotation species and read counts\
    3) summary file of all annotation classes and read counts\
    ')


    parser.add_argument('-s ', dest='sam_filenames', type=str,
                    help='name of input file')
    parser.add_argument('-a ', dest='anno', type=str,
                    help='name of output file')
    parser.add_argument('--o ', dest='output_file', type=str,
                    help='path of directory that contains input file')
    parser.add_argument('-PriorityAnno ', dest='priority', type=str,
                        help='file that contains the priority of the annotation classes, each line is an annotation class name followed by a tab with priority number')
    parser.add_argument('--bowtie ', action='store_true',
                    help='alignment type is from a bowtie file')
    parser.add_argument('--sam ', action='store_true',
                    help='alignment output type is a sam file')
    parser.add_argument('--bam ', action='store_true',
                    help='alignment output type is a bam file')
    # parser.add_argument('-e ', dest='exp', type=str,
                    # help='a directory containing the experimental files from my predictions ')

    args = parser.parse_args()
    # if not args.sam_filenames or not args.anno or not args.output_file:
    #     print parser.error("must input a file with all samfiles")

    if len(sys.argv) < 6:
        parser.print_help()
        exit()
    fin = args.sam_filenames
    annotations = args.anno
    output_dir = args.output_file
    if output_dir == None:
        output_dir = ''
    priority = args.priority
    d= vars(args)
    alignment_type = ''
    if d['bowtie ']==True:
	    alignment_type = "bowtie"
    if d['sam ']==True:
	    alignment_type= "sam"
    if d['bam ']==True:
	    alignment_type= "bam"
    if alignment_type != "bowtie" and alignment_type != "sam" and alignment_type != "bam":
    	print 'you need to specify --bowtie, --bam, or --sam as the alignment output\n'.format(alignment_type)
        parser.print_help()
        exit()

    return fin, annotations, alignment_type, output_dir, priority

def run_command(cmd, stdin=0):
    if not isinstance(cmd, list):
        cmd = cmd.split(' ')
    if stdin==0:
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        stdout,stderr= p.communicate()
    else:
        p= Popen(cmd, stdout=PIPE, stderr=PIPE, stdin= PIPE)
        stdout, stderr= p.communicate(stdin)
    if p.returncode != 0:
        sys.stderr.write('command failed')
        sys.stderr.write(stdout)
        sys.stderr.write(stderr)
        sys.exit(1)

    return stdout,stderr

def annotation_dict(annot_filenames):

    num_annotations = 0
    annotation_dictionary = {}
    #key= (chromosome,strand,binstart) tuple
    #value = (start:end:ID:species)

    with open(annot_filenames, 'r') as fin:
        for anno_file in fin:
            anno_file = anno_file.strip('\n')
            if anno_file.startswith('#'):
                continue
            with open(anno_file, 'r') as fin2:
                print 'opening annotation file: {}'.format(anno_file)
                fin2.readline()
                for line in fin2:

                    # if line.startswith('#') or line.startswith('chr'):
                        # continue
                    line = line.strip().split()

                    name = line[5] #hsa-miR-424-3p
                    species = line[4] #mirna or pirna

                    chrom, strand, start, end = line[0], line[1], int(line[2]), int(line[3])

                    # for the dictionary, we want the key (look up term) to be the chromosome start and end and strand.
                    # the key should be the name of the annotation.
                    ######## adding dictionaries #########

                    binstart = math.floor((start)/1000)
                    binend = math.floor((end)/1000)

                    while binend >= binstart:
                        #key= (chromosome,strand,binstart) tuple
                        key = (chrom, strand, binstart) #(chr1,-,1120000)
                        #value = (start:end:ID:species)
                        value = ':'.join([str(start), str(end), name, species]) #'1103293:1103320:miR-200a-3p:mirna'
                        if key not in annotation_dictionary:
                            annotation_dictionary[key] = [value]
                        else:
                            annotation_dictionary[key].append(value)
                        binstart +=1

    print 'successfully created annotation dictionary'

    return annotation_dictionary



########################################################################
################# Sam alignment file type ###########################
########################################################################

def ID_count(ID):
    if len(ID.split('|')) == 1:
        return 1
    else:
        count = ID.split('|')[1]
        return float(count)

def check_stderr(stderr):
    if stderr != '':
        sys.stderr.write(stderr)
        sys.exit(1)
    return

def mmap_counter(alignment_file):

    c1 = ['sed', '/^@/d', alignment_file]

    p1 = subprocess.Popen(c1,stdout=subprocess.PIPE)
    p2 = subprocess.Popen('cut -f 1'.split(),stdin=p1.stdout, stdout=subprocess.PIPE)
    # p2_1 = subprocess.Popen(['sed', '-i ', '/^\s*$/d'],stdin=p2.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen('uniq -c'.split(),stdin=p2.stdout, stdout=subprocess.PIPE)
    p4 = subprocess.Popen(['sed', 's/^ *//g'],stdin=p3.stdout, stdout=subprocess.PIPE)
    output = p4.communicate()[0]
    mmap_count_dict = {}
    output = output.split()
    total = len(output)
    for i in range(0,total,2 ):
        count = output[i]
        ID = output[i+1]
        ID = ID.strip()
        mmap_count_dict[ID] = float(count)
    # print mmap_count_dict
    return mmap_count_dict


def make_priority_dict(priority_file):
    #priority file will look like the following
    # 1   miRNA
    # 2   trf
    # 3   pirna
    #key = RNA species
    #value = rank number
    rank_dict = {} #{'miRNA':1, 'trf':2, 'piRNA':3}
    with open(priority_file, 'r') as f_in:
        count = 0
        for line in f_in:
            if line.startswith('#'):
                continue
            else:
                line= line.strip().split('\t')
                if len(line) == 2:

                    try:
                        rank_num = int(line[0])
                        rank_dict[line[1]] = rank_num
                    except ValueError:
                        print 'priority file is incorrectly formated. First column must be rank number and second column is the species annotation name'
                        sys.exit()

                else: #the case where the numbers are missing but the line number represents the ranking
                    rank_dict[line[0]] = count
                    count +=1
    print rank_dict
    return rank_dict

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

def update_nte_field(chrom, strand, reference, position, read_nte):
    """
    update nte_dictionary with key ('position:ref>nte') and value (read_count) so we can give it back to the annotation object

    """
    position = str(position)
    modification_profile = ''.join([reference, '>', read_nte])
    key = ':'.join([chrom, strand, position, modification_profile])  # '22102508:TT>AA'
    # print key

    return key


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


def add_modification_position_fields(anno_chrom, anno_strand, read_start, read_end, read_nte5, read_nte3, cigar, read_seq, read_count):
    '''
    Add 3 more NTE fields, that show the chrom:strand:position:ref>minor, for each modification in read

    '''

    cigar = cigar.split(':')[2]  # MD:Z:0A0C20 --> '0A0C20'
    x = re.split("[ACTG]", cigar)  # ex: ['0', '20']  #2 5' ex: ['0', '0', '20']
    ref = re.sub("[^a-zA-Z]", "", cigar)  # 'AC'
    original_ref = ref

    reference, position = '', 0

    nte3_profile = '0'
    nte5_profile = '0'
    internal_profile = '0'
    internal = '0'
    if read_nte3 != "0":
        # if the read has a 3'NTE, then find out the 1) positions of the mismatches 2) NTE nucleotides 3) the reference nucleotides and
        reference, position = get_nte3_position(ref, read_nte3, read_start, read_end, anno_strand)
        # updating nte field profile
        nte3_profile = update_nte_field(anno_chrom, anno_strand, reference, position, read_nte3)

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

        # updating nte5 profile
        nte5_profile = update_nte_field(anno_chrom, anno_strand, reference, position, read_nte5)

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

        index_counter = 0
        internal_modification_list = []
        # print x
        minor_allele_list = []
        for m in range(len(x_temp) - 1):  # positive strand = [0,0,4,12]
            minor_allele = ''
            # print read_seq, int(x_temp[m]), shift_start, index_counter
            position = int(read_start) + int(x_temp[m]) + shift_start
            minor_allele = read_seq[int(x_temp[m]) + shift_start + index_counter]  # find the minor allele
            reference = ref[m]
            if anno_strand == "-":  # strand is negative need to get compl.
                reference = reverse_complement(reference)
                minor_allele = reverse_complement(minor_allele)
            minor_allele_list.append(minor_allele)
            # updating dictionary
            internal_modification = update_nte_field(anno_chrom, anno_strand, reference, position, minor_allele)
            internal_modification_list.append(internal_modification)
            index_counter += 1
        internal = '|'.join(minor_allele_list)
        internal_profile = '|'.join(internal_modification_list)
    new_fields = [nte5_profile, nte3_profile, internal_profile, internal]

    return new_fields


def sam_alignment_parser(sam_alignments, annotation_dictionary, priority_file):
    unannotated = 0
    # priority =  {'trf-1':1, 'trf-3':2, 'trf-5':3, 'mirna': 4, 'piRNA': 5 }
    priority_anno_dict = make_priority_dict(priority_file)
    total_read_count = 0
    summary_output = open(sam_alignments+'.summary', 'w')
    annotation_counts_file = open(sam_alignments+'.anno_counts', 'w')
    matched_annotation = [] #this matched annotation will contain a list of 'ID|abundance'
    cmd = 'wc -l {}'.format(sam_alignments)
    num_alignments, stderr = run_command(cmd)
    num_alignments= int(num_alignments.lstrip().split(' ')[0])
    # print num_alignments
    # sys.exit()
    percent = [ i/100.0 for i in range(0,100,10)]
    p2 = [int(i*num_alignments) for i in percent]
    print p2
    line_count = 0

    with open(sam_alignments, 'r') as alignments:
        print 'opening alignment file: {}'.format(sam_alignments)
        prev_alignment_name = 'Empty' #init
        prev_matched= False #init as uniquely mapped
        all_annotations = {} # dictionary that keeps count of all of the annotations (this will be used to make scatter plot)
        mmap_count_dict = mmap_counter(sam_alignments) #dictionary that keeps track of the number of times a read is multimapped;
        #mmap_count_dict = {'ID3|count': 1, 'ID1|count': 10, 'ID2|count': 32} #value is # of mapping locations
        #in the end, to assign read count, we can distribute the count by 1/n, n=number of mapping locations

        # initialize all_annotation dictionary
        for key in annotation_dictionary: #each value is a list of multiple annotations belonging to that bin/strand
            for i in annotation_dictionary[key]:
                i = i.split(':') #start:end:name:species #'1103293:1103320:miR-200a-3p:mirna'
                annotationID, species =  i[2], i[3]
                all_annotations[(annotationID,species)] = 0

        species_dict = {} #keeps count of the number of reads aligning to sa specific species (i.e species in priority dictionary)

##CHANGE
        for key in priority_anno_dict:
            species_dict[key] = 0
##CHANGE

        alignments.readline()

        for line in alignments:
            line_count += 1

            if line_count in p2:
                print "{}% Completed".format(line_count*100/num_alignments)
            if line.startswith('@'):
                continue

            line = line.strip().split()
            nte5 = line[7]
            nte3 = line[8]
            if line[1] == '4':
                continue
            try:
                mismatch = line[13]
            except IndexError:
                # print line
                mismatch = [s for s in line if "NM:i:" in s]

            try:
                md_flag= line[12]
            except IndexError:
                # print line
                md_flag = [s for s in line if "MD:Z:" in s]


            ID, chrom, start, sequence  = line[0], line[2], int(line[3]), line[9]

            if ID != prev_alignment_name:
                total_read_count += ID_count(ID)
                prev_alignment_name = ID


            if line[1] == '16':
                strand = '-'
            if line[1] =='0':
                strand= '+'

            end = start + len(line[9])-1
            binstart = math.floor(start/1000) #used as a key to look up bin in dictionary
            binend = math.floor(end/1000)
            key = (chrom, strand, binstart)
            annotation_options = ''
            counted = False #This is a flag to make sure that if the alignment goes over 2 bins, and it was already
                                #annotated, then it will not double count
            species_rank1 = ''
            species_rank2 = ''
            while binstart <= binend:
                # if counted == True:
                #     break
                if key in annotation_dictionary:
                    annotations = annotation_dictionary[key] #gets all the possible annotations for that bin
                    # print annotations
                    for i in annotations:
                        original = i #annotation info 'start:end:name:species'
                        i = i.split(':')

                        anno_start, anno_end, anno_name, anno_species = int(i[0]), int(i[1]), i[2], i[3]

                    ##### version 8 edit start ##############################
                        if abs(anno_start - anno_end)+1 > 50:

                            if (start>=anno_start-3) and (end<=anno_end+3):
                                # print 'long annotation {}'.format(i)

                                # print 'HIT'
                                ########### HIT!!! READ WITH SAME COORDINATES AS ANNOTATION #################
                                if annotation_options == '':
                                    annotation_options = original #holds the result of the annotation hit
                                    if anno_species in priority_anno_dict:
                                        species_rank1 = priority_anno_dict[anno_species]
                                    else:
                                         priority_anno_dict[anno_species] = len(priority_anno_dict) +1
                                         species_rank1 = len(priority_anno_dict) +1
                                else:#another annotation already exists for this read

                                    if anno_species in priority_anno_dict:
                                        species_rank2 = priority_anno_dict[anno_species]
                                    else:
                                        species_rank2 = len(priority_anno_dict) + 1
                                        priority_anno_dict[anno_species] = len(priority_anno_dict) +1
                                    # compare priority, if species rank2 has a lower number (aka higher rank, then change the annotation)
                                    old = annotation_options
                                    if species_rank2 < species_rank1:
                                        annotation_options=original
                                    # print 'conflict: {}:{}: new = {}, old = {}\twinner: {}'.format(chrom, strand, old, original, annotation_options)
                        else: #annotation is 50 nucletoides or less ##### version 8 edit
                            if (start>=anno_start-3 and start<=anno_start+3) and (end<=anno_end+3 and end>=anno_end-3): ##### version 8 edit
                            # print 'HIT'
                            ########### HIT!!! READ WITH SAME COORDINATES AS ANNOTATION #################
                                if annotation_options == '':
                                    annotation_options = original #holds the result of the annotation hit
                                    if anno_species in priority_anno_dict:
                                        species_rank1 = priority_anno_dict[anno_species]
                                    else:
                                         priority_anno_dict[anno_species] = len(priority_anno_dict) +1
                                         species_rank1 = len(priority_anno_dict) +1
                                else:#another annotation already exists for this read

                                    if anno_species in priority_anno_dict:
                                        species_rank2 = priority_anno_dict[anno_species]
                                    else:
                                        species_rank2 = len(priority_anno_dict) + 1
                                        priority_anno_dict[anno_species] = len(priority_anno_dict) +1
                                    # compare priority, if species rank2 has a lower number (aka higher rank, then change the annotation)
                                    old = annotation_options
                                    if species_rank2 < species_rank1:

                                        annotation_options=original
                                    # print 'conflict: {}:{}: new = {}, old = {}\twinner: {}'.format(chrom, strand, old, original, annotation_options)
                binstart += 1

            # {0: 'mirna', 1: 'tRf', 2: 'piRNA', 'piRNA': 5, 'tRf': 6, 'mirna': 4}?

            # i = annotation_options
            if annotation_options != '':
                # print annotation_options
                annotation_options= annotation_options.split(':')

                anno_start, anno_end, anno_name, anno_species = int(annotation_options[0]), int(annotation_options[1]), annotation_options[2], annotation_options[3]


                all_annotations[(anno_name, anno_species)] += (ID_count(ID)/mmap_count_dict[ID])

                if anno_species in species_dict:

                    species_dict[anno_species] += (ID_count(ID)/mmap_count_dict[ID])
                else:
                    species_dict[anno_species] = (ID_count(ID)/mmap_count_dict[ID])

                output_string = [chrom, strand, str(start), str(end), ID, nte5, nte3, mismatch, md_flag, sequence, anno_name, chrom, str(anno_start), str(anno_end), anno_species]

                read_count = ID.split('|')[1]
                added_NTE_fields = add_modification_position_fields(chrom, strand, str(start), str(end), nte5, nte3, md_flag, sequence, int(read_count))
                output_string.extend(added_NTE_fields)

                string = '\t'.join(output_string)
                matched_annotation.append(string)

############# Summary file ##############
# only output total read counts in species
    print priority_anno_dict
    summary_output.write('file:\t{}\n'.format(sam_alignments))
    summary_output.write('total number of reads:\t{}\n'.format(total_read_count))
    summary_output.write('breakdown of reads aligned to annotations\n')
    for key in species_dict:
        summary_output.write('{}\t{}\n'.format(key, species_dict[key]))

#header

    annotation_counts_file.write("annotationID\tspecies\tcount\tRPM\n")

    for i in sorted(all_annotations, key=all_annotations.get, reverse=True):
        if all_annotations[i]==0:
            rpm=0
        else:
            rpm = float(all_annotations[i])*1000000/total_read_count
            rpm = "{0:.2f}".format(round(rpm,2))

        output_line = '\t'.join([i[0], i[1], str(all_annotations[i]), str(rpm)])
        annotation_counts_file.write(output_line)
        annotation_counts_file.write('\n')
    return matched_annotation


#################### FIX SCRIPT TO REMOVE PRIORITY AND JUST ALLOCATE READS TO ALL ALIGNMENTS ######################

def bowtie_alignment(bowtie_alignments, annotation_dictionary):
    unannotated = 0
    print " bowtie format"
    summary_output = open(bowtie_alignments+'.summary', 'w')
    annotation_counts_file = open(bowtie_alignments+'.anno_counts', 'w')
    matched_annotation = [] #this matched annotation will contain a list of 'ID|abundance'

    with open(bowtie_alignments, 'r') as alignments:
        total_read_count = 0
        print 'opening alignment file: {}'.format(bowtie_alignments)
        prev_alignment_name = 'Empty' #init
        prev_matched= False
        all_annotations = {} # dictionary that keeps count of all of the individual annotations names and the counts for each species (this will be used to make scatter plot)
        mmap_count_dict = mmap_counter(bowtie_alignments)

        readID_dict = defaultdict(list) # dictionary to hold the annotations for each read ID, this is to ensure double counting of multimapped reads does not occur
        # initialize all_annotation dictionary
        for key in annotation_dictionary: #each value is a list of multiple annotations belonging to that bin/strand
            for i in annotation_dictionary[key]:
                i = i.split(':') #start:end:ID:species
                annotationID, species =  i[2], i[3]
                all_annotations[(annotationID,species)] = 0
        species_dict = {} #keeps count of the number of reads aligning to the specific class (i.e class in priority dictionary)
##CHANGE
        for key in priority_anno_dict:
            species_dict[key] = 0
##CHANGE
        for line in alignments:
            if line.startswith('@'):
                continue
            line = line.strip().split()
            if len(line) == 8: #if there are 8 fields then there is a mismatch, and the mismatch profile is added at the end of hte alignemnt file. i.e 10:A>G  (position:ref>alt)
                mismatch_line = line[7].split(',')
                mismatch_position= []
                mismatch_type= []
                for i in mismatch_line:
                    mp = i.split(':')[0]
                    mismatch_position.append(mp)
                    mt = i.split(':')[1]
                    mismatch_type.append(mt)
            else:
                mismatch_position = ['NA']
                mismatch_type = ['NA']
            #adding 1 because bowtie alignment starts are 0 based
            ID, strand, chrom, start, sequence = line[0], line[1], line[2], int(line[3])+1, line[4]

            if ID != prev_alignment_name:
                total_read_count += ID_count(ID)
                prev_alignment_name = ID


            end = start + len(line[4])-1
            binstart = math.floor(start/1000) #used as a key to look up bin in dictionary
            binend = math.floor(end/1000)
            key = (chrom, strand, binstart) #key for variable annotation_dictionary
            annotation_options = [] #dictionary of multiple options of annotations in the bin


            counted = False #This is a flag to make sure that if the alignment goes over 2 bins, and it was already
                                #annotated, then it will not double count
            species_rank = ''
            while binstart <= binend:
                # if counted == True:
                #     break
                if key in annotation_dictionary: #search the read's coordinate bin in the annotation dictionary (chrom,strand,binstart)
                    annotations = annotation_dictionary[key] #look up all the possible annotations for that bin


                    for i in annotations: #loop through all annotations within the bin key=(chrom, strand, binstart)
                        original = i #each annotation info 'ann_start:ann_end:ann_name:ann_species'
                        i = i.split(':')

                        anno_start, anno_end, anno_name, anno_species = int(i[0]), int(i[1]), i[2], i[3]

                        ########### HIT!!! READ WITH SAME COORDINATES AS ANNOTATION #################
                        if (start>=anno_start-3) and (end<=anno_end+3): #if the alignment coordinate is within the annotation coordinate
                            if annotation_options == []:
                                annotation_options.append(original) #holds the result of the annotation hit
                                species_rank1 = priority_dict[anno_species]
                            else:
                                species_rank2 = priority_dict[anno_species]
                                if species_rank2 < species_rank1:
                                    annotation_options=[original]
                            counted ==True
                binstart += 1


            if not annotation_options: #checking to see if dictionary is empty, if so there are no annotations within this alingnment bin,
                continue    #since no annotations for this alignment, can move to next alignment
            else:

                # chosen_annotation  = min(annotation_options.iteritems(), key=operator.itemgetter(1))[0] #picks the annotation with the highest rank (i.e pick priority 1 over priority 4)

                for i in annotation_options:
                    i = i.split(':')
                    anno_start, anno_end, anno_name, anno_species = int(i[0]), int(i[1]), i[2], i[3]
                    all_annotations[(anno_name, anno_species)] += (ID_count(ID)/mmap_count_dict[ID])

                    if anno_species in species_dict:
                        # print ID_count(ID)/mmap_count_dict[ID]
                        species_dict[anno_species] += (ID_count(ID)/mmap_count_dict[ID])
                    else:
                        species_dict[anno_species] = (ID_count(ID)/mmap_count_dict[ID])

                    output_string = [chrom, strand, str(start), str(end), ID, ','.join(mismatch_position), ','.join(mismatch_type),  sequence, anno_name, chrom, str(anno_start), str(anno_end), anno_species]
                    string = '\t'.join(output_string)
                    matched_annotation.append(string)



############# Summary file ##############
# only output total read counts in species
    summary_output.write('file:\t{}\n'.format(bowtie_alignments))
    summary_output.write('total number of reads:\t{}\n'.format(total_read_count))
    summary_output.write('Reads mapped to annotations\n')
    for key in species_dict:
        summary_output.write('{}\t{}\n'.format(key, species_dict[key]))

#header

    annotation_counts_file.write("annotationID\tspecies\tcount\tRPM\n")

    for i in sorted(all_annotations, key=all_annotations.get, reverse=True):
        if all_annotations[i]==0:
            rpm=0
        else:
            rpm = float(all_annotations[i])*1000000/total_read_count
            rpm = "{0:.2f}".format(round(rpm,2))

        output_line = '\t'.join([i[0], i[1], str(all_annotations[i]), str(rpm)])
        annotation_counts_file.write(output_line)
        annotation_counts_file.write('\n')

    return matched_annotation

# def __unicode__(self):
   # return unicode(self.some_field) or u''

def main():
    start_time = datetime.now()
    sam_file, annotation_files, alignment_type, output_dir, priority = get_arguments()
    if output_dir != '':
        if not output_dir.endswith('/'):
            output_dir= output_dir+'/'
        if output_dir == './':
            output_dir=''

    with open(annotation_files, 'r') as anno:
        anno_dict = annotation_dict(annotation_files)

    if alignment_type =="bam":
        # with open(sam_file, 'r') as fin:
        #     for bam_file in fin:
        bam_file = sam_file.strip('\n')
        output_fn = bam_file+'.anno_overlap'
        output_file = open(output_fn, 'w')
        output_file.write('#chr\tstrand\tstart\tend\tID\t5NTE\t3NTE\tmismatch_position\tmismatch_type\tsequence\tann_name\tann_chr\tann_start\tann_end\tspecies\tnte5_profile\tnte3_profile\tinternal_profile\tINTERNAL\n') #header
        cmd = 'samtools view -o temp.sam {}'.format(bam_file)
        cmd = cmd.split(' ')
        ps = Popen(cmd,stdout=PIPE,stderr=PIPE)
        stdout, stderr = ps.communicate()
        check_stderr(stderr)
        sam_file = bam_file.split('.')[0]
        sam_file = sam_file+'.sam'
        cmd = 'mv temp.sam {}'.format(sam_file)
        ps = Popen(cmd,shell=True, stdout=PIPE,stderr=PIPE)
        stdout, stderr = ps.communicate()
        check_stderr(stderr)
        matches = sam_alignment_parser(sam_file, anno_dict)
        print 'writing output for {}'.format(bam_file)
        cmd = 'rm {}'.format(sam_file)
        for i in matches:
			output_file.write(i)
			output_file.write('\n')
# SAM alignment
    if alignment_type =="sam":

        # if sam_file.startswith('#'):
        #     continue
        sam_file = sam_file.strip('\n')
        # sam_name = sam_file.split('/')[-1]
        # output_fn = output_dir + sam_name + '.anno_overlap'
        output_fn = sam_file + '.anno_overlap'
        output_file = open(output_fn, 'w')
        output_file.write('#chr\tstrand\tstart\tend\tID\t5NTE\t3NTE\tmismatch_position\tmismatch_type\tsequence\tann_name\tann_chr\tann_start\tann_end\tspecies\tnte5_profile\tnte3_profile\tinternal_profile\tINTERNAL\n') #header
        matches = sam_alignment_parser(sam_file, anno_dict, priority)
        print 'writing output for {}'.format(output_fn)
        for i in matches:
        	output_file.write(i)
    		output_file.write('\n')
# bowtie alignment format
    if alignment_type =="bowtie":
    #if combine all alignment files into one large file.

        # with open(sam_filenames, 'r') as fin:
        #     for sam_file in fin:
        #         sam_file = sam_file.strip('\n')
        #         matches = bowtie_alignment(sam_file, anno_dict)
        #         print 'writing output for {}'.format(sam_file)
        #         for i in matches:
        #     		output_file.write(i)
        #   	    	output_file.write('\n')
        # with open(sam_filenames, 'r') as fin:
        #     for sam_file in fin:
        sam_file = sam_file.strip('\n')
        # sam_name = sam_file.split('/')[-1]
        output_fn = sam_name+'.anno_overlap'
        output_file = open(output_fn, 'w')
        output_file.write('#chr\tstrand\tstart\tend\tID\t5NTE\t3NTE\tmismatch_position\tmismatch_type\tsequence\tann_name\tann_chr\tann_start\tann_end\tspecies\tnte5_profile\tnte3_profile\tinternal_profile\tINTERNAL\n') #header
        matches = bowtie_alignment(sam_file, anno_dict)
        print 'writing output for {}'.format(sam_filenames)
        for i in matches:
            output_file.write(i)
            output_file.write('\n')
    end_time=datetime.now()
    print 'Duration: {}'.format(end_time-start_time)
    return

def test():
    sam_file = '~/mount/mismatch_modification_project/scripts/trimming_pipe_test/Sample_3S20.NTE.sam'
    annotation_files = 'test.paths'
    alignment_type= 'sam'
    anno_dict = annotation_dict(annotation_files)
    # count =0

    # for key in annotation_dict:
    #     print '{}:{}'.format(key, annotation_dict[key])
    #     count +=1
    #     if count ==10:
    #         sys.exit()

    # matches = sam_alignment_parser(sam_file, anno_dict)
if __name__ == "__main__":
	# test()
    main()
