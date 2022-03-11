#!/bin/python

'''
Author: Kikuye Koyano
Contant: kkoyano@ucla.edu
Date: 06/04/2018

* * * * UPDATES * * * *
06/22/2018:
    make the output file to reflect the sample and recount threshold
    <outputfn>.credible-r${read_count}-N${sample_count}

06/26/2018:
in function filter_modfications_under_threshold(), if the total number of samples is less than the sample number threshold, then change the sample number threshold to be the total number of samples.



Description:
#Part 1: Get all of the credible modification sites
    #   - go through each of the <sample_name>.anno_summary
    #   - for each modification type, add if the modification in the sample had at least 2 reads, add the modification key to dictionary 'chr:strand:position:ref>alt' and add +1 to sample count (value)
    #   - do not record modification if the position only has < read_coverage_threshold
    # end result:
    # 3 dictionaries with all modification profiles, and the number of samples they are observed in

# Part 2:
    # go back through each 3 dictionaries, and remove key if the profile only occurs in 1 sample
    # this will leave us with a dictionary containing a credible set of modifications

# Part 3:
    # go back through annotation summary file. if the modification in the read is not considered a credible modification, then throw out the modification for that read

To be run on the *.anno_summary file
Filtering for credible modifications
    - modifications that have at least -r <int> reads supporting modification
    - present in at least -t <int> number of samples within group


Input:
     -r <int>
     -t <int>
     -i <file> # list of sample_names in group

Output:
    <sample_name1>.anno_summary.credible
    ...
    ...
    <sample_nameN>.anno_summary.credible

Notes:
    Usage:  python ~/mount/scripts/mismatch_profile/filter_credible_modifications.py -r 2 -t 2 -o fout.txt -i input_files.txt
    cat input_files.txt # tab delimited file with <sample>.anno_summary file in col1, col2 = <sample>.anno_overlap.nte_fields
        /Users/kiku/mount/scripts/mismatch_profile/test/credible_filter/test1.anno_summary  /Users/kiku/mount/scripts/mismatch_profile/test/credible_filter/test1.anno_overlap.nte_fields
        /Users/kiku/mount/scripts/mismatch_profile/test/credible_filter/test2.anno_summary  /Users/kiku/mount/scripts/mismatch_profile/test/credible_filter/test2.anno_overlap.nte_fields
'''
import argparse
import os
import subprocess
import re
import gzip
from datetime import datetime

# f = gzip.open(filename, 'rb') if re.search(r'\.gz$', filename) else open(filename, 'r')


def get_arguments():

    parser = argparse.ArgumentParser(description='Get credible set of modifications from the sample')

    parser.add_argument('-i ', dest='input_file', type=str,
                        help='name of samples in input file group, 2 tab-delimited columns: col1 = anno_summary file, col2 = anno_overlap file')
    parser.add_argument('-r ', dest='reads_threshold', type=int,
                        help='minimum number of reads that support modification')
    parser.add_argument('-t ', dest='threshold_sample', type=int,
                        help='minimum number of samples the modification must occur in')
    parser.add_argument('-o ', dest='output_group_filename', type=str,
                        help='name of output file to summarize modification sample counts')
    parser.add_argument('--output_anno_summary', action='store_true',
                        help='flag that indicates the output should be the anno-summary rather than the NTE.sam format')
    # parser.add_argument('-e ', dest='exp', type=str,
    # help='a directory containing the experimental files from my predictions ')

    args = parser.parse_args()
    # if len(sys.argv) < 2:
    # parser.print_help()
    # exit()
    fin = args.input_file
    fout = args.output_group_filename
    reads_threshold = args.reads_threshold
    threshold_sample = args.threshold_sample
    output_anno_summary = args.output_anno_summary
    # bed_file = args.bed

    return fin, fout, reads_threshold, threshold_sample, output_anno_summary


def check_file_exist(file_path):
    if os.path.exists(file_path) is False:
        print 'ERROR!!!'
        print '{} does not exist'.format(file_path)
    return

def update_modification_dict(chrom, strand, modification_field, modification_dict, reads_threshold):

    if modification_field == '0':
        return modification_dict
    else:
        #  unique modifications are separated by a pipe
        modification_field = modification_field.strip().split('|')  # ['58218440:G>T:12', '58218439:T>A:1', '58218440:G>A:6']
        # print modification_field
        for i in modification_field:

            i = i.split(':')  # ['58218440', 'G>T', '12']

            # if the reads for the modification are less than the designated reads_threshold, then continue. Not enough evidence that modification exists for this sample and cannot add to the credible set (will filter this out later).
            if int(i[2]) < reads_threshold:
                continue
            else:
                key = ':'.join([chrom, strand, i[0], i[1]])  # 'chr:strand:position:Ref>minor'
                if key in modification_dict:
                    modification_dict[key] += 1
                else:
                    modification_dict[key] = 1

    return modification_dict


def check_modification_in_set(chrom, strand, total_nte_count, modification_field, modification_set, reads_threshold, total_read_count):
    '''
    chrom, strand, total_nte5_count, nte5_profile, nte5_credible_modification_set, reads_threshold

    checks to see if all the modificaiton in the modification field are in the credible modification set, if they are not, modify the anno_modification_count (ex: total number of 3'NTEs for that given ncRNA annotation)

    return the updated anno_modification_type_count and modification field
    '''
    new_nte_count = 0

    credible_modifications = []
    if modification_field == '0':

        return total_nte_count, modification_field, total_read_count  # return the modification profile as '0', and keep the original number of anno_modification_type_count for annotation
    else:
        # unique modifications are separated by a pipe
        modification_field = modification_field.strip().split('|')  # ['58218440:G>T:12', '58218439:T>A:1', '58218440:G>A:6']

        for i in modification_field:
            original = i
            i = i.split(':')  # ['58218440', 'G>T', '12']
            modification_count = int(i[2])
            # if the reads for the modification are less than the designated reads_threshold, then don't count. Not enough evidence that modification exists and cannot add to the credible set (will filter this out later).
            if modification_count < reads_threshold:
                # total_read_count -= modification_count
                continue
            else:
                key = ':'.join([chrom, strand, i[0], i[1]])  # 'chr:strand:position:Ref>minor'
                if key in modification_set:  # if the modification is deemed credible, add to credible list
                    new_nte_count += modification_count
                    credible_modifications.append(original)
                # else:
                    # total_read_count -= modification_count
        # join all credible modifications
        credible_modifications = '|'.join(credible_modifications)
        if credible_modifications == '':
            credible_modifications = '0'
        return str(new_nte_count), credible_modifications, total_read_count

def filter_modfications_under_threshold(fout_obj, group_label, total_num_samples, modification_dict, threshold_sample):

    '''
    INPUT FILES: takes in modification dictionary, with key being the unique modification location, and the value being the number of samples that were found to contain that modification
    RETURNS: set of modifications that were in at least N <int> number of samples (-t 'threshold_sample')
    '''
    credible_modificatio_set = set()
    output_lines = ''

    if total_num_samples < threshold_sample:
        print 'total number of samples = {} is less than threshold samples ({}), adjust the threshold to be the total number of samples in this fluid'.format(total_num_samples, threshold_sample)
        threshold_sample = total_num_samples

    for key in modification_dict:
        output_line = '{}\t{}\t{}\t{}'.format(key, str(modification_dict[key]), str(total_num_samples), group_label)
        output_lines = '\n'.join([output_lines, output_line])
        if modification_dict[key] >= threshold_sample:
            credible_modificatio_set.add(key)
    fout_obj.write(output_lines)

    return credible_modificatio_set


def check_perfect_match(mismatch_flag):
    '''
    returns True is the alignemnt has NO mismatches
    returns False is the alignemnt has mismatches
    '''
    mismatch_flag = mismatch_flag.split(':')[2]
    # no mismatches in the read's alignment
    if mismatch_flag == '0':
        return True
    else:
        return False


def filter_anno_summary_for_credible_modifications(fin, reads_threshold, nte5_credible_modification_set, nte3_credible_modification_set, internal_credible_modification_set, threshold_sample):

    with open(fin, 'r') as fin_obj:
        for sample_fn in fin_obj:
            sample_fn = sample_fn.strip().split()[0]
            # # opening  <sample_name>.anno_summary  file
            # print sample_fn
            header = ''
            # annotation_name   species chromosome  strand  start   end anno_sequence   total_read_count    perfect_match_count anno_mismatch_count RPM total_nte5_count    total_nte3_count    total_internal_count    nte5    nte3    internal_mismatch
            print 'parsing modifications in {}'.format(sample_fn)
            output_lines_all = ''

            output_fn = sample_fn + '.credible-r{}-N{}'.format(reads_threshold, threshold_sample)
            with open(sample_fn, 'r') as sample_fn_obj:
                fout_obj = open(output_fn, 'w')
                for i in sample_fn_obj:
                    # skip if a header line
                    if i.startswith('#') or i.startswith('anno'):
                        output_lines_all = i.strip()
                        continue
                    i = i.strip().split()
                    chrom, strand, total_nte5_count, total_nte3_count, total_internal_count, nte5_profile, nte3_profile, internal_profile, total_read_count = i[2], i[3], i[11], i[12], i[13], i[14], i[15], i[16], int(i[7])
                    total_nte5_count, nte5_profile, total_read_count = check_modification_in_set(chrom, strand, total_nte5_count, nte5_profile, nte5_credible_modification_set, reads_threshold, total_read_count)
                    total_nte3_count, nte3_profile, total_read_count = check_modification_in_set(chrom, strand, total_nte3_count, nte3_profile, nte3_credible_modification_set, reads_threshold, total_read_count)
                    total_internal_count, internal_profile, total_read_count = check_modification_in_set(chrom, strand, total_internal_count, internal_profile, internal_credible_modification_set, reads_threshold, total_read_count)
                    # THE TOTAL READ COUNT HERE MAY BE INACCURATE BECAUSE THE READ MAY CONTAIN BOTH AN INTERNAL MODIFICATION AND 3'NTE (OR 5'NTE), SO BY SUBTRACTING THE TOTAL_READ_COUNT BY THE MODIFICATION COUNT, YOU MAY BE SUBTRACTING THE TOTAL READ COUNT TWICE
                    i[7], i[11], i[12], i[13], i[14], i[15], i[16] = total_read_count, total_nte5_count, total_nte3_count, total_internal_count, nte5_profile, nte3_profile, internal_profile

                    output = [str(j) for j in i]
                    output = '\t'.join(output)
                    output_lines_all = '\n'.join([output_lines_all, output])
                fout_obj.write(output_lines_all)
                fout_obj.close()
    return


def find_credible_modifications_MAIN(fin, fout, reads_threshold, threshold_sample, output_anno_summary):
    # check to see if the file is gziped
    fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')

    # modification dictionary, key = 'chrom:strand:position:ref>minor', value = number of samples with this exact modifcation
    nte3_modification_sample_count_dict = {}
    nte5_modification_sample_count_dict = {}
    internal_modification_sample_count_dict = {}

    #Part 1: Get all of the credible modification sites
    #   - go through each of the <sample_name>.anno_summary file in the input_file list and record if the modification was found
    #   - do not record modification if the position only has < read_coverage_threshold
    total_num_samples = 0
    for sample_fn in fin_obj:
        sample_fn = sample_fn.strip().split()[0]

        # # opening  <sample_name>.anno_summary  file

        print 'parsing modifications in {}'.format(sample_fn)
        with open(sample_fn , 'r') as sample_fn_obj:
            for i in sample_fn_obj:
                # skip if a header line
                if i.startswith('#') or i.startswith('anno'):
                    continue
                i = i.strip().split()
                # print i

                chrom, strand, nte5, nte3, internal_mismatch = i[2], i[3], i[14], i[15], i[16]
                nte3_modification_sample_count_dict = update_modification_dict(chrom, strand, nte3, nte3_modification_sample_count_dict, reads_threshold)
                nte5_modification_sample_count_dict = update_modification_dict(chrom, strand, nte5, nte5_modification_sample_count_dict, reads_threshold)
                internal_modification_sample_count_dict = update_modification_dict(chrom, strand, internal_mismatch, internal_modification_sample_count_dict, reads_threshold)
            total_num_samples += 1
    fin_obj.close()

    # CHECKPOINTL
    # checking dictionaries # these are prefiltered the N number of minimum samples
    # print 'nte3_sample_count: {}'.format(nte3_modification_sample_count_dict)
    # print 'nte5_sample_count: {}'.format(nte5_modification_sample_count_dict)
    # print 'internal_sample_count: {}'.format(internal_modification_sample_count_dict)
    # sys.exit()

    #Part 2: Go through the dictionary and filter out any modifications that do not occur in more than N individuals (-t threshold sample), also print out the results to later plot, see how many modifications are being removed here.
    # this output file contains the number of samples that had each modification, will use for histogram
    output_fn2 = fout+'.sample_count-r{}-N{}.sample_count'.format(reads_threshold, threshold_sample)
    fout_obj = open(output_fn2, 'w')
    fout_obj.write('modification\tsample_count\ttotal_num_samples_tested\tgroup_label')

    print 'parsing nte3 credible set'
    nte3_credible_modification_set= filter_modfications_under_threshold(fout_obj, 'nte3', total_num_samples,nte3_modification_sample_count_dict, threshold_sample)
    print 'parsing nte5 credible set'
    nte5_credible_modification_set= filter_modfications_under_threshold(fout_obj, 'nte5', total_num_samples,nte5_modification_sample_count_dict, threshold_sample)
    print 'parsing internal credible set'
    internal_credible_modification_set= filter_modfications_under_threshold(fout_obj, 'internal', total_num_samples,internal_modification_sample_count_dict, threshold_sample)
    fout_obj.close()

    # checkpoint to make sure the set only contains positions that have at least N number of samples also with modification position
    # print 'nte3_credible_set: {}'.format(nte3_credible_modification_set)
    # print 'nte5_credible_set: {}'.format(nte5_credible_modification_set)
    # print 'internal_credible_set: {}'.format(internal_credible_modification_set)
    # sys.exit()


    if output_anno_summary is True:
        filter_anno_summary_for_credible_modifications(fin, reads_threshold, nte5_credible_modification_set, nte3_credible_modification_set, internal_credible_modification_set, threshold_sample)
    else:
    #Part 3: Go back through each sample's SAM file and only print out the reads that are in the credible set, also output if the modification is perfectly matched, then this new output file with all credible modifications will be the used to do *.anno_summary and --> mismatch_profile
        fin_obj = gzip.open(fin, 'rb') if re.search(r'\.gz$', fin) else open(fin, 'r')
        for sample_fn in fin_obj:

            sample_fn = sample_fn.strip().split()[1]
            print 'parsing reads in {}'.format(sample_fn)
        # opening  <sample_name>.NTE.sam file
            with open(sample_fn , 'r') as sample_fn_obj:
                output_lines = ''
                for i in sample_fn_obj:
                    if i.startswith('#'):
                        continue
                    i = i.strip().split()
                    readID_count = i[4].split('|')[1]

                    # if the read is a perfect match, then add to output lines
                    # also add '0' to INTERNAL field (this field outputs the minor allele w.r.t. RNA, and if None exist, then output '0')
                    if check_perfect_match(i[7]) is True:
                        i.append('0')
                        i = '\t'.join(i)
                        output_lines = '\n'.join([output_lines, i])
                        continue

                    # check to see if the modification profile of the read appear in the credible set, if they do we need to output the line,
                    # if one of the fields is not in the credible set, set it to '0'
                    else:

                        nte5_profile, nte3_profile, internal_profile = i[15], i[16], i[17]

                        # if modification is in credible set then
                        output_flag = False
                        if nte5_profile in nte5_credible_modification_set: output_flag = True
                        else: i[5], i[15] = '0', '0'
                        if nte3_profile in nte3_credible_modification_set: output_flag = True
                        else: i[6], i[16] = '0' , '0'
                        all_internal_list = internal_profile.split('|')

## Changes
                        if len(all_internal_list) > 1:
                            internal_profile_out = []
                            internal_mod_out = []
                            for internal_mod in all_internal_list:
                                if internal_profile in internal_credible_modification_set:
                                    output_flag = True
                                    ref_to_minor = internal_profile.split(':')[3]
                                    minor = ref_to_minor.split('>')[1]
                                    internal_mod_out.append(minor)
                                    internal_profile_out.append(internal_profile)


                                else:
                                    internal_mod_out.append('0')
                                    internal_profile_out.append('0')
                                    # i.append('0') #field 18
                            i[17] = '|'.join(internal_profile_out)
                            i[18] = '|'.join(internal_mod_out)
                        else:
                            if internal_profile in internal_credible_modification_set:
                                output_flag = True
                                ref_to_minor = internal_profile.split(':')[3]
                                # minor = ref_to_minor.split('>')[1]
                                # i.append(minor)
                            else:
                                i[17] = '0'
                                i[18] = '0'
                                # i.append('0') #field 18
                    if output_flag is True:
                        i = '\t'.join(i)
                        output_lines = '\n'.join([output_lines, i])
            header_line = '#chr\tstrand\tstart\tend\tID\t5NTE\t3NTE\tmismatch_position\tmismatch_type\tsequence\tann_name\tann_chr\tann_start\tann_end\tspecies\tnte5_profile\tnte3_profile\tinternal_profile\tINTERAL'
            # print header_line
            # print output_lines

            output_fn3 = sample_fn + '.credible-r{}-N{}'.format(reads_threshold, threshold_sample)
            with open(output_fn3, 'w') as output_obj:
                output_obj.write(header_line)
                output_obj.write(output_lines)


    return


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


def main():
    start_time = datetime.now()

    # Testing files and parameters
    # fin = '/Users/kiku/mount/scripts/mismatch_profile/test/credible_filter/input_files.txt'
    # fout = '/Users/kiku/mount/scripts/mismatch_profile/test/credible_filter/fout.txt'
    # reads_threshold = 2
    # threshold_sample = 2
    # output_anno_summary = True

    fin, fout, reads_threshold, threshold_sample, output_anno_summary = get_arguments()
    find_credible_modifications_MAIN(fin, fout, reads_threshold, threshold_sample, output_anno_summary)
    end_time = datetime.now()
    print 'runtime = {}'.format(end_time - start_time)
    print 'DONE!'
    return


if __name__ == "__main__":
    main()
