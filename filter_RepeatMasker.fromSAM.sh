#!/bin/bash
# Author: Kikuye Koyano
# Contant: kkoyano@ucla.edu
# Date: 06/20/2018

# Description: Script to take in <sample_name>.anno_summary file (rows = miRNA species) and overlap with a RepeatMasker.bed file. Script will output a new <sample_name>.anno_summary.rmsk file that has all annotations overlapping with a repeat removed.

#
# Notes:
#   - bedtools intersect default the number of bases needed to intersect is 1bp
# Usage:
# filter_RepeatMasker.sh directory sample_name.NTE.sam RepeatMasker.bed

# for testing
# RepeatMasker_bed=/u/home/k/kkoyano/nobackup-gxxiao3/annotations/hg19/RepeatMasker/rmsk.bed
# anno_summary_file=U87_sRNA_R3_S3_L002_R1_001.NTE.sam.anno_overlap.t0.5.canonical_end.nte_fields.mirna.ALL.psl.best_score.bed.results.lvl3.anno_summary.nosnp.credible


# directory containing the anno_summary file
dir=$1
anno_summary_file=$2
RepeatMasker_bed=$3
# directory with intermediate file of miRNA species that intersect with regions from RepeatMasker
intermediate_output_dir=intermediate/RepeatMasker/

# directories to put the intermediate files.
mkdir -p ${dir}/intermediate
mkdir -p ${dir}/intermediate/RepeatMasker/

# Checking if the directory, anno_summary file, and repeat masker file exist.
if [ ! -d "${dir}" ]; then echo "Dir ${dir} does not exist"; fi
if [ ! -f "${dir}/${anno_summary_file}" ]; then echo "File ${dir}/${anno_summary_file} does not exist"; fi
if [ ! -f "${RepeatMasker_bed}" ]; then echo "Repeat Masker bed file ${RepeatMasker_bed} does not exist"; fi


# reformat columns in SAM file <sample_name>.anno_summary file conform to BED format
# also for the start position, make it 0-based by subtracting 1 from the start position.
# 1) from anno_overlap make bed file format = chrom, read_start, read_end, name, score, strand
# 2) intersect with repeatmasker bed file
# 3) get readID names that fall into repetitive regions
awk '{OFS="\t"} { print $1, $3-1, $4, $5, $11, $2 }' ${dir}/${anno_summary_file} | sed '1d' \
| bedtools intersect -a stdin -b ${RepeatMasker_bed} -wo | cut -f 4 | sort -u > ${dir}/${intermediate_output_dir}/${anno_summary_file}.intersect.readID

# filter out miRNA species from the previous output
awk 'NR == FNR {a[$0]=$0; next}; !($1 in a) {print $0}' ${dir}/${intermediate_output_dir}/${anno_summary_file}.intersect.readID ${dir}/${anno_summary_file} > ${dir}/${anno_summary_file}.rmsk

rm ${dir}/${intermediate_output_dir}/${anno_summary_file}.intersect.readID
