#!/bin/bash

dir=$1
sam_filename=$2
fluid=$3
annoSummary_suffix=$4
sam_suffix=.anno_overlap.t0.5.canonical_end.nte_fields.mirna.ALL.psl.best_score.bed.results.lvl1

jobarray_name=filter_credible_modifications.${fluid}.jobarray
credible_input_fn=input_file_${fluid}.txt
output_group_fn=credible_modifications_${fluid}.txt

min_sample_modification_coverage=$5
min_sample_occurance=$6

# min_sample_modification_coverage=2
# min_sample_occurance=2
# min_sample_occurance=5

# if the col1 and col2 files already exist, remove them to start a new input file.
if [ -f ${dir}/col1_anno_summary.txt ]; then
    rm ${dir}/col1_anno_summary.txt
fi

if [ -f ${dir}/col2_anno_overlap.txt ]; then
    rm ${dir}/col2_anno_overlap.txt
fi

while read sample_name; do echo ${dir}/${sample_name}${annoSummary_suffix} >> ${dir}/col1_anno_summary.txt; echo ${dir}/${sample_name}${sam_suffix} >> ${dir}/col2_anno_overlap.txt; done < ${sam_filename}

paste ${dir}/col1_anno_summary.txt ${dir}/col2_anno_overlap.txt > ${dir}/${credible_input_fn}
python /u/home/k/kkoyano/nobackup-gxxiao3/scripts/mismatch_profile/filter_credible_modifications.py -i ${dir}/${credible_input_fn} -r ${min_sample_modification_coverage} -t ${min_sample_occurance} -o ${dir}/${output_group_fn} --output_anno_summary

