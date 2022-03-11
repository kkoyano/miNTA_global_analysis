#!/bin/bash

# this is the full pipeline for identifying end modification in small RNAseq
# Usage:
# ./end_modificication_pipeline-FULL.V4.sh Sample_1S12 mapping_opt_unique.txt /u/nobackup/gxxiao2/apps/genomes/hg19/bowtieindex/all.chr.fa annotation_paths_mir_pir_trf_rfam.txt priority_file.txt 1 0.5 /home/data/plasma/disease /home/data/plasma/raw 30 /home/data/plasma/raw


# BOWTIE MAPPING PARAMETERS
# cat mapping_opt_unique.txt
    # bowtie -v 1 -m 1 --best --strata -S --sam-nohead

# SMALL NON-CODING RNA ANNOTATIONS
# cat annotation_paths_mir_pir_trf_rfam.txt
    # /u/home/k/kkoyano/nobackup-gxxiao3/annotations/mirbase/hg19.gff3.mature.mirs-for_parse.bowtie.py.txt
    # /u/home/k/kkoyano/nobackup-gxxiao3/annotations/pirnabank/hsa_piRNABank_hg18.to.hg19liftover.txt
    # /u/home/k/kkoyano/nobackup-gxxiao3/annotations/rna.genes/reformat.rfam.human.nomirna.notrna.txt
    # /u/home/k/kkoyano/nobackup-gxxiao3/annotations/hg19/tRfdb/tRfs_all_annotation.bed

# Usage:
# test
# sample_name=T_lymphocyte
# id_abundance_ext=noadapter.id.abundance.fa
# id_abundance_path=/u/home/k/kkoyano/nobackup-gxxiao3/data/cell_type_srna/sra/id.abundance.fa
# mapping_parameter_file=mapping_opt_unique.txt
# output_dir=/u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/unique_reads/data/cell_type
# annotation_paths_file=/u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/scripts/annotation_paths_mir_pir_trf_rfam.txt
# bowtieindex_path=/u/nobackup/gxxiao2/apps/genomes/hg19/bowtieindex/all.chr.fa
# PriorityAnno_file=/u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/scripts/priority_file.txt


# input parameters
sample_name=$1


mapping_parameter_file=$2
bowtieindex_path=$3
annotation_paths_file=$4
PriorityAnno_file=$5
chrom_fa_path=/u/nobackup/gxxiao2/apps/genomes/hg19/CHR
level=$6
canonical_threshold=$7
output_dir=$8
id_abundance_path=$9
id_abundance_ext=id.abundance.fa
snp_dir_path=/u/home/k/kkoyano/nobackup-gxxiao3/annotations/hg19/snps/hg19.comprehensive.7.14.2017/hg19.comprehensive.7.14.2017.snp.
min_AF=100
max_AF=100
phred_threshold=${10}
fq_path=${11}
fq_extension=fastq.gz

# step 1: Main script to identify end modifications. fa file --> sam file with known NTEs

python /u/nobackup/gxxiao3/kkoyano/scripts/trimming_pipeline/NTE_trimming_pipeline.V4.py -i ${id_abundance_path}/${sample_name}.${id_abundance_ext} --fasta -o ${output_dir}/${sample_name} -b ${mapping_parameter_file} -x ${bowtieindex_path}
echo "STEP 1 finished"

if [[ ! -f ${output_dir}/${sample_name}.NTE.sam ]] ; then
    echo "File ${output_dir}/${sample_name}.NTE.sam is not there, aborting."
    exit
fi

# step 2: Annotating alignments from the sam file
python /u/home/k/kkoyano/nobackup-gxxiao3/scripts/get.reads.overlap.long_annotation.V9.trimming.py -s ${output_dir}/${sample_name}.NTE.sam -a ${annotation_paths_file} -PriorityAnno ${PriorityAnno_file} --sam
echo "STEP 2 finished"


if [[ ! -f ${output_dir}/${sample_name}.NTE.sam.anno_overlap ]] ; then
    echo "File ${output_dir}/${sample_name}.NTE.sam.anno_overlap is not there, aborting."
    exit
fi
# step3: Remove alignments that overlap with repeat masker
./filter_RepeatMasker.fromSAM.sh ${output_dir} ${sample_name}.NTE.sam.anno_overlap /u/home/k/kkoyano/nobackup-gxxiao3/annotations/hg19/RepeatMasker/rmsk.bed
echo "STEP 3 finished"

# step4: BLAT all mapped reads. remove any that may be incorrectly mapped, this step will may take up to 30min-1hr
mkdir -p ${output_dir}/intermediate
mkdir -p ${output_dir}/intermediate/blat
grep 'mirna' ${output_dir}/${sample_name}.NTE.sam.anno_overlap.rmsk | awk '{print">"$5"\n"$10}' > ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.all_mirna.fa
for i in {1..22} X Y M; do blat -stepSize=1 -minScore=0 -minIdentity=50 -fine ${chrom_fa_path}/chr${i}.fa ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.all_mirna.fa ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.chr${i}.psl; done
echo "STEP 4 finished"

# step5: consolidate BLAT results to only select miRNA that have same alignment as BLAT
cat ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.chr*.psl > ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl
echo "STEP 5 finished"

# step6: get the best scoring alignment and convert into a bed file
/u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/unique_reads/scripts/psl_to_bed_best_score.pl ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed
echo "STEP 6 finished"

# step7: Convert the original sam file into a bed file. Then compare the blat alignment bed file and the original bed file, and only output the alignments that are concordinate between the two alignment methods.
python /u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/unique_reads/scripts/NTE_sam_to_bed.py -b ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed -s ${output_dir}/${sample_name}.NTE.sam.anno_overlap.rmsk
echo "STEP 7 finished"

# step 8: Apply Blat stringency filter, and only keep alignments that match exactly, the level input parameter indicates the type of stringency is used. Level 1 by default.
python /u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/unique_reads/scripts/apply_BLAT_stringency.py -level ${level} -s ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed.results
echo "STEP 8 finished"

# step 9: for each species, select miRNA reads that have the canonical end position
python /u/home/k/kkoyano/nobackup-gxxiao3/scripts/mismatch_profile/get_canonical_species.py -i ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed.results.lvl${level} -s ${output_dir}/${sample_name}.NTE.sam.summary -t ${canonical_threshold}
echo "STEP 9 finished"

# step 10: convert SAM file (rows = readID|abundance) --> annotation summary file (rows = miRNA species), this summarizes the results of all the reads mapping to a miRNA species
python /u/nobackup/gxxiao3/kkoyano/scripts/mismatch_profile/make_anno_overlap_summary_ntes.py -anno ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed.results.lvl${level}.t${canonical_threshold}.canonical_end -summary ${output_dir}/${sample_name}.NTE.sam.summary
echo "STEP 10 finished"

# step 11: remove annotated SNPs
suffix=NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed.results.lvl${level}.t${canonical_threshold}.canonical_end.anno_summary

head -n 1 ${output_dir}/intermediate/blat/${sample_name}.${suffix} > ${output_dir}/intermediate/blat/header.${suffix}
sort -k 3,3 ${output_dir}/intermediate/blat/${sample_name}.${suffix} > ${output_dir}/intermediate/blat/${sample_name}.${suffix}.temp
cat ${output_dir}/intermediate/blat/header.${suffix} ${output_dir}/intermediate/blat/${sample_name}.${suffix}.temp > ${output_dir}/intermediate/blat/${sample_name}.${suffix}.temp2
rm ${output_dir}/intermediate/blat/${sample_name}.${suffix}.temp
python /u/nobackup/gxxiao3/kkoyano/scripts/mismatch_profile/filter_nte_snps.py -anno_summary ${output_dir}/intermediate/blat/${sample_name}.${suffix}.temp2 -o ${output_dir}/intermediate/blat/${sample_name}.${suffix}.nosnp -SNP_DIR ${snp_dir_path}
echo "STEP 11 finished"

# step 12: remove homozygous modifications SNPs
python /u/home/k/kkoyano/nobackup-gxxiao3/nte_global_analysis_2017-09/unique_reads/scripts/filter_unannotated_SNPs.py -i ${output_dir}/intermediate/blat/${sample_name}.${suffix}.nosnp -m ${min_AF} -M ${max_AF}
echo "STEP 12 finished"

# step 13: remove modifications that have low phred score
python /u/home/k/kkoyano/nobackup-gxxiao3/scripts/mismatch_profile/filter_phred_score.V2.py -a ${output_dir}/intermediate/blat/${sample_name}.${suffix}.nosnp.snv-m${min_AF}M${max_AF} -s ${output_dir}/intermediate/blat/${sample_name}.NTE.sam.anno_overlap.rmsk.mirna.ALL.psl.mirna.ALL.psl.best_score.bed.results.lvl${level} -q ${fq_path}/${sample_name}.${fq_extension} -phred ${phred_threshold}
echo "STEP 13 finished"

# step 14: get credible modifications (this relies on all the datasets being finished), requires each modification to occur in at least 5 samples of the same fluid, and also to have at least 2 read counts for the sample)

# step 15: for all samples, get the 3'NTE composition and also the 5'NTE composition
# ==> see end_modification_pipeline-FULL.V4.pt2.sh


