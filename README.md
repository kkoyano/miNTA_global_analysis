# miNTA_global_analysis
all scripts used for the project "nte_global_analysis_2017-09" 

# Owner: Kiku Koyano
# email: kkoyano@ucla.edu
# date: 08/30/2018

# Description: packet of scripts for Ling. All together, the script will start from raw fastq files, convert them into FASTA format so that the readID is the name|abundance, this is commonly used for small RNAseq to consolidate the sequencing files since many of the same sequences will appear more than once.



# step1, run the pipeline for each sample. Here I just gave you all the files for one sample --> Sample_1S12
# This pipeline will map the FASTA reads and then go through multiple steps to identify 3'-end, 5'-end and internal modifications. Multiple steps are performed in order to remove false positive modifications arising from mapping errors, sequencing errors or SNPs.
./end_modificication_pipeline-FULL.V4.sh Sample_1S12 mapping_opt_unique.txt /u/nobackup/gxxiao2/apps/genomes/hg19/bowtieindex/all.chr.fa annotation_paths_mir_pir_trf_rfam.txt priority_file.txt 1 0.5 /home/data/plasma/disease /home/data/plasma/raw 30 /home/data/plasma/raw

# at this step we will further filter for only credible modifications. these are modifications that are trustable because they are in at least 5 samples (all modifications within samples from a group are searched), and the modification must have at least 2 read coverage in the sample. The number of samples and read coverage are user defined, but in our analysis we use 5 samples and 2 read coverage.

./get_credible_modifications-2.sh ${blat_dir} ${nte_sam_filenames} ${fluid} .${suffix} ${min_read_modification_coverage} ${min_sample_occurance}"


# 4. Getting the 3',  5' and internal modifcation profiles ${fluid} in ${dir}
python /u/nobackup/gxxiao3/kkoyano/scripts/mismatch_profile/summarize_mir_ntes.py -i ${blat_dir}/${f}.${suffix}.credible-r${min_read_modification_coverage}-N${min_sample_occurance} -o ${blat_dir}/nte_levels/

# make internal modification composition files, input is the anno_summary file.
python /u/nobackup/gxxiao3/kkoyano/scripts/mismatch_profile/internal_mutations_profile.py -i ${blat_dir}/${f}.${suffix}.credible-r${min_read_modification_coverage}-N${min_sample_occurance} -o ${blat_dir}/nte_levels/

# other files to transfer 
# GenomeFetch.py 
# helper_functions.py (in the trimming pipeline) 
# /Users/kiku/mount/scripts/mismatch_profile/helper_functions_mismatch_profile.py
# /Users/kiku/mount/scripts/mismatch_profile/helper_functions.py
#
