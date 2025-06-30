#!/bin/bash

#$ -l mem=8G,time=3:55:
#$ -S /bin/bash
#$ -N merge_v0_infernal_bvbrc
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

# NOTE: this is just one of the infernal searches done. 
# For detailed parameters of every infernal search, you can check README or the header of infernal output, where input files and parameters are specified.

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/infernal

RE_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/revision_bvbrc/infernal/bvbrc"

CM_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/merge_v0.cm"
DB_PATH="/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/BV_BRC/BVBRC_genome.20240423.fa.gz"

SEARCH_RE_PATH="${RE_PATH}/merge_v0.cmsearch.out"
SEARCH_RE_D_PATH="${RE_PATH}/merge_v0.cmsearch.detail.out"

mkdir -p $RE_PATH

# search
# reuse CM in the initial NCBI search (v1_2)
cmsearch -T 0 --tblout $SEARCH_RE_PATH -o $SEARCH_RE_D_PATH $CM_PATH $DB_PATH 
