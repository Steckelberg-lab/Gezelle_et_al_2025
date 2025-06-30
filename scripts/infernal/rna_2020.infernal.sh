#!/bin/bash

#$ -l mem=8G,time=3:55:
#$ -S /bin/bash
#$ -N rna_2020_infernal
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/infernal

RE_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/rna_2020_green_plants"

STO_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/rna_2020_paper.edited.sto"
CM_PATH="${RE_PATH}/rna_2020.cm"
# Note: ST9 (NC_004045) is not in RNAVirus or Viridiplantae subset (possibly due to lack of metadata)
DB_PATH="/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/ncbi_virus/ncbi_virus_Viridiplantae.fa.gz"

BUILD_RE_PATH="${RE_PATH}/rna_2020.cmbuild.out"
CAL_RE_PATH="${RE_PATH}/rna_2020.cmcalibrate.out"
SEARCH_RE_PATH="${RE_PATH}/rna_2020.cmsearch.out"
SEARCH_RE_D_PATH="${RE_PATH}/rna_2020.cmsearch.detail.out"

mkdir -p $RE_PATH

# 1. build the cm
cmbuild $CM_PATH $STO_PATH > $BUILD_RE_PATH
# 2. calibrate
# cmcalibrate --forecast --nforecast 1 $CM_PATH
cmcalibrate $CM_PATH > $CAL_RE_PATH
# search
cmsearch -T 0 --tblout $SEARCH_RE_PATH -o $SEARCH_RE_D_PATH $CM_PATH $DB_PATH 