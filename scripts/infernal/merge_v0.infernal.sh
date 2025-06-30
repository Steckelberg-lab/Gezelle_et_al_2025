#!/bin/bash

#$ -l mem=8G,time=3:55:
#$ -S /bin/bash
#$ -N merge_v0_infernal
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

# NOTE: this is just one of the infernal searches done. 
# For detailed parameters of every infernal search, you can check README or the header of infernal output, where input files and parameters are specified.

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/infernal

RE_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v2_2_pk/two_steps"

STO_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/merge_v0.align.sto"
CM_PATH="${RE_PATH}/merge_v0.cm"
DB_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v2_2_pk/merge_v0.cmsearch.out.fa.gz"

BUILD_RE_PATH="${RE_PATH}/merge_v0.cmbuild.out"
CAL_RE_PATH="${RE_PATH}/merge_v0.cmcalibrate.out"
SEARCH_RE_PATH="${RE_PATH}/merge_v0.cmsearch.out"
SEARCH_RE_D_PATH="${RE_PATH}/merge_v0.cmsearch.detail.out"

mkdir -p $RE_PATH

# 1. build the cm
cmbuild $CM_PATH $STO_PATH > $BUILD_RE_PATH
# 2. calibrate
# cmcalibrate --forecast --nforecast 1 $CM_PATH
cmcalibrate $CM_PATH > $CAL_RE_PATH
# search
cmsearch -T 0 --tblout $SEARCH_RE_PATH -o $SEARCH_RE_D_PATH $CM_PATH $DB_PATH 