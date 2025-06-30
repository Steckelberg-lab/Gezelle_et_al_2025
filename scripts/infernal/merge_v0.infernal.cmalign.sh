#!/bin/bash

#$ -l mem=8G,time=:10:
#$ -S /bin/bash
#$ -N merge_v0_infernal
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/infernal

RE_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign_pk"

FA_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk/merge_v0.cmsearch.manual_pk.anchor_added.filter.fa"
STO_PATH="/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/merge_v0.align.dedup.pk.sto"
CM_PATH="${RE_PATH}/merge_v0.cm"
BUILD_RE_PATH="${RE_PATH}/merge_v0.cmbuild.out"
NEW_STO_PATH="${RE_PATH}/merge_v0.align.mpkadded.sto"

mkdir -p $RE_PATH

# 1. build the cm
cmbuild $CM_PATH $STO_PATH > $BUILD_RE_PATH
# 2. align sequences to a covariance model
cmalign --mapali $STO_PATH --mapstr --noprob -o $NEW_STO_PATH $CM_PATH $FA_PATH
