#!/bin/bash

#$ -l mem=20G,time=3:00:
#$ -S /bin/bash
#$ -N nofold_score
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/NoFold

cd /ifs/data/as6282_gp/fy2306/tools/nofold/src/

python score_and_normalize.py \
/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold/merge_v0_clustering.align.dedup.clustering.NoFoldID.fasta \
--cpus=1 \
--cm-folder=/ifs/data/as6282_gp/fy2306/tools/nofold/models/rfam_cms/
