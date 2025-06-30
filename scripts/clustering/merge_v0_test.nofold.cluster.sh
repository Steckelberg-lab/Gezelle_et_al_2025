#!/bin/bash

#$ -l mem=20G,time=3:00:
#$ -S /bin/bash
#$ -N nofold_cluster
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/NoFold

cd /ifs/data/as6282_gp/fy2306/tools/nofold/src/

python nofold_pipeline.py \
/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold/merge_v0_clustering.zNorm.pcNorm100.zNorm.bitscore \
/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold/merge_v0_clustering.align.dedup.clustering.NoFoldID.fasta \
--cpus=1 \
--bounds-file=../thresh/bounds_100seq.txt \
--verbose
