#!/bin/bash

#$ -l mem=2G,time=2:00:
#$ -S /bin/bash
#$ -N retrieve_taxo
#$ -cwd
#$ -t 1
#$ -o /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.out
#$ -e /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/logs/$JOB_NAME.err

source /ifs/data/as6282_gp/fy2306/miniconda3/bin/activate /ifs/data/as6282_gp/fy2306/miniconda3/envs/py3.12

python /ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/script/merge_v0/metagene/accession_info.py