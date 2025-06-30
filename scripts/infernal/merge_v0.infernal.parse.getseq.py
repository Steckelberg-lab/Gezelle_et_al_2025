# Parse the cmsearch --tblout output
# intended for 2-step cmsearch
# get the entire seq for each result and save as fasta

import pandas as pd
from Bio import SeqIO, AlignIO
import gzip

# Config
re_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v2_2_pk/merge_v0.cmsearch.out"
db_path = "/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/ncbi_virus/ncbi_virus_All.fa.gz"
re_fa_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v2_2_pk/merge_v0.cmsearch.out.fa.gz"

# 1. parse the cmsearch output

re_df_cols = ["target_name", "accession", "query_name", "accession1", "mdl", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc", "pass", "gc", "bias", "score", "E_value", "inc", "description_of_target"]

# Initialize a list to hold all the processed rows
data = []

# Open the file and read line by line
with open(re_path, 'r') as file:
    for line in file:
        # Ignore lines that are comments
        if line.startswith('#'):
            continue
        
        # Split the line at the first occurrence of '|'
        parts = line.split('|', 1)
        
        if len(parts) == 2:
            # The part before '|' is further split by whitespace
            first_parts = parts[0].strip().split()
            # The part after '|' is taken as is, with leading and trailing whitespaces stripped
            last_part = parts[1].strip()
            # Combine the split parts and the last part
            row = first_parts + [last_part]
            data.append(row)
        else:
            # Handle cases where there is no '|' character
            data.append(parts[0].strip().split())

# Create DataFrame
re_df = pd.DataFrame(data, columns=re_df_cols)
print(len(re_df))
accession_list = re_df['target_name'].unique().tolist()
print(len(accession_list))



# 2. retrieve the seqs according to target_name
# 3. Save as fasta in re_fa_path
def filter_accessions(record, accessions):
    return any(acc in record.description.split('|')[0] for acc in accessions)

# Reading the gzip FASTA file, filtering, and writing to a new gzip FASTA file
with gzip.open(db_path, 'rt') as original_fasta, gzip.open(re_fa_path, 'wt') as filtered_fasta:
    sequences = SeqIO.parse(original_fasta, 'fasta')
    filtered_sequences = (seq for seq in sequences if filter_accessions(seq, accession_list))
    SeqIO.write(filtered_sequences, filtered_fasta, 'fasta')