# Following the manual pk search
# Filter the results and generate a fasta for alignment

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

out_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk/merge_v0.cmsearch.manual_pk.anchor_added.tsv"
re_fa_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk/merge_v0.cmsearch.manual_pk.anchor_added.filter.fa"

df = pd.read_csv(out_path, sep="\t")
# sort & dedup
df = df.sort_values(by='E_value', ascending=True)
df = df.drop_duplicates(subset="sequence", keep="first")

# filter on e-value
print(len(df))
filtered_df = df[df['E_value'] <= 10]
print(len(filtered_df))

# filter on seq length (to make sure there's stem-loop)
filtered_df['sequence_length'] = filtered_df['sequence'].apply(len)
filtered_df = filtered_df[filtered_df['sequence_length'] >= 40]
print(len(filtered_df))

def make_fasta_record(row):
    # TODO: check seq_to. Currevtly it's the seq_to of infernal results, instead of seq_to after manual pk search and trimming.
    header = f"{row['target_name']}.{row['seq_from']}-{row['seq_to']}.{row['strand']}.infernal"
    return SeqRecord(Seq(row['sequence']), id=header, description='')

seq_records = filtered_df.apply(make_fasta_record, axis=1).tolist()
SeqIO.write(seq_records, re_fa_path, "fasta")
