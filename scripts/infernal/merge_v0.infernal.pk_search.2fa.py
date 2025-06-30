# Following the pk search
# Filter the results and generate a fasta for alignment

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

out_path = "/path/to/merge_v0.cmsearch.manual_pk.anchor_added.tsv"
re_fa_path = "/path/to/merge_v0.cmsearch.manual_pk.anchor_added.filter.fa"

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
    header = f"{row['target_name']}.{row['seq_from']}-{row['seq_to']}.{row['strand']}.infernal"
    return SeqRecord(Seq(row['sequence']), id=header, description='')

seq_records = filtered_df.apply(make_fasta_record, axis=1).tolist()
SeqIO.write(seq_records, re_fa_path, "fasta")
