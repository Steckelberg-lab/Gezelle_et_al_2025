import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict

# Define paths
sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.dedup.sto"
df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/taxonomy/merge_v0.align.mpkadded.edited.confident.dedup.taxo.tsv"

df = pd.read_csv(df_path, sep='\t')

alignment = AlignIO.read(sto_path, "stockholm")

sequences_by_taxo = defaultdict(list)

# quick lookup of taxo by target_name
taxo_dict = df.set_index('seqID')['taxo'].to_dict()

# Classify sequences based on taxo
for record in alignment:
    target_name = record.id
    if target_name in taxo_dict:
        taxo = taxo_dict[target_name]
        sequences_by_taxo[taxo].append(record)
print("Taxo collected")

# Write classified sequences into separate STO files
for taxo, records in sequences_by_taxo.items():
    output_path = f"/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/taxonomy/merge_v0.align.mpkadded.edited.confident.dedup.{taxo}.sto"
    alignment = MultipleSeqAlignment(records)
    with open(output_path, "w") as output_handle:
        AlignIO.write(alignment, output_handle, "stockholm")
    print(f"{taxo} written.")

print("Sequences have been classified and written to separate STO files.")
