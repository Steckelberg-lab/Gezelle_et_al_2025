# Duplication check
# intended for R-scape

from Bio import AlignIO
import pandas as pd

sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.dedup.sto"

sto = AlignIO.read(sto_path, 'stockholm')
sto_data = [{'seq_name': record.id, 'sequence': str(record.seq).replace('-', '')} for record in sto]
sto_df = pd.DataFrame(sto_data)

duplicates = sto_df.groupby('sequence')['seq_name'].agg(list)
duplicates = duplicates[duplicates.apply(len) > 1]

for _, names in duplicates.items():
    print(f"Duplicate sequence found in sequences: {', '.join(names)}")
    
print(len(sto_df))