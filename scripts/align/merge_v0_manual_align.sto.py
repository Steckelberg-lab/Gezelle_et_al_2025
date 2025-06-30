import pandas as pd

# load in the alignment tsv
align_path = "/path/to/merge_v0.align.3.tsv"
df = pd.read_csv(align_path, delimiter='\t')

# merge seq
def check_sequence_match(row):
    columns_to_merge = [
        'end_5', 'P1_1', 'J1_2', 'P2_1', 'L2_1', 'L2_pk', 'L2_2', 
        'P2_2_1', 'P2_2_L', 'P2_2_2', 'J2_1', 'P1_2', 'J1_3', 
        'J1_3_3nt', 'J1_3_pk', 'end_3'
    ]

    merged_string = ''.join(row[col] for col in columns_to_merge if not pd.isna(row[col]))
    merged_string = merged_string.replace('.', '-')
    
    return merged_string

df['sto_seq'] = df.apply(check_sequence_match, axis=1)

# write
with open('/path/to/merge_v0.align.sto', 'w') as file:
    for index, row in df.iterrows():
        file.write(f"{row['seqname']} {row['sto_seq']}\n")
