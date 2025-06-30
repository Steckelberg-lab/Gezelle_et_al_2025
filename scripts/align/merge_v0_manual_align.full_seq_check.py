import pandas as pd

# load in the alignment tsv
align_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/merge_v0.align.3.tsv"
df = pd.read_csv(align_path, delimiter='\t')

# merge seq
def check_sequence_match(row):
    columns_to_merge = [
        'end_5', 'P1_1', 'J1_2', 'P2_1', 'L2_1', 'L2_pk', 'L2_2', 
        'P2_2_1', 'P2_2_L', 'P2_2_2', 'J2_1', 'P1_2', 'J1_3', 
        'J1_3_3nt', 'J1_3_pk', 'end_3'
    ]

    merged_string = ''.join(row[col] for col in columns_to_merge if not pd.isna(row[col]))
    merged_string = merged_string.replace('.', '')

	# check
    if merged_string.upper() != row['full_seq'].upper():
        return row['seqname']
    return None

df['mismatch_seqname'] = df.apply(check_sequence_match, axis=1)

mismatched_sequences = df[df['mismatch_seqname'].notnull()]['mismatch_seqname']

print("Sequence names with mismatches:")
print(mismatched_sequences)

