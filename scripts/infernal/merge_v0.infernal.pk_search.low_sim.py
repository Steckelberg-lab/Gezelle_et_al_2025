# following manual pk search - cmalign step
# take advantage of cmalign results and filter for low-similarity seqs
# Note: This step is not considered as necessary anymore.

import pandas as pd
from Bio import AlignIO

def low_sim(sto_path):
	# load the sto file
	sto = AlignIO.read(sto_path, 'stockholm')
	sto_data = [{'seq_name': record.id, 'sequence': str(record.seq).replace('.', '-')} for record in sto]	# uniform insertions
	sto_df = pd.DataFrame(sto_data)

	# random shuffle
	# because we do pairwise comparison in the order of row index
	# for highly similar row1 and row2, if we remove row1 first, we don't have to remove row2
	# whether a row, highly similar to one or more rows, is removed is related to its index
	# seed = 1
	# sto_df = sto_df.sample(frac=1, random_state=seed)
	# sto_df = sto_df.reset_index(drop=True)
      
	def calculate_dissimilarity(seq1, seq2):
		return sum(1 for a, b in zip(seq1, seq2) if a.upper() != b.upper())

	rows_to_remove = []

	for i, row in sto_df.iterrows():
		sequence = row['sequence']
		dissimilarities = [
			calculate_dissimilarity(sequence, other_row['sequence'])
			for j, other_row in sto_df.iterrows() if i != j
		]
		if min(dissimilarities) < 1:
			rows_to_remove.append(i)

	sto_df = sto_df.drop(rows_to_remove).reset_index(drop=True)

	return sto_df

def write_new_sto(sto_path, sto_df, new_sto_path):
	with open(new_sto_path, 'w') as file:
		file.write("# STOCKHOLM 1.0\n")
		file.write(f"#=GF SQ {len(sto_df)}\n")
		for _, row in sto_df.iterrows():
			file.write(f"{row['seq_name']} {row['sequence']}\n")
		with open(sto_path, 'r') as original_sto_file:
			for line in original_sto_file:
				if line.startswith('#=GC SS_cons'):
					file.write(line)
					break
		file.write("//")

sto_path_1 = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.sto"
sto_df_1 = low_sim(sto_path_1)
print(len(sto_df_1))

new_sto_path_1 = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.low_sim.sto"
write_new_sto(sto_path_1, sto_df_1, new_sto_path_1)