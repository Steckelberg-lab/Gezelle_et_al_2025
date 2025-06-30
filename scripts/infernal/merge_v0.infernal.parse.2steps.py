# Parse the cmsearch --tblout output
# intended for 2-step cmsearch

import pandas as pd
from Bio import SeqIO, AlignIO
import gzip

meta_path = "/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/ncbi_virus/ncbi_virus_Viridiplantae.metadata.csv"
sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/merge_v0.align.sto"

E_val_thresh_step2 = 10

shared_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2"
# sl: stem-loop, pk: pseudoknot
# sl
out_1_step1 = f"{shared_path}/merge_v0.cmsearch.out"
fa_1_step1 = f"{shared_path}/merge_v0.cmsearch.out.fa.gz"
# pk
out_1_step2 = f"{shared_path}/two_steps/merge_v0.cmsearch.out"

# pk
out_2_step1 = f"{shared_path}_pk/merge_v0.cmsearch.out"
fa_2_step1 = f"{shared_path}_pk/merge_v0.cmsearch.out.fa.gz"
# sl
out_2_step2 = f"{shared_path}_pk/two_steps/merge_v0.cmsearch.out"

final_re_path = f"{shared_path}/merge_v0.cmsearch.2steps.out"


# load in the cmsearch tbl out
def load_df(re_path):
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
	re_df['E_value'] = pd.to_numeric(re_df['E_value'], errors='raise')
	re_df['seq_from'] = pd.to_numeric(re_df['seq_from'], errors='raise')
	re_df['seq_to'] = pd.to_numeric(re_df['seq_to'], errors='raise')
	return re_df

# check if two results overlap
def check_overlap(row1, row2, stype1, stype2):
    pk_row = row1 if stype1 == 'pk' else row2
    sl_row = row1 if stype1 == 'sl' else row2

    # Calculate differences based on pk and sl
    seq_from_diff = (pk_row['seq_from'] - sl_row['seq_from']) if row1["strand"]=="+" else (sl_row['seq_from'] - pk_row['seq_from'])
    seq_to_diff = (pk_row['seq_to'] - sl_row['seq_to']) if row1["strand"]=="+" else (sl_row['seq_to'] - pk_row['seq_to'])

    # Check if both differences satisfy the conditions
	# add abs?
    if 0 <= abs(seq_from_diff) <= 20 and 0 <= abs(seq_to_diff) <= 20:
        return True
    return False


# retrieve the seq from fasta
def extract_sequences(fasta_path, df):
    # Load the compressed FASTA file
    with gzip.open(fasta_path, "rt") as handle:
        # Create a dictionary from the FASTA file
         seq_dict = {record.id.split(' |')[0]: record for record in SeqIO.parse(handle, "fasta")}

    # Function to fetch sequence based on row in DataFrame
    def fetch_sequence(row):
        # Extract relevant sequence information
        seq_id = row['target_name']
        start = min(int(row['seq_from'])-1, int(row['seq_to'])-1)
        end = max(int(row['seq_from']), int(row['seq_to']))
        strand = row['strand']

        # Retrieve sequence and apply range and strand
        try:
            seq = seq_dict[seq_id].seq[start:end]
            if strand == '-':
                seq = seq.reverse_complement()
            seq = str(seq).upper()
            seq = seq.replace('T', 'U')
            return seq
        except KeyError:
            return "Sequence not found" 

    # Apply function across DataFrame rows
    df['sequence'] = df.apply(fetch_sequence, axis=1)
    return df


# main function
def steps2seq(step1_out_path, step2_out_path, step1_fa_path, step1_stype, step2_stype):

	# 1. load 2 steps df
	step1_df = load_df(step1_out_path)
	step2_df = load_df(step2_out_path)

	# 2. filters
	# step1_df: no filters
	# step2_df: E-value <= 10
	step2_df = step2_df[step2_df["E_value"] <= E_val_thresh_step2]
	# sort based on E_value (should be the default)
	step1_df.sort_values(by='E_value', ascending=True, inplace=True)
	step1_df.sort_values(by='E_value', ascending=True, inplace=True)
	print(f"After E-value filtering, # of seqs in step1_df: {len(step1_df)}")
	print(f"After E-value filtering, # of seqs in step2_df: {len(step2_df)}")

	# 3. overlap of two steps
	results = []

	for index2, row2 in step2_df.iterrows():
		matches = step1_df[(step1_df['target_name'] == row2['target_name']) & (step1_df['strand'] == row2['strand'])]
		for index1, row1 in matches.iterrows():
			if check_overlap(row1, row2, step1_stype, step2_stype):
				sl_row = row2 if step2_stype == 'sl' else row1
				pk_row = row1 if step1_stype == 'pk' else row2
				result = {
					'target_name': row1['target_name'],
					'seq_from': min(sl_row['seq_from'], pk_row['seq_from']) if row1["strand"] == "+" else max(sl_row['seq_from'], pk_row['seq_from']),
					'seq_to': max(pk_row['seq_to'], sl_row['seq_to']) if row1["strand"] == "+" else min(pk_row['seq_to'], sl_row['seq_to']),
					'strand': row1['strand']
				}
				results.append(result)
				break

	steps_df = pd.DataFrame(results)
	print(f"Overlapping results in 2 steps: {len(steps_df)}")


	# 4. retrieve the seqs
	steps_df = extract_sequences(step1_fa_path, steps_df)
 
	# 5. optional: dedup
	steps_df = steps_df.drop_duplicates(subset='sequence', keep='first')
	print(f"Deduplicated overlapping results in 2 steps: {len(steps_df)}")

	return steps_df

print("Dealing with sl-pk")
df_1 = steps2seq(out_1_step1, out_1_step2, fa_1_step1, "sl", "pk")
print(f"df_1 row #: {len(df_1)}")
print("-----------------------------------------------------------------------------------")
print("Dealing with pk-sl")
df_2 = steps2seq(out_2_step1, out_2_step2, fa_2_step1, "pk", "sl")
print(f"df_2 row #: {len(df_2)}")

# merge
re_df = pd.concat([df_1, df_2], ignore_index=True)

# optional: dedup
re_df = re_df.drop_duplicates(subset='sequence', keep='first')
# re_df = re_df.drop_duplicates(subset=None, keep='first')
print(f"Deduplication on the merged df: {len(re_df)}")

# optional: dedup on sto file
sto = AlignIO.read(sto_path, 'stockholm')
sto_data = [{'seq_name': record.id, 'sequence': str(record.seq).replace('-', '')} for record in sto]
sto_df = pd.DataFrame(sto_data)
re_df = re_df[~re_df["sequence"].isin(sto_df["sequence"].tolist())]
print(f"Deduplication based on sto input seqs: {len(re_df)}")

# add metadata
meta_df = pd.read_csv(meta_path, sep=",", usecols=['Accession', 'Organism_Name', 'Species', 'Molecule_type'])
meta_df.columns = ["target_name", "Organism_Name", "Species", "Molecule_type"]
re_df = pd.merge(re_df, meta_df, on='target_name', how='inner')
print(f"Final seq #: {len(re_df)}")

# # optional: filter all that occur in green plants
# gp_path = "/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/ncbi_virus/ncbi_virus_Viridiplantae.metadata.csv"
# gp_df = pd.read_csv(gp_path, sep=",", usecols=['Accession'])
# gp_accession = gp_df["Accession"].tolist()
# re_df = re_df[~re_df['target_name'].isin(gp_accession)]
# print(f"Final seq #: {len(re_df)}")
# print(f"Final species #: {len(re_df.drop_duplicates(subset=["sequence", "Species"], keep="first"))}")

re_df.to_csv(final_re_path, sep="\t", index=False)