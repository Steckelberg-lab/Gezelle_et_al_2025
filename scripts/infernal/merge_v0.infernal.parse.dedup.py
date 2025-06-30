# Parse the cmsearch --tblout output
# intended for 1-step cmsearch
# get the exact seq and dedup

import pandas as pd
from Bio import SeqIO, AlignIO
import gzip

# Config
re_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2_pk/merge_v0.cmsearch.out"
interm_re_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2_pk/merge_v0.cmsearch.parsed.undedup.tsv"
final_re_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2_pk/merge_v0.cmsearch.parsed.dedup.tsv"
db_path = "/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/ncbi_virus/ncbi_virus_Viridiplantae.fa.gz"
db_meta_path = "/ifs/data/as6282_gp/fy2306/projects/xrRNA_search/data/ncbi_virus/ncbi_virus_Viridiplantae.metadata.csv"
sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/merge_v0.align.pk.sto"



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




# 2. retrieve the seqs according to target_name
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

re_df = extract_sequences(db_path, re_df)

# quick check
# print("LN865081.1, positive strand:")
# print(re_df[re_df['target_name']=="LN865081.1"]["sequence"].iloc[0])
# print("MF425857.1, negative strand:")
# print(re_df[re_df['target_name']=="MF425857.1"]["sequence"].iloc[0])



# 3. retrieve the metadata according to target_name
meta_df = pd.read_csv(db_meta_path, sep=",", usecols=['Accession', 'Species', 'Molecule_type'])
meta_df.columns = ["target_name", "Species", "Molecule_type"]
print(len(re_df))
re_df = pd.merge(re_df, meta_df, on='target_name', how='inner')
print(len(re_df)) # though we use inner-join, it should stay the same

# save this interm file
re_df.to_csv(interm_re_path, sep="\t", index=False)

# 4. Deduplication based on seqs
# remove the exactly same ones
re_df = re_df.drop_duplicates(subset='sequence', keep='first')
print(f"Deduplication based on seqs: {len(re_df)}")
# TODO: remove those that are subseq of others?

# 5. Deduplication based on input seqs
sto = AlignIO.read(sto_path, 'stockholm')
sto_data = [{'seq_name': record.id, 'sequence': str(record.seq).replace('-', '')} for record in sto]
sto_df = pd.DataFrame(sto_data)
re_df = re_df[~re_df["sequence"].isin(sto_df["sequence"].tolist())]
print(f"Deduplication based on sto input seqs: {len(re_df)}")

re_df.to_csv(final_re_path, sep="\t", index=False)