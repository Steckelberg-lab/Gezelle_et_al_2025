import pandas as pd
from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/merge_v0.align.mpkadded.edited.confident.dedup.info.tsv"
df_save_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/merge_v0.align.mpkadded.edited.confident.dedup.complete_cds.tsv"
df = pd.read_csv(df_path, sep="\t")

# complete cds
df = df[(df["complete"] == "complete") & (df["genome"] == "cds")]

Entrez.email = "fy2306@columbia.edu"

def fetch_genbank_record(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

def extract_gene_structure(record):
    gene_structure = {
        'cds_start': None,
        'cds_end': None,
        'utr5_start': None,
        'utr5_end': None,
        'utr3_start': None,
        'utr3_end': None
    }
    
    cds_start_positions = []
    cds_end_positions = []

    for feature in record.features:
        if feature.type == "CDS":
            cds_start_positions.append(int(feature.location.start))
            cds_end_positions.append(int(feature.location.end))
    
	# Multiple cds
    if cds_start_positions and cds_end_positions:
        gene_structure['cds_start'] = min(cds_start_positions)
        gene_structure['cds_end'] = max(cds_end_positions)

    # Define UTRs
    if gene_structure['utr5_start'] is None:
        gene_structure['utr5_start'] = 0  
        gene_structure['utr5_end'] = gene_structure['cds_start'] - 1

    if gene_structure['utr3_start'] is None:
        gene_structure['utr3_start'] = gene_structure['cds_end'] + 1
        gene_structure['utr3_end'] = len(record.seq)  

    return gene_structure


gene_structures = {}
for accession in df["accession"].unique():
    record = fetch_genbank_record(accession)
    gene_structure = extract_gene_structure(record)
    gene_structures[accession] = gene_structure

# Normalize
# Note! overlapping definition
def normalize_position(row, gene_structure):
    relative_positions = []
    start = row["seq_from"]
    end = row["seq_to"]
    
    # 5'UTR
    if start < gene_structure["cds_start"]:
        relative_start = start / gene_structure["cds_start"]
        if end < gene_structure["cds_start"]:	# entirely in 5'UTR
            relative_end = end / gene_structure["cds_start"]
            relative_positions.append(("5'UTR", relative_start, relative_end))
        else:	# overlap with CDS
            relative_positions.append(("5'UTR", relative_start, 1.0))
            start = gene_structure["cds_start"]
    
    # CDS
    if gene_structure["cds_start"] <= start <= gene_structure["cds_end"]:
        relative_start = (start - gene_structure["cds_start"]) / (gene_structure["cds_end"] - gene_structure["cds_start"])
        if end <= gene_structure["cds_end"]:	# entirely in CDS
            relative_end = (end - gene_structure["cds_start"]) / (gene_structure["cds_end"] - gene_structure["cds_start"])
            relative_positions.append(("CDS", relative_start, relative_end))
        else:	# overlap with 3'UTR
            relative_positions.append(("CDS", relative_start, 1.0))
            start = gene_structure["cds_end"] + 1
    
    # 3'UTR
    if start > gene_structure["cds_end"]:
        relative_start = (start - gene_structure["cds_end"]) / (gene_structure["utr3_end"] - gene_structure["cds_end"])
        relative_end = (end - gene_structure["cds_end"]) / (gene_structure["utr3_end"] - gene_structure["cds_end"])
        relative_positions.append(("3'UTR", relative_start, relative_end))

    return relative_positions

df["relative_positions"] = df.apply(lambda row: normalize_position(row, gene_structures[row["accession"]]), axis=1)

# Save the new DataFrame to a TSV file
df.to_csv(df_save_path, sep='\t', index=False)


# Aggregate

flat_relative_positions = []
for _, row in df.iterrows():
    for each in row["relative_positions"]:
        region, start, end = each
        flat_relative_positions.append((region, start, end))

relative_df = pd.DataFrame(flat_relative_positions, columns=["region", "start", "end"])


# Visualize

visualization_data = []

# Expand each range into individual positions
for _, row in relative_df.iterrows():
    start = row["start"]
    end = row["end"]
    region = row["region"]
    positions = np.linspace(start, end, num=int((end - start) * 100))
    visualization_data.extend([(region, position) for position in positions])

visualization_df = pd.DataFrame(visualization_data, columns=["region", "position"])

# def adjust_position(row):
#     if row["region"] == "CDS":
#         return row["position"] + 1
#     elif row["region"] == "3'UTR":
#         return row["position"] + 2
#     else:
#         return row["position"]

# # Apply the function to each row
# visualization_df["position"] = visualization_df.apply(adjust_position, axis=1)
visualization_df = visualization_df[visualization_df["region"] == "3'UTR"]

plt.figure(figsize=(8, 6))
sns.histplot(data=visualization_df, x="position", bins=100, kde=False, element="step")
# plt.axvline(x=1, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=2, color='red', linestyle='--', linewidth=1)
plt.title("Metagene Analysis")
plt.xlabel("3'UTR")
plt.ylabel("Frequency")

output_filename = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/complete_cds_3UTR.png"
plt.savefig(output_filename, dpi=300)