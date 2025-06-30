import pandas as pd
from Bio import Entrez, SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/merge_v0.align.mpkadded.edited.confident.dedup.complete_genome.tsv"
df = pd.read_csv(df_path, sep="\t")

df = df[(df["complete"] == "complete") & (df["genome"] == "genome")]
df = df[(df['feature'] == "intergenic") | (df['feature'] == "misc_feature,intergenic")]

Entrez.email = "fy2306@cumc.columbia.edu"


def fetch_genome_features(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    features = []
    cds_features = []
    for feature in record.features:
        if feature.type == "CDS":
            cds_features.append({
                "type": feature.type,
                "start": int(feature.location.start),
                "end": int(feature.location.end)
            })
    

    # Sort CDS features by start position
    cds_features.sort(key=lambda x: x["start"])
    
    # Merge overlapping CDS features
    merged_cds = []
    current_cds = cds_features[0]
    
    for next_cds in cds_features[1:]:
        if next_cds["start"] <= current_cds["end"]:
            current_cds["end"] = max(current_cds["end"], next_cds["end"])
        else:
            merged_cds.append(current_cds)
            current_cds = next_cds
    
    merged_cds.append(current_cds)
    
    # Add intergenic regions
    for i in range(len(merged_cds) - 1):
        intergenic_start = merged_cds[i]["end"] + 1
        intergenic_end = merged_cds[i + 1]["start"] - 1
        if intergenic_start <= intergenic_end:
            features.append({
                "type": "intergenic",
                "start": intergenic_start,
                "end": intergenic_end
            })

    return features

gene_structures = {}
for accession in df["accession"].unique():
    gene_structure = fetch_genome_features(accession)
    gene_structures[accession] = gene_structure
    
def normalize_position(row, gene_structure):
    relative_positions = []
    start = row["seq_from"]
    end = row["seq_to"]
    # Note! overlapping definition
    for feature in gene_structure:
        if feature["start"] <= start <= feature["end"] and feature["end"] >= end >= feature["start"]:
            relative_start = (start - feature["start"]) / (feature["end"] - feature["start"])
            relative_end = (end - feature["start"]) / (feature["end"] - feature["start"])
            relative_positions.append(("intergenic", relative_start, relative_end))
            break
            
    return relative_positions
    
df["relative_positions"] = df.apply(lambda row: normalize_position(row, gene_structures[row["accession"]]), axis=1)


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

plt.figure(figsize=(8, 6))
sns.histplot(data=visualization_df, x="position", bins=100, kde=False, element="step")
plt.title("Metagene Analysis")
plt.xlabel("complete genome - intergenic region")
plt.ylabel("Frequency")

output_filename = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/complete_genome_intergenic.png"
plt.savefig(output_filename, dpi=300)