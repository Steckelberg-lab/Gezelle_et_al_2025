import pandas as pd
from Bio import Entrez, SeqIO


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
    
def length(row, gene_structure):
    length = []
    start = row["seq_from"]
    end = row["seq_to"]
    # Note! overlapping definition
    for feature in gene_structure:
        if feature["start"] <= start <= feature["end"] and feature["end"] >= end >= feature["start"]:
            if feature["type"] == "intergenic":
                length.append(("intergenic", (abs(feature["end"]-feature["start"])+1)))
                break
            
    return length
    
df["intergenic_length"] = df.apply(lambda row: length(row, gene_structures[row["accession"]]), axis=1)

df.to_csv("/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/merge_v0.align.mpkadded.edited.confident.dedup.complete_genome.inter_length.tsv", sep="\t", index=False)
