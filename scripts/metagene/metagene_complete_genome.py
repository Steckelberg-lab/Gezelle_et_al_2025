import pandas as pd
from Bio import Entrez, SeqIO

df_path = "/path/to/merge_v0.align.mpkadded.edited.confident.dedup.info.tsv"
df_save_path = "/path/to/merge_v0.align.mpkadded.edited.confident.dedup.complete_genome.tsv"
df = pd.read_csv(df_path, sep="\t")

# complete genome
df = df[(df["complete"] == "complete") & (df["genome"] == "genome")]

Entrez.email = "name@email.com"

def fetch_genome_features(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    features = []
    cds_features = []
    for feature in record.features:
        if feature.type != "gene" and feature.type != "source":
            features.append({
                "type": feature.type,
                "start": int(feature.location.start),
                "end": int(feature.location.end)
            })
        if feature.type == "CDS":
            cds_features.append({
                "type": feature.type,
                "start": int(feature.location.start),
                "end": int(feature.location.end)
            })
    
    try:
        min_cds_start = min(f["start"] for f in features if f["type"] == "CDS")
        max_cds_end = max(f["end"] for f in features if f["type"] == "CDS")

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
            
    except:
        min_cds_start = 0
        max_cds_end = len(record.seq)
    
    features.append({"type": "GenomeStart", "start": 0, "end": min_cds_start - 1})
    features.append({"type": "GenomeEnd", "start": max_cds_end + 1, "end": len(record.seq)})

    return features

def determine_feature(seq_start, seq_end, features):
    overlapping_features = []
    # Note! overlapping definition
    for feature in features:
        if feature["start"] <= seq_start <= feature["end"] or feature["end"] >= seq_end >= feature["start"] or seq_start <= feature["start"] <= seq_end or seq_start <= feature["end"] <= seq_end:
            overlapping_features.append(feature["type"])
    
    if not overlapping_features:
        return "unknown"
    return ",".join(overlapping_features)


results = []
for index, row in df.iterrows():
    accession = row["accession"]
    seq_from = row["seq_from"]
    seq_to = row["seq_to"]
    
    features = fetch_genome_features(accession)
    feature_info = determine_feature(seq_from, seq_to, features)
    
    results.append(feature_info)
    print(f"{accession} + {feature_info}")

df["feature"] = results

df.to_csv(df_save_path, sep='\t', index=False)