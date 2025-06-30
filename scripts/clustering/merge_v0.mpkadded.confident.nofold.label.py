import pandas as pd

df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/taxonomy/merge_v0.align.mpkadded.edited.confident.dedup.taxo.tsv"
df = pd.read_csv(df_path, sep="\t", usecols=["accession","taxo"])

df["accession"] = df["accession"].str.replace('_', '-')

df.columns = ["seqID", "label"]

save_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold_cm_all//merge_v0.align.mpkadded.edited.confident.dedup.label.tsv"
df.to_csv(save_path, sep="\t", index=False)