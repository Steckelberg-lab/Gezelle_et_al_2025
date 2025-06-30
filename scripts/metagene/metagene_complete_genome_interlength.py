import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ast

df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/merge_v0.align.mpkadded.edited.confident.dedup.complete_genome.inter_length.tsv"

df = pd.read_csv(df_path, sep="\t")

def parse_intergenic_length(intergenic_str):
    return ast.literal_eval(intergenic_str)

df['intergenic_length'] = df['intergenic_length'].apply(parse_intergenic_length)

def extract_lengths(intergenic_data):
    return [length for feature, length in intergenic_data]

df['lengths'] = df['intergenic_length'].apply(extract_lengths)
lengths_series = df['lengths'].explode().astype(int)

plt.figure(figsize=(10, 6))
sns.histplot(lengths_series, bins=100, kde=False)
plt.xlabel('Intergenic Region Length')
plt.ylabel('Frequency')
plt.title('Distribution of Intergenic Lengths')
plt.grid(False)

output_path = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/metagene/complete_genome_intergenic_length.png'
plt.savefig(output_path, dpi=300)
