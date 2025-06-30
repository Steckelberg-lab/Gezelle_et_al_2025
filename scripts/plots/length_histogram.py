import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO

sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/plots/merge_v0.align.mpkadded.edited.confident.dedup.length.sto"

alignment = AlignIO.read(sto_path, 'stockholm')

sequences = []
for record in alignment:
    seq = str(record.seq).replace('-', '')
    sequences.append(seq)
    if len(seq) >= 15:
        print(f"{record.id}: {len(seq)}")

seq_lengths = [len(seq) for seq in sequences]

df = pd.DataFrame({'sequence_length': seq_lengths})
df.to_csv("/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/plots/merge_v0.align.mpkadded.edited.confident.dedup.length.tsv", sep="\t", index=False)

min_length = df['sequence_length'].min()
max_length = df['sequence_length'].max()
all_lengths = pd.DataFrame({'sequence_length': range(min_length, max_length + 1)})

df = all_lengths.merge(df.value_counts('sequence_length').reset_index(name='count'), on='sequence_length', how='left').fillna(0)

sns.set(style='ticks')
color = sns.color_palette("Set2")[0]

plt.figure(figsize=(10, 6))
sns.barplot(x='sequence_length', y='count', data=df, color=color)

plt.xlabel('Sequence Length')
plt.ylabel('Frequency')
plt.title('Distribution of Sequence Lengths')

output_file = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/plots/sequence_length_histogram.svg'
plt.savefig(output_file, format='svg', dpi=600)

sub_df = df[df["sequence_length"] >= 15]
print(sub_df)