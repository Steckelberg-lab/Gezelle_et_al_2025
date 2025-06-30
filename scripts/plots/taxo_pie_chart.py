import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/taxonomy/merge_v0.align.mpkadded.edited.confident.dedup.taxo.tsv"
df = pd.read_csv(df_path, sep="\t")

taxo_counts = df['taxo'].value_counts()

taxo_counts = taxo_counts.sort_values(ascending=False)
total = taxo_counts.sum()
percentages = [(count / total) * 100 for count in taxo_counts]

# Create labels with percentages
labels = [f"{taxo_counts.index[i]} ({percentages[i]:.1f}%)" for i in range(len(taxo_counts))]

# Set Seaborn style
sns.set(style="whitegrid")

# Custom colors using Seaborn palette
colors = sns.color_palette("Set2")[0:len(taxo_counts)]

# Create a pie chart with custom colors without labels
plt.figure(figsize=(6, 8))
patches, texts = plt.pie(
    taxo_counts, 
    labels=None, 
    startangle=140, 
    colors=colors,
	radius=0.2
)
plt.axis('equal')  # Equal aspect ratio ensures that the pie is drawn as a circle.

# Add a legend below the pie chart
plt.legend(patches, labels, loc='upper center', bbox_to_anchor=(0.5, 0.12), ncol=2, frameon=False)

# Save the plot to a file
output_file = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/plots/taxo_pie_chart.svg'
plt.savefig(output_file, format='svg', dpi=600)