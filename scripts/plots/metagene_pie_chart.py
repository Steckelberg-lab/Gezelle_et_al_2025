import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_path = "/path/to/merge_v0.align.mpkadded.edited.confident.dedup.info.region.tsv"
df = pd.read_csv(df_path, sep="\t")

metagene_dict = {
	"3'UTR,GenomeEnd": "3'UTR",
	"CDS,CDS": "CDS",
	"CDS,CDS,intergenic": "CDS",
	"CDS,intergenic": "CDS",
	"misc_feature,intergenic": "intergenic",
	"intergenic": "intergenic",
	"CDS": "CDS",
	"GenomeEnd": "3'UTR",
	"GenomeStart": "5'UTR",
	"unknown": "unknown"
}

df['feature'] = df['feature'].map(metagene_dict).fillna(df['feature'])

feature_counts = df['feature'].value_counts()
print(feature_counts)

feature_counts = feature_counts.sort_values(ascending=False)
total = feature_counts.sum()
percentages = [(count / total) * 100 for count in feature_counts]

# Create labels with percentages
labels = [f"{feature_counts.index[i]} ({percentages[i]:.1f}%)" for i in range(len(feature_counts))]

# Set Seaborn style
sns.set(style="whitegrid")

# Custom colors using Seaborn palette
colors = sns.color_palette("Set2")[0:len(feature_counts)]

# Create a pie chart with custom colors without labels
plt.figure(figsize=(6, 8))
patches, texts = plt.pie(
    feature_counts, 
    labels=None, 
    startangle=140, 
    colors=colors,
	radius=0.2
)
plt.axis('equal')  # Equal aspect ratio ensures that the pie is drawn as a circle.

# Add a legend below the pie chart
plt.legend(patches, labels, loc='upper center', bbox_to_anchor=(0.5, 0.12), ncol=2, frameon=False)

# Save the plot to a file
output_file = '/path/to/metagene_all_pie_chart.svg'
plt.savefig(output_file, format='svg', dpi=600)