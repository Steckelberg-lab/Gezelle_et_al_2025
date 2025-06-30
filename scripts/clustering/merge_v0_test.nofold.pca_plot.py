import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_label_path = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold_cm_all/merge_v0.align.mpkadded.edited.confident.dedup.label.tsv'
df_pca_path = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold_cm_all/merge_v0.bitscore.edited'

df_label = pd.read_csv(df_label_path, sep='\t')
print(df_label.head())
df_pca = pd.read_csv(df_pca_path, sep='\t', index_col=None)
print(df_pca.head())

df_merged = pd.merge(df_pca, df_label, on='seqID')

# Create a scatter plot using seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(x='PC1', y='PC2', hue='label', data=df_merged, palette=sns.color_palette("Set2", n_colors=len(df_merged['label'].unique())), s=15, alpha=0.8)
plt.title('PCA (NoFold)')
plt.xlabel('CM built on stem-loop structure')
plt.ylabel('CM built on pseudoknot structure')
plt.legend(title='Label')
plt.grid(False)

plt.savefig('/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold_cm_all/pca_plot.png')