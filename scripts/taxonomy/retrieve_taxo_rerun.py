import pandas as pd
import re
from Bio import Entrez, SeqIO, AlignIO

df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/taxonomy/merge_v0.align.mpkadded.edited.confident.dedup.taxo.tsv"
df = pd.read_csv(df_path, sep="\t")


Entrez.email = "fy2306@columbia.edu"

def fetch_taxonomy_info(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    
        # Extract the taxonomy information
        taxonomy_info = record.annotations['taxonomy']
        
        family = next((taxon for taxon in taxonomy_info if taxon.endswith("viridae")), taxonomy_info[-1])
        print(f"{accession} fetched")
        return taxonomy_info, family
    except Exception as e:
        print(f"Error fetching taxonomy for {accession}: {e}")
        return [], "Unknown"
    
for index, row in df[df['taxo'] == 'Unknown'].iterrows():
    accession = row['accession']
    if accession:
        taxonomy_info, family = fetch_taxonomy_info(accession)
        df.at[index, 'taxo_long'] = taxonomy_info
        df.at[index, 'taxo'] = family

df.to_csv(df_path, sep='\t', index=False)