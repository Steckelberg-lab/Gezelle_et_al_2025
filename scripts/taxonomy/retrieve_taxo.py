import pandas as pd
import re
from Bio import Entrez, SeqIO, AlignIO

# load the sto
sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.dedup.sto"
alignment = AlignIO.read(sto_path, "stockholm")

# extract the Genbank ID from sequence identifier

# Regex pattern to match GenBank IDs
genbank_pattern = re.compile(r"(?:[A-Z]{2}\d{6}|NC_\d{6}|[A-Z]\d{5})")
# two letters + 6 digit num
# NC_ + 6 digit num

results = []

for record in alignment:
    seq_id = record.id
    match = genbank_pattern.search(seq_id)
    if match:
        genbank_id = match.group(0)
    else:
        genbank_id = None
    results.append([seq_id, genbank_id])


df = pd.DataFrame(results, columns=["seqID", "accession"])

has_none = df['accession'].isnull().any()

if has_none:
    print(df[df['accession'].isnull()])
else:
    print("All seqID contain accession number.")



# retrieve taxo info
# need manual check
# some associated viruses will be wrongly labeled as Polerovirus or not labeled

Entrez.email = "feiyue_yang@outlook.com"

def fetch_taxonomy_info(accession):
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession)
        record = Entrez.read(handle)
        if not record["IdList"]:
            return [], "Unknown"
        nucleotide_id = record["IdList"][0]
        
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="xml")
        record = Entrez.read(handle)
        taxonomy = record[0]["GBSeq_taxonomy"].split("; ")
        
        family = next((taxon for taxon in taxonomy if taxon.endswith("viridae")), taxonomy[-1])
        print(f"{accession} fetched")
        return taxonomy, family
    except Exception as e:
        print(f"Error fetching taxonomy for {accession}: {e}")
        return [], "Unknown"

# Add columns for taxonomy information
df['taxo_long'] = None
df['taxo'] = None

# Fetch and update taxonomy information for each accession number
for index, row in df.iterrows():
    accession = row['accession']
    if accession:
        taxonomy, family = fetch_taxonomy_info(accession)
        df.at[index, 'taxo_long'] = taxonomy
        df.at[index, 'taxo'] = family


df_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/taxonomy/merge_v0.align.mpkadded.edited.confident.dedup.taxo.tsv"

df.to_csv(df_path, index=False, sep="\t")
