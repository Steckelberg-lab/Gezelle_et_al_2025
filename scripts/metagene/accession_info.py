import pandas as pd
import re
from Bio import Entrez, SeqIO, AlignIO


# load the sto
sto_path = "/path/to/merge_v0.align.mpkadded.edited.confident.dedup.IDcurated.sto"
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


# retrieve seq_from, seq_to, and strand
seq_pattern = re.compile(r'.*?(\d+)-(\d+).*')

# Functions to extract seq_from, seq_to, and strand
def extract_seq_from_to(seqID):
    match = seq_pattern.match(seqID)
    if match:
        seq1 = int(match.group(1))
        seq2 = int(match.group(2))
        return min(seq1, seq2), max(seq1, seq2)
    return None, None

def extract_strand(seqID):
    parts = seqID.split('.')
    for part in parts:
        if part == '+' or part == '-':
            return part
    return "+"	# no strand info: seqs from RNA,2020

# Apply the functions to extract information and create new columns
df['seq_from'], df['seq_to'] = zip(*df['seqID'].apply(extract_seq_from_to))
df['strand'] = df['seqID'].apply(extract_strand)

# retrieve description
Entrez.email = "name@email.com"

# Function to fetch description information
def fetch_description(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        description = record.description
        print(f"{accession} fetched.")
        return description
    except Exception as e:
        print(f"Error fetching description for {accession}: {e}")
        return None

# Fetch and update description information
df['description_of_target'] = df['accession'].apply(fetch_description)


df['complete'] = df['description_of_target'].apply(
    lambda x: 'partial' if isinstance(x, str) and 'partial' in x.lower() else
              'complete' if isinstance(x, str) and 'complete' in x.lower() else
              'unknown' if isinstance(x, str) else None
)
df['genome'] = df['description_of_target'].apply(
    lambda x: 'cds' if isinstance(x, str) and 'cds' in x.lower() else
              'genome' if isinstance(x, str) and 'genome' in x.lower() else
              'unknown' if isinstance(x, str) else None
)


df_save_path = "/path/to/merge_v0.align.mpkadded.edited.confident.dedup.info.tsv"
df.to_csv(df_save_path, sep="\t", index=False)
