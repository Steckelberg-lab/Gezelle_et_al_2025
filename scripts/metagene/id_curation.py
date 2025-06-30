# Due to manual curation of sequences, 
# the seq coordinates in seq identifiers may not align well with the actual sequences in the sto file
import pandas as pd
import re
from Bio import Entrez, SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

input_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.dedup.sto"
output_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.dedup.IDcurated.sto"

alignment = AlignIO.read(input_path, "stockholm")

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
    sequence = str(record.seq).replace("-", "")
    sequence = sequence.upper()
    sequence = sequence.replace("U", "T")
    results.append([seq_id, genbank_id, sequence, str(record.seq)])


df = pd.DataFrame(results, columns=["seqID", "accession", "sequence", "sto_seq"])

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

# bug with - strand. manually curate them for now.
def extract_strand(seqID):
    parts = seqID.split('.')
    for part in parts:
        if part == '+' or part == '-':
            return part
    return "+"	# no strand info: seqs from RNA,2020

# Apply the functions to extract information and create new columns
df['seq_from'], df['seq_to'] = zip(*df['seqID'].apply(extract_seq_from_to))
df['strand'] = df['seqID'].apply(extract_strand)


# check the seq

Entrez.email = "fy2306@columbia.edu"

def fetch_sequence(accession, seq_from, seq_to):

    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", seq_start=seq_from, seq_stop=seq_to)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq).upper()
    except Exception as e:
        print(f"Error fetching sequence for {accession}: {e}")
        return None
    
def fetch_nucleotide_sequence(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq).upper()

def get_sequence_coordinates(accession, strand, original_sequence):
    nucleotide_sequence = fetch_nucleotide_sequence(accession)
    if strand == '-':
        original_sequence = str(Seq(original_sequence).reverse_complement())
    start_pos = nucleotide_sequence.find(original_sequence)
    # if start_pos == -1:
    #     raise ValueError(f"Sequence {original_sequence} not found in {accession} on {strand} strand")
    # Convert to 1-based coordinates
    seq_start = start_pos + 1
    seq_end = seq_start + len(original_sequence) - 1
    
    return seq_start, seq_end


def check_from_to(row):
    accession = row["accession"]
    seq_from = row["seq_from"]
    seq_to = row["seq_to"]
    strand = row['strand']
    original_sequence = row["sequence"]
    
    correct_sequence = fetch_sequence(accession, seq_from, seq_to)
    
    if correct_sequence and original_sequence != correct_sequence:
        print(f"Mismatch in seqID {row['seqID']}:")
        new_seq_from, new_seq_to = get_sequence_coordinates(accession, strand, original_sequence)
        if new_seq_from == 0:
            print(f"Correction failed: {accession}")
        print("------------------------------------------------")
    elif correct_sequence and original_sequence == correct_sequence:
        new_seq_from, new_seq_to = row["seq_from"], row["seq_to"]
        print(f"Correct match! {row['seqID']}")
        
    return new_seq_from, new_seq_to

df[['seq_from', 'seq_to']] = df.apply(lambda row: check_from_to(row), axis=1, result_type='expand')


# save        
def create_seq_id(row):
    return f"{row['accession']}.{row['strand']}.{row['seq_from']}-{row['seq_to']}"

# Create SeqRecord objects
seq_records = []
for index, row in df.iterrows():
    seq_id = create_seq_id(row)
    seq = Seq(row['sto_seq'])
    seq_record = SeqRecord(seq, id=seq_id, description="")
    seq_records.append(seq_record)

alignment = MultipleSeqAlignment(seq_records)

with open(output_path, "w") as output_handle:
    AlignIO.write(alignment, output_handle, "stockholm")