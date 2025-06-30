import re
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# load the sto
sto_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/manual_pk_cmalign/merge_v0.align.mpkadded.edited.confident.dedup.sto"
fasta_output_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/NoFold_cm_all/merge_v0.align.mpkadded.edited.confident.dedup.NoFoldID.fasta"

# Read the Stockholm file
alignment = AlignIO.read(sto_path, "stockholm")

# Regex pattern to match GenBank IDs
genbank_pattern = re.compile(r"(?:[A-Z]{2}\d{6}|NC_\d{6}|[A-Z]\d{5})")

results = []

# Extract GenBank ID and process sequences
for record in alignment:
    seq_id = record.id
    match = genbank_pattern.search(seq_id)
    if match:
        genbank_id = match.group(0)
    else:
        genbank_id = None
    # Remove gaps from the sequence
    sequence = str(record.seq).replace("-", "")
    genbank_id = genbank_id.replace("_","-")
    results.append(SeqRecord(Seq(sequence), id=genbank_id if genbank_id else seq_id, description=""))

# Save to FASTA file
SeqIO.write(results, fasta_output_path, "fasta")