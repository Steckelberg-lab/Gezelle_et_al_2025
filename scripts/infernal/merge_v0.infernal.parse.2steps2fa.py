# Parse the cmsearch --tblout output - 2steps.out
# intended for 2-step cmsearch
# generate a non-redundant low-similarity fasta file for alignment

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

out_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/infernal/v1_2/merge_v0.cmsearch.2steps.out"
fa_path = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/data/fasta/merge_v0.fasta"
re_fa_path = f"{out_path}.low_sim.fa"

df = pd.read_csv(out_path, sep="\t")
# dedup (just to make sure)
df = df.drop_duplicates(subset="sequence", keep="first")

fasta_sequences = list(SeqIO.parse(fa_path, 'fasta'))

# get the min number of unaligned bases for each row
def get_min_unaligned_bases(seq, fasta_sequences):
    print("aligning......")
    min_unaligned_bases = 0
    max_aligned_bases = 0
    for fasta_seq in fasta_sequences:
        # pairwise alignment
        # params from ST9_pairwise_align.py
        alignments = pairwise2.align.globalms(seq, str(fasta_seq.seq), 2, -0.5, -1.5, -1)
        best_alignment = max(alignments, key=lambda x: x[2])
        # Find the alignment with the maximum number of aligned bases
        seq1, seq2, score, begin, end = best_alignment
        aligned_bases = sum((1 for a, b in zip(seq1, seq2) if a == b and a != '-'))
        if aligned_bases > max_aligned_bases:
            max_aligned_bases = aligned_bases
    min_unaligned_bases = len(seq) - max_aligned_bases
    return min_unaligned_bases

df['min_unaligned_bases'] = df['sequence'].apply(get_min_unaligned_bases, fasta_sequences=fasta_sequences)
# TODO: align those seqs to each other and remove highly similar ones

print(len(df))
filtered_df = df[df['min_unaligned_bases'] >= 3]
print(len(filtered_df))

def make_fasta_record(row):
    header = f"{row['target_name']}.{row['seq_from']}-{row['seq_to']}.{row['strand']}.infernal"
    return SeqRecord(Seq(row['sequence']), id=header, description='')

seq_records = filtered_df.apply(make_fasta_record, axis=1).tolist()
SeqIO.write(seq_records, re_fa_path, "fasta")
