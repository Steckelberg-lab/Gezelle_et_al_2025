# Tried to use pairwise alignment to accelerate manual alignment
# However, it's not very accurate
from Bio import SeqIO, pairwise2
import random

filename = "/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/data/fasta/merge_v0.fasta"
sequences = list(SeqIO.parse(filename, "fasta"))

def calculate_details(alignment):
    seq1, seq2, score, begin, end = alignment
    matches = sum((1 for a, b in zip(seq1, seq2) if a == b))
    length = max(len(seq1), len(seq2))
    identity_score = matches / length
    return matches, identity_score, score

def shuffle_str(s):
    chars = list(s)
    random.shuffle(chars)
    return ''.join(chars)

# File to save the summary of the results
with open('/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/ST9_merge_v0.tsv', 'w') as summary_file, open('/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/ST9_merge_v0.details.txt', 'w') as details_file:
   # Write headers for the TSV file
    summary_file.write("seq1\tseq2\tseq1_length\tseq2_length\tmax_length\tmatches\tidentity\n")
    
    seq1 = sequences[0]  # ST9
    for seq_record in sequences[:]:  
        # Perform alignment and sort by the highest score
        # identical
        # mismatch
        # gap
        # extension
        alignments = pairwise2.align.globalms(seq1.seq, seq_record.seq, 2, -0.5, -1.5, -1)
        if alignments:
            # Select the alignment with the highest score
            best_alignment = max(alignments, key=lambda x: x[2])
            matches, identity_score, best_score = calculate_details(best_alignment)
            
            # Write to summary TSV file
            summary_file.write(f"{seq1.id}\t{seq_record.id}\t{len(seq1.seq)}\t{len(seq_record.seq)}\t{max(len(seq1.seq), len(seq_record.seq))}\t{matches}\t{identity_score:.4f}\n")
            
            # Write the best detailed alignment to a separate file
            details_file.write(f"Best alignment between {seq1.id} and {seq_record.id} with matches {matches} and identity {identity_score:.4f}:\n")
            details_file.write(f"{best_alignment[0]}\n{best_alignment[1]}\n\n")

print("Alignment data has been saved.")

# randomized input
with open('/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/ST9_merge_v0.random.tsv', 'w') as summary_file, open('/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/align/ST9_merge_v0.random.details.txt', 'w') as details_file:
   # Write headers for the TSV file
    summary_file.write("seq1\tseq2\tseq1_length\tseq2_length\tmax_length\tmatches\tidentity\n")
    
    seq1 = sequences[0]  # ST9
    for seq_record in sequences[:]:  
        # Perform alignment and sort by the highest score
        # identical
        # mismatch
        # gap
        # extension
        # shuffle the seqs
        seq_record_seq = shuffle_str(seq_record.seq)
        alignments = pairwise2.align.globalms(seq1.seq, seq_record_seq, 2, -0.5, -1.5, -1)
        if alignments:
            # Select the alignment with the highest score
            best_alignment = max(alignments, key=lambda x: x[2])
            matches, identity_score, best_score = calculate_details(best_alignment)
            
            # Write to summary TSV file
            summary_file.write(f"{seq1.id}\t{seq_record.id}\t{len(seq1.seq)}\t{len(seq_record_seq)}\t{max(len(seq1.seq), len(seq_record_seq))}\t{matches}\t{identity_score:.4f}\n")
            
            # Write the best detailed alignment to a separate file
            details_file.write(f"Best alignment between {seq1.id} and {seq_record.id} with matches {matches} and identity {identity_score:.4f}:\n")
            details_file.write(f"{best_alignment[0]}\n{best_alignment[1]}\n\n")

print("Alignment data has been saved. (Randomized)")