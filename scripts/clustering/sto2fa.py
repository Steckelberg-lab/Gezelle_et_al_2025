from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sto_path = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/merge_v0.align.dedup.clustering.sto'
alignment = AlignIO.read(sto_path, 'stockholm')

records = []
for record in alignment:
    # Remove gaps from the sequence
    sequence = str(record.seq).replace("-", "")
    # Create a new SeqRecord
    new_record = SeqRecord(Seq(sequence), id=record.id, description="")
    records.append(new_record)
    
fasta_path = '/ifs/scratch/as6282_gp/fy2306/projects/xrRNA_search/results/clustering/merge_v0.align.dedup.clustering.fasta'
SeqIO.write(records, fasta_path, 'fasta')