import re
from Bio import SeqIO,pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import numpy as np

filename = "/path/to/merge_v0.fasta"
sequences = list(SeqIO.parse(filename, "fasta"))

def iupac_to_regex(sequence):
    # Maps IUPAC ambiguous base symbols to regex character sets
    iupac_code = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U',
        'R': '[AG]', 'Y': '[CTU]', 'S': '[GC]', 'W': '[ATU]',
        'K': '[GTU]', 'M': '[AC]', 'B': '[CGTU]', 'D': '[AGTU]',
        'H': '[ACTU]', 'V': '[ACG]', 'N': '[ACGTU]'
    }
    regex_pattern = ''.join(iupac_code.get(base, base) for base in sequence)
    return regex_pattern

def find_subseq(seq, pattern):
    regex = iupac_to_regex(pattern)
    for match in re.finditer(regex, seq):
        yield match.start(), match.end()

# This is very stringent
# pk should occur at the end of seq1 and at the start of seq2

def can_pair(base1, base2):
    # Expand the base pairing rules to include ambiguous bases
    # include wobble
    pairings = {
        'A': 'UT', 'U': 'AG', 'T': 'AG', 'C': 'G', 'G': 'CUT',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'
    }
    # Convert bases to sets of possible matches
    def expand_base(base):
        iupac_expansions = {
            'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U',
            'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
            'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
        }
        return iupac_expansions.get(base, base)

    base1_expanded = expand_base(base1)
    base2_expanded = expand_base(base2)
    return any(b1 in pairings and b2 in expand_base(pairings[b1]) for b1 in base1_expanded for b2 in base2_expanded)

def find_pk(seq1, seq2):
    min_length = min(len(seq1), len(seq2))
    reversed_seq1 = seq1[-min_length:][::-1]
    pairs = []

    for i in range(min_length):
        if can_pair(reversed_seq1[i], seq2[i]):
            pairs.append((len(seq1) - 1 - i, i))
        else:
            break
            
    # Now pairs: [(seq1_index,seq2_index)]
    # be stringent on pair length to avoid occurence of pairs by chance
    return pairs[-1] if len(pairs)>3 else ()


# initialize results df
columns = ['seqname', 'L2_pk_up', 'L2_pk', 'L2_2', 'L2_2_dn', 'J1_3_3nt', 'J1_3_pk', 'end_3', 'full_seq']
df = pd.DataFrame(columns=columns)

for seq in sequences[:]:
    
	# Search for the loop UUGSAA
    L2_2_pos = list(find_subseq(str(seq.seq), "UUGSAA"))
    if L2_2_pos:
          L2_2_pos = L2_2_pos[0]
          # extract L2_2
          L2_2_seq = str(seq.seq)[L2_2_pos[0]:L2_2_pos[1]]
          
          L2_2_dn_seq = str(seq.seq)[L2_2_pos[1]:]
          L2_2_up_seq = str(seq.seq)[:L2_2_pos[0]]
          # Search for AGY & AAC (also searching for AAU using this regex)
          J1_3_3nt_pos = list(find_subseq(L2_2_dn_seq, "ARY"))
          if J1_3_3nt_pos:
                for i in range(0,len(J1_3_3nt_pos)):
                    J1_3_3nt_pos_i = J1_3_3nt_pos[i]
                    J1_3_3nt_seq = L2_2_dn_seq[J1_3_3nt_pos_i[0]:J1_3_3nt_pos_i[1]]
                    J1_3_3nt_dn_seq = L2_2_dn_seq[J1_3_3nt_pos_i[1]:]
                    # Search for PK
                    pk_pos = find_pk(L2_2_up_seq, J1_3_3nt_dn_seq)
                    # pk_pos = (seq1_start, seq2_end)
                    J1_3_3nt_pos_to_choose = J1_3_3nt_pos[i]
                    if pk_pos != ():
                        break
                    else:
                        continue
                J1_3_3nt_seq = L2_2_dn_seq[J1_3_3nt_pos_to_choose[0]:J1_3_3nt_pos_to_choose[1]]
                J1_3_3nt_dn_seq = L2_2_dn_seq[J1_3_3nt_pos_to_choose[1]:]
                pk_pos = find_pk(L2_2_up_seq, J1_3_3nt_dn_seq)
                if pk_pos != ():
                    L2_pk_seq = L2_2_up_seq[pk_pos[0]:]
                    J1_3_pk_seq = J1_3_3nt_dn_seq[:pk_pos[1]+1]
                    end_3_seq = J1_3_3nt_dn_seq[pk_pos[1]+1:]
                    
                    # Search for stems and bulding loops (need to optimize)
                    # using pairwise2.align from biopython: cannot generate ideal alignment results
                    
                    # record the L2-pk, loop UUGSAA, AGY, J1_3_pk, end_3 / full seq
                    new_row = pd.DataFrame(
						{
							'seqname': [str(seq.id)], 
							'L2_pk_up': [L2_2_up_seq[:pk_pos[0]]], 
							'L2_pk': [L2_pk_seq], 
							'L2_2': [L2_2_seq], 
							'L2_2_dn': [L2_2_dn_seq[:J1_3_3nt_pos_to_choose[0]]], 
							'J1_3_3nt': [J1_3_3nt_seq], 
							'J1_3_pk': [J1_3_pk_seq], 
							'end_3': [end_3_seq], 
							'full_seq': [str(seq.seq)]
						}
							)

                else:
                     # record the loop UUGSAA, AGY / full seq
                     new_row = pd.DataFrame(
						{
							'seqname': [str(seq.id)], 
							'L2_pk_up': [L2_2_up_seq], 
							'L2_pk': [None], 
							'L2_2': [L2_2_seq], 
							'L2_2_dn': [L2_2_dn_seq[:J1_3_3nt_pos_to_choose[0]]], 
							'J1_3_3nt': [J1_3_3nt_seq], 
							'J1_3_pk': [None], 
							'end_3': [J1_3_3nt_dn_seq], 
							'full_seq': [str(seq.seq)]
						}
							)
                
          else:
                # record the loop UUGSAA / full seq
                new_row = pd.DataFrame(
					{
						'seqname': [str(seq.id)], 
						'L2_pk_up': [L2_2_up_seq], 
						'L2_pk': [None], 
						'L2_2': [L2_2_seq], 
						'L2_2_dn': [None], 
						'J1_3_3nt': [None], 
						'J1_3_pk': [None], 
						'end_3': [L2_2_dn_seq], 
						'full_seq': [str(seq.seq)]
					}
						)
                
          

    else:
        # record the full seq
        new_row = pd.DataFrame(
            {
                'seqname': [str(seq.id)], 
                'L2_pk_up': [None], 
                'L2_pk': [None], 
                'L2_2': [None], 
                'L2_2_dn': [None], 
                'J1_3_3nt': [None], 
                'J1_3_pk': [None], 
                'end_3': [None], 
                'full_seq': [str(seq.seq)]
            }
                )
    
    df = pd.concat([df, new_row], ignore_index=True)
        

df.to_csv('/path/to/merge_v0.align.0.tsv', sep='\t', index=False)