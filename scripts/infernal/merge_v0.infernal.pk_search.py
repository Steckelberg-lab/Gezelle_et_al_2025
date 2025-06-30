# Following intial step: cmsearch - SS_cons stem-loop
# Search for pk using anchor seqs

import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import gzip
import re

# Config
shared_path = "/path/to/shared_path"

re_path = f"{shared_path}/merge_v0.cmsearch.out"
interm_re_path = f"{shared_path}/manual_pk/merge_v0.cmsearch.manual_pk.interm.tsv"
final_re_path = f"{shared_path}/manual_pk/merge_v0.cmsearch.manual_pk.tsv"

db_path = "/path/to/ncbi_virus_Viridiplantae.fa.gz"
db_meta_path = "/path/to/ncbi_virus_Viridiplantae.metadata.csv"
sto_path = "/path/to/merge_v0.align.sto"    # initial sto

# ----------------------------------------------------------------------------------------------------------------------

# Generate the extended seq2table from step1 cmsearch output

# 1. parse the cmsearch output

re_df_cols = ["target_name", "accession", "query_name", "accession1", "mdl", "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc", "pass", "gc", "bias", "score", "E_value", "inc", "description_of_target"]

# Initialize a list to hold all the processed rows
data = []

# Open the file and read line by line
with open(re_path, 'r') as file:
    for line in file:
        # Ignore lines that are comments
        if line.startswith('#'):
            continue
        
        # Split the line at the first occurrence of '|'
        parts = line.split('|', 1)
        
        if len(parts) == 2:
            # The part before '|' is further split by whitespace
            first_parts = parts[0].strip().split()
            # The part after '|' is taken as is, with leading and trailing whitespaces stripped
            last_part = parts[1].strip()
            # Combine the split parts and the last part
            row = first_parts + [last_part]
            data.append(row)
        else:
            # Handle cases where there is no '|' character
            data.append(parts[0].strip().split())

# Create DataFrame
re_df = pd.DataFrame(data, columns=re_df_cols)


# 2. retrieve the seqs according to target_name
def extract_sequences(fasta_path, df):
    # Load the compressed FASTA file
    with gzip.open(fasta_path, "rt") as handle:
        # Create a dictionary from the FASTA file
         seq_dict = {record.id.split(' |')[0]: record for record in SeqIO.parse(handle, "fasta")}

    # Function to fetch sequence based on row in DataFrame
    def fetch_sequence(row):
        # Extract relevant sequence information
        seq_id = row['target_name']
        start = min(int(row['seq_from'])-1, int(row['seq_to'])-1)
        end = max(int(row['seq_from']), int(row['seq_to']))
        strand = row['strand']

        # Retrieve sequence and apply range and strand
        # Important: extend 40nt after the seq to include pk region
        # 40nt: chosen due to Lena's 2020 plrv paper consensus structure
        try:
            if strand == '+':
                seq = seq_dict[seq_id].seq[start:(end + 40)]
            if strand == '-':
                seq = seq_dict[seq_id].seq[(start - 40):end]
                seq = seq.reverse_complement()
            seq = str(seq).upper()
            seq = seq.replace('T', 'U')
            return seq
        except KeyError:
            return "Sequence not found"  

    # Apply function across DataFrame rows
    df['sequence'] = df.apply(fetch_sequence, axis=1)
    return df

re_df = extract_sequences(db_path, re_df)


# 3. retrieve the metadata according to target_name
meta_df = pd.read_csv(db_meta_path, sep=",", usecols=['Accession', 'Organism_Name', 'Species', 'Molecule_type'])
meta_df.columns = ["target_name", "Organism_Name", "Species", "Molecule_type"]
print(len(re_df))
re_df = pd.merge(re_df, meta_df, on='target_name', how='inner')
print(len(re_df)) # though we use inner-join, it should stay the same

# 4. Deduplication based on seqs
# remove the exactly same ones
re_df = re_df.drop_duplicates(subset='sequence', keep='first')
print(f"Initial deduplication based on seqs: {len(re_df)}")


# --------------------------------------------------------------------------------------------------------------

# Search for pk and trim the seq to 2nt after the pk

def iupac_to_regex(sequence):
    # Maps IUPAC ambiguous base symbols to regex character sets
    iupac_code = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'U',
        'R': '[AGR]', 'Y': '[CTUY]', 'S': '[GCS]', 'W': '[ATUW]',
        'K': '[GTUK]', 'M': '[ACM]', 'B': '[CGTUB]', 'D': '[AGTUD]',
        'H': '[ACTUH]', 'V': '[ACGV]', 'N': '[ACGTUN]'
    }
    regex_pattern = ''.join(iupac_code.get(base, base) for base in sequence)
    return regex_pattern

def find_subseq(seq, pattern):
    regex = iupac_to_regex(pattern)
    for match in re.finditer(regex, seq):
        yield match.start(), match.end()
        
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
    # return True if can pair, else False
    return any(b1 in pairings and b2 in expand_base(pairings[b1]) for b1 in base1_expanded for b2 in base2_expanded)


# pk should occur at the end of seq1

def find_pk(up_seq, dn_seq, name):
	# take the 4nt seq at the end of seq1
	initial_pk = up_seq[-4:]
	# reverse complement
	initial_pk_complement = Seq(initial_pk).reverse_complement()
	initial_pk_complement = initial_pk_complement.replace('T', 'U')
	# search in seq after UUGSAA
	initial_pk_found = list(find_subseq(dn_seq, initial_pk_complement))

	if initial_pk_found:
		pass
	else:
		# wobble base pairs
		# U/T - A/G (R)
		initial_pk_complement = initial_pk_complement.replace('A', 'R')
		# G - C/U (Y)
		initial_pk_complement = initial_pk_complement.replace('C', 'Y')
		initial_pk_found = list(find_subseq(dn_seq, initial_pk_complement))

	if initial_pk_found:
		initial_pk_found = initial_pk_found[0]
		pk_length = 4
		extend_pos = -5
		# if found, try to extend the pk
		while True:
			try:
				if can_pair(up_seq[extend_pos], dn_seq[initial_pk_found[0]+pk_length]):
					pk_length = pk_length + 1
					extend_pos = extend_pos - 1
				else:
					break
			except:
				print(up_seq, dn_seq, extend_pos, name)
				break 
		pk_2nt_pos = initial_pk_found[0] + pk_length + 2	# in dn_seq, 1-based
	else:
		pk_length = 0
		pk_2nt_pos = len(dn_seq)
	return pk_length, pk_2nt_pos
        

# For each row, find pk
def find_pk_row(row):
    
    # Search for the loop UUGSAA
    seq = row['sequence']
    L2_2_pos = list(find_subseq(seq, "UUGSRA"))
    
    # anchor seq found: return pk length (no pk : return 0)
    if L2_2_pos:
        L2_2_pos = L2_2_pos[0]
        L2_2_seq = seq[L2_2_pos[0]:L2_2_pos[1]]
        L2_2_dn_seq = seq[L2_2_pos[1]:]
        L2_2_up_seq = seq[:L2_2_pos[0]]

        pk_re, _ = find_pk(L2_2_up_seq, L2_2_dn_seq, row['target_name'])

	# no anchor seq found: return -1
    else:
        pk_re = -1

    return pk_re

re_df['pk_re'] = re_df.apply(find_pk_row, axis=1)

# trim the seq to 2nt after the pk
def trim_row(row):
    # Search for the loop UUGSAA
    seq = row['sequence']
    L2_2_pos = list(find_subseq(seq, "UUGSRA"))
    
    # anchor seq found: return pk length (no pk : return 0)
    if L2_2_pos:
        L2_2_pos = L2_2_pos[0]
        L2_2_seq = seq[L2_2_pos[0]:L2_2_pos[1]]
        L2_2_dn_seq = seq[L2_2_pos[1]:]
        L2_2_up_seq = seq[:L2_2_pos[0]]

        _, trim_dn_seq_pos = find_pk(L2_2_up_seq, L2_2_dn_seq, row['target_name'])

        trim_seq = L2_2_up_seq + L2_2_seq + L2_2_dn_seq[:trim_dn_seq_pos]

	# no anchor seq found: return -1
    else:
        trim_seq = seq

    return trim_seq

re_df['sequence'] = re_df.apply(trim_row, axis=1)

# Dedup
print("Final dedup")
re_df = re_df.drop_duplicates(subset='sequence', keep='first')
print(f"Deduplication based on seqs: {len(re_df)}")

sto = AlignIO.read(sto_path, 'stockholm')
sto_data = [{'seq_name': record.id, 'sequence': str(record.seq).replace('-', '')} for record in sto]
sto_df = pd.DataFrame(sto_data)
re_df = re_df[~re_df["sequence"].isin(sto_df["sequence"].tolist())]
print(f"Deduplication based on sto input seqs: {len(re_df)}")

# Save results
re_df.to_csv(interm_re_path, sep="\t", index=False)

re_df = re_df[re_df["pk_re"] > 0]
re_df.to_csv(final_re_path, sep="\t", index=False)