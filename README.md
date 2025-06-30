# Gezelle_et_al_2025
# xrRNA_search

**Gezelle et al. 2025. README for covariation analysis and search of new class 3 xrRNAs.**  
Based on the scripts and pipeline from Feiyue Yang and utilizing Infernal (E. P. Nawrocki and S. R. Eddy, Infernal 1.1: 100-fold faster RNA homology searches, *Bioinformatics* 29:2933–2935 (2013).)

This repository contains scripts for identifying and analyzing viral XRN1-resistant RNAs (xrRNAs), including alignment, structure-based search with Infernal, clustering, taxonomy classification, and metagene context analysis.

## Input Sequences

Start with:
- ST9a minimal sequence
- ST9a highly similar sequences (from BLASTn)
- Previously reported plant xrRNAs (see Steckelberg et al., RNA 2020: https://doi.org/10.1261/rna.076224.120)

## Step 1: Generate the STO File (Structured Alignment)

Scripts: `scripts/align/`

1. Align pseudoknots using anchor sequences:
   ```bash
   python merge_v0_manual_align.py input.fasta > output.tsv
   ```

2. Manually edit the `.tsv` file to align unanchored regions. Each row = sequence; each column = structural feature.

3. Convert the `.tsv` to a Stockholm `.sto` file:
   ```bash
   python merge_v0_manual_align.sto.py aligned.tsv > aligned.sto
   ```

4. Verify full sequence coverage:
   ```bash
   python merge_v0_manual_align.full_seq_check.py aligned.sto
   ```

## Step 2: Search with Infernal

Scripts: `scripts/infernal/`

1. Download viral genomes from NCBI Virus or use a local database (e.g., `data/ncbi_virus/`).

2. Run cmsearch with Infernal:
   ```bash
   ./merge_v0.infernal.sh
   ```

3. Parse results and extract hits using anchor-based pk search:
   ```bash
   python merge_v0.infernal.pk_search.py cmsearch_output.tbl
   ```

4. Manually review `.tsv` output and filter sequences:
   ```bash
   python merge_v0.infernal.pk_search.2fa.py filtered.tsv > output.fasta
   ```

5. Align filtered sequences with `cmalign`:
   ```bash
   ./merge_v0.infernal.cmalign.sh
   ```

6. Manually check alignment and adjust sequences as needed.

7. Finalize sequence IDs:
   ```bash
   python ../metagene/id_curation.py corrected_alignment.sto
   ```

## Step 3: Taxonomy Annotation

Scripts: `scripts/taxonomy/`

1. Retrieve taxonomy info from NCBI:
   ```bash
   python retrieve_taxo.py sequence_list.txt
   ```

2. Split `.sto` file into subsets for R-scape:
   ```bash
   python split_sto.py alignment.sto taxonomy.tsv
   ```

## Step 4: Structural Clustering (NoFold)

Scripts: `scripts/clustering/`  
NoFold GitHub: https://github.com/kimpenn/nofold

1. Convert `.sto` to `.fasta` for NoFold:
   ```bash
   python merge_v0.mpkadded.confident.nofold.id.py input.sto > input.fasta
   ```

2. Run NoFold scoring:
   ```bash
   ./merge_v0.mpkadded.confident.nofold.score_cm.sh
   ```

3. Create taxonomy labels for PCA plotting:
   ```bash
   python merge_v0.mpkadded.confident.nofold.label.py
   ```

4. Generate PCA plots:
   ```bash
   python merge_v0_test.nofold.pca_plot.py
   ```

## Step 5: Metagene Context

Scripts: `scripts/metagene/`

1. Retrieve NCBI metadata:
   ```bash
   python accession_info.py
   ```

2. Pull gene context annotations:
   ```bash
   python metagene_all.py
   python metagene_complete_genome.py
   python metagene_complete_cds.py
   ```

## Directory Structure

```
scripts/
├── align/
├── clustering/
├── infernal/
├── metagene/
├── taxonomy/
data/
results/
```

## Requirements

- Python 3+
- Infernal: http://eddylab.org/infernal/
- NoFold: https://github.com/kimpenn/nofold
- Biopython, NumPy, Pandas, Matplotlib

## Notes

- All paths in this guide assume relative locations.
- Manual review steps are noted where required.
- For questions, open an issue or contact the Steckelberg Lab.

