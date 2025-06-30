# Gezelle_et_al_2025
# xrRNA_search

**Gezelle et al. 2025. README for covariation analysis and search of new class 3 xrRNAs.**  
Based on the scripts and pipeline from Feiyue Yang and utilizing Infernal (E. P. Nawrocki and S. R. Eddy, Infernal 1.1: 100-fold faster RNA homology searches, *Bioinformatics* 29:2933–2935 (2013).)

This repository contains scripts for identifying and analyzing viral XRN1-resistant RNAs (xrRNAs), including alignment, structure-based search with Infernal, taxonomy classification, and metagene context analysis.

## Input Sequences

Start with:
- ST9a minimal sequence
- ST9a highly similar sequences (from BLASTn)
- Previously reported plant xrRNAs (see Steckelberg et al., RNA 2020: https://doi.org/10.1261/rna.076224.120)

## Step 1: Generate the initial STO File (Structured Alignment)

1. Align pseudoknots using anchor sequences:
   ```bash
   ./scripts/align/merge_v0_manual_align.py
   ```

2. Manually edit the `.tsv` file to align unanchored regions. Each row = sequence; each column = structural feature.

3. Convert the `.tsv` to a Stockholm `.sto` file:
   ```bash
   ./scripts/align/merge_v0_manual_align.sto.py
   ```

4. Verify full sequence coverage:
   ```bash
   ./scripts/align/merge_v0_manual_align.full_seq_check.py
   ```

## Step 2: Search with Infernal

1. Download the database from NCBI virus / BV-BRC

2. Build CM model with initial alighment created and calibrate (`cmbuild`, `cmcalibrate`), run `cmsearch` with the parameter -T 0

3. Parse results and extract hits using anchor-based pk search:
   ```bash
   ./scripts/infernal/merge_v0.infernal.pk_search.py
   ```

4. Manually review `.tsv` output, filter sequences, and convert to fasta format for alignment:
   ```bash
   ./scripts/infernal/merge_v0.infernal.pk_search.2fa.py
   ```

5. Align filtered sequences with `cmalign`, using the initial CM as reference

6. Manually check alignment and adjust sequences as needed.

7. Finalize sequence IDs:
   ```bash
   ./scripts/metagene/id_curation.py
   ```

## Step 3: Taxonomy Annotation

1. Retrieve taxonomy info from NCBI:
   ```bash
   ./scripts/taxonomy/retrieve_taxo.py
   ```

## Step 4: Metagene Context

1. Retrieve NCBI metadata:
   ```bash
   ./scripts/metagene/accession_info.py
   ```

2. Pull gene context annotations:
   ```bash
   ./scripts/metagene/metagene_all.py
   ```
   ```bash
   ./scripts/metagene/metagene_complete_genome.py
   ```

## Notes

- All paths in this guide assume relative locations.
- Manual review steps are noted where required.
- For questions, open an issue or contact the Steckelberg Lab.

