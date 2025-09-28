# Serotype and genomic diversity of dengue virus during the 2023 outbreak in Pakistan reveals the circulation of genotype III of DENV-1 and cosmopolitan genotype of DENV-2

Reproducible code and workflow that mirror the analysis in the manuscript from Journal of Medical Virology (2024). DOI: 10.1002/jmv.29727

## Program summary
End to end analysis of DENV-1 and DENV-2 using metagenomic NGS. Steps match the study design and are fully scripted to allow reviewers and collaborators to reproduce the results.

1) Inputs
   - Paired-end FASTQ from serum RNA libraries sequenced on Illumina MiSeq 2x150.

2) Quality control and trimming
   - FastQC for initial QC.
   - Trimmomatic for adapter and quality trimming with the parameters reported in the paper.
   - Picard MarkDuplicates to remove PCR duplicates.

3) De novo assembly and contig validation
   - Assemble with SPAdes.
   - Contig QC: length, percent Ns, GC, N50.
   - BLAST remote against NCBI nt to identify closest matches.
   - Download selected GenBank sequences for context trees and for mapping references.

4) Reference-based mapping and masked consensus
   - For each sample, choose the best reference from BLAST hits and map with BWA-MEM.
   - Sort, mark duplicates, compute per-base depth.
   - Mask sites below a coverage threshold to N and call variants with bcftools.
   - Build a masked consensus per sample. The default minimum depth is 10, matching the study.

5) Phylogeny
   - MAFFT alignment of sample consensus plus context genomes.
   - IQ-TREE with 1000 ultrafast bootstraps.
   - Default substitution model is GTR+G+I to match the manuscript. You can switch to ModelFinder with `-m MFP` in config if desired.
   - Trees saved in Newick format with logs for review.

6) Mutation analysis by protein
   - Trim to the single polyprotein ORF using the annotated reference and translate.
   - Compare amino acids against reference DENV-1 polyprotein AAW64436.1 and DENV-2 reference NC_001474.1.
   - Report differences per protein region including E, M, C, NS1 to NS5. Output is a tidy TSV.

7) Outputs
   - Per sample consensus: `results/consensus/<sample>.fa`
   - Combined consensus: `results/consensus/all_consensus.fasta`
   - Alignments: `results/aln/wg_alignment.fasta`
   - Phylogeny: `results/iqtree/wg.treefile` and `.log`
   - Mutation tables: `results/mutations/<sample>_aa_diffs.tsv`
   - QC, BLAST tables, and reference selection under `results/`

8) Repro and compliance
   - MIT LICENSE, CITATION.cff.
   - CI sanity check for plotting.
   - Never commit restricted or clinical data. Set NCBI_EMAIL once for Entrez-based steps.

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake for the full pipeline
- BLAST+ with remote enabled, Entrez Direct

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One-command end to end run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate denv-2023-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Example:
```yaml
pairs:
  - sample: DENV1_001
    r1: data-private/DENV1_001_R1.fastq.gz
    r2: data-private/DENV1_001_R2.fastq.gz
  - sample: DENV2_001
    r1: data-private/DENV2_001_R1.fastq.gz
    r2: data-private/DENV2_001_R2.fastq.gz

reference_preference:
  denv1_ref_accessions:
    - OR936753.1
    - OR936754.1
    - OR936755.1
    - OR936756.1
    - OR936757.1
    - OR936758.1
    - OR936759.1
    - OR936760.1
    - OR936761.1
    - OR936762.1
    - OR936763.1
    - OR936764.1
    - OR936765.1
    - OR936766.1
    - OR936767.1
    - OR936768.1
    - OR936769.1
    - OQ445883.1
    - ON908217.1
    - ON873938.1
    - OP898557.1
    - OP898558.1
    - NC_001477.1

  denv2_ref_accessions:
    - OR936752.1
    - OP811977.1
    - OP811984.1
    - MK543472.1
    - MW186240.1
    - NC_001474.1

context_accessions:
  - OR936753.1
  - OR936754.1
  - OR936755.1
  - OR936756.1
  - OR936757.1
  - OR936758.1
  - OR936759.1
  - OR936760.1
  - OR936761.1
  - OR936762.1
  - OR936763.1
  - OR936764.1
  - OR936765.1
  - OR936766.1
  - OR936767.1
  - OR936768.1
  - OR936769.1
  - OR936752.1
  - ON908217.1
  - ON873938.1
  - OP898557.1
  - OP898558.1
  - MK543472.1
  - MW186240.1

params:
  threads: 4
  trim_adapters: env/TruSeq3-PE.fa
  min_depth_consensus: 10
  min_qual: 20
  iqtree_model_wg: GTR+G+I
  bootstrap: 1000
  max_blast_hits: 50
```
