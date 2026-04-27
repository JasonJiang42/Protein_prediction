# Protein Homology Analysis Pipeline

A unified computational pipeline for **homolog prediction** using sequence similarity, protein domain architecture, and 3D structure modeling.

---

## Overview

This pipeline identifies and characterizes homologous proteins through a multi-stage analysis:

1. **Clustering** — Pan-proteome clustering with MMseqs2 or CD-HIT
2. **Sequence** — DIAMOND blastp search with identity/coverage filtering
3. **Domain** — Pfam domain annotation via hmmscan
4. **Domain Filter** — Restrict hits to those overlapping query domain regions
5. **Structure** — ColabFold structure prediction + Foldseek alignment
6. **Visualization** — BLAST coverage plots, domain architecture, Foldseek scatter plots

---

## Quick Start

### Prerequisites

```bash
# Required tools
diamond           # DIAMOND blastp
hmmscan           # HMMER profile search
colabfold_batch   # ColabFold API MSA + structure prediction
foldseek          # Structural alignment
Rscript           # R environment for plotting

# Required databases
Pfam-A.hmm        # Pfam database (download from https://pfam.xfam.org/)
                  # Must be hmmpress-compressed

# Python libraries
pandas            # Data manipulation
matplotlib        # Plotting (for --plot)
numpy             # Numerical operations
```

### Installation

```bash
# Install Python dependencies
pip install pandas matplotlib numpy

# Download and compress Pfam database
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### Basic Usage

```bash
# Run all local steps (clustering → sequence → domain → structure prep)
python3 prot_homology.py run-all \
  -i proteins_by_species/all_species.faa \
  -q query.faa

# Or run individual modules
python3 prot_homology.py clustering -i proteins.faa
python3 prot_homology.py sequence -q query.faa --db protein_cluster/panproteome_db
python3 prot_homology.py domain --subject <clustered_faa> -q query.faa
```

---

## Modules

### 1. Clustering

**Purpose**: Build a pan-proteome database from subject sequences.

```bash
python3 prot_homology.py clustering \
  -i proteins_by_species/all_species.faa \
  --tool mmseqs \                    # mmseqs | cdhit
  -c 0.9 \                           # Sequence identity (0–1)
  --clust-qcov 0.8 \                 # Query coverage
  --clust-scov 0.8 \                 # Subject coverage
  -T 16                              # CPU threads
```

**Output**:
- `protein_cluster/mmseqs_c0.9_cov0.8.renamed.faa` — Clustered sequences
- `protein_cluster/panproteome_db.dmnd` — DIAMOND database

---

### 2. Sequence

**Purpose**: DIAMOND blastp search of query sequences against pan-proteome.

```bash
python3 prot_homology.py sequence \
  -q query.faa \
  --db protein_cluster/panproteome_db \
  --seq-output query_vs_pan.tsv \
  --identity 20.0 \                  # Min % identity
  --max-hits 0 \                     # Max targets per query (0 = unlimited)
  -T 16
```

**Output**:
- `query_vs_pan.tsv` — Full BLAST results
- `query_vs_pan_id20.tsv` — Identity-filtered results

---

### 3. Domain

**Purpose**: Pfam hmmscan annotation of query and subject sequences.

```bash
python3 prot_homology.py domain \
  --subject protein_cluster/mmseqs_c0.9_cov0.8.renamed.faa \
  -q query.faa \
  --pfam /path/to/Pfam-A.hmm \
  --dom-outdir domain/ \
  -T 16
```

**Output**:
- `domain/pfam_results.tsv` — Combined GA-filtered domain hits
- `domain/subject.domtblout`, `domain/query.domtblout` — Raw hmmscan output

---

### 4. Domain Filter

**Purpose**: Restrict BLAST hits to those overlapping a user-specified Pfam domain in queries.

```bash
python3 prot_homology.py domain-filter \
  --domains Peptidase_M64 \           # or accession: PF03009
  --pfam-tsv domain/pfam_results.tsv \
  --blast query_vs_pan_id20.tsv \
  --out domain/function_domain_filter.tsv
```

**Output**:
- `domain/function_domain_filter.tsv` — Filtered BLAST hits with matched domain info

---

### 5. Structure — MSA

**Purpose**: Generate multiple sequence alignments using ColabFold EBI API (run on login node).

```bash
python3 prot_homology.py structure-msa \
  --fasta structure_prediction/query_filtered.faa \
  --output structure_prediction/
```

**Output**:
- `structure_prediction/*/a3m/` — MSA files for each sequence

---

### 6. Structure — Predict

**Purpose**: Submit ColabFold structure prediction job to HPC cluster.

```bash
python3 prot_homology.py structure-predict \
  --input structure_prediction/ \
  --output /scratch/jiangshuo/colabfold_output/ \
  --weights /path/to/colabfold_weights \
  --partition NV_4090D \
  --gpus 1 \
  --cpus 8
```

**Workflow**:
1. Generates SLURM script: `submit_colabfold.sh`
2. Submits job automatically
3. Monitor: `squeue -j <job_id>`

---

### 7. Structure — Align

**Purpose**: Foldseek structural alignment of predicted PDB structures.

```bash
python3 prot_homology.py structure-align \
  --query-dir structure_prediction/output/ \
  --targets pdb_database/ \
  --output structure_prediction/foldseek_results.tsv \
  --threshold 0.5 \                  # TM-score cutoff
  -T 8
```

**Output**:
- `structure_prediction/foldseek_results.tsv` — Structural alignment results

---

### 8. Integrate

**Purpose**: Merge BLAST (sequence) + Foldseek (structure) results.

```bash
python3 prot_homology.py integrate \
  --blast query_vs_pan_id20.tsv \
  --foldseek structure_prediction/foldseek_results.tsv \
  --output results_figure/integrated_results.tsv
```

**Output**:
- `results_figure/integrated_results.tsv` — Combined BLAST + Foldseek + domain info

---

### 9. Plot

**Purpose**: Generate publication-quality figures.

```bash
# BLAST alignment coverage (one figure per query)
python3 prot_homology.py plot \
  --blast query_vs_pan_id20.tsv \
  --out results_figure/

# Domain architecture (requires --subject)
python3 prot_homology.py plot \
  --domain domain/ \
  --subject protein_cluster/mmseqs_c0.9_cov0.8.renamed.faa \
  --out results_figure/

# Foldseek scatter (sequence identity vs structural similarity)
python3 prot_homology.py plot \
  --foldseek results_figure/integrated_results.tsv \
  --blast-threshold 30.0 \            # BLAST identity line
  --align-threshold 0.5 \             # Foldseek TM-score for star marker
  --lddt-threshold 0.5 \              # pLDDT for star marker
  --out results_figure/

# All plots at once
python3 prot_homology.py plot --all --subject <panproteome.faa>
```

**Output**:
- `results_figure/blast_alignment_plot_*.{pdf,png}` — BLAST coverage bars
- `results_figure/domain_architecture.{pdf,png}` — gggenes domain plots
- `results_figure/foldseek_scatter.{pdf,png}` — Integrated scatter plot

---

## Workflow Example

### Local Steps (Single Machine)

```bash
# 1. Cluster pan-proteome
python3 prot_homology.py clustering \
  -i proteins_by_species/all_species.faa \
  --tool mmseqs -c 0.9

# 2. DIAMOND blastp
python3 prot_homology.py sequence \
  -q query.faa \
  --db protein_cluster/panproteome_db \
  --identity 20

# 3. Pfam annotation
python3 prot_homology.py domain \
  --subject protein_cluster/mmseqs_c0.9_cov0.8.renamed.faa \
  -q query.faa

# 4. Filter by domain overlap
python3 prot_homology.py domain-filter \
  --domains Peptidase_M64 \
  --pfam-tsv domain/pfam_results.tsv \
  --blast query_vs_pan_id20.tsv

# 5. Generate MSAs for ColabFold
python3 prot_homology.py structure-msa \
  --fasta domain/query_filtered.faa \
  --output structure_prediction/
```

### HPC Steps (GPU cluster)

```bash
# 6. Submit ColabFold job
python3 prot_homology.py structure-predict \
  --input structure_prediction/ \
  --output colabfold/ \
  --weights /path/to/weights \
  --partition NV_4090D --gpus 1 --cpus 8

# [Monitor on cluster, download results]

# 7. Run Foldseek alignment
python3 prot_homology.py structure-align \
  --query-dir structure_prediction/output/ \
  --targets pdb_database/ \
  --output structure_prediction/foldseek_results.tsv

# 8. Integrate results
python3 prot_homology.py integrate \
  --blast query_vs_pan_id20.tsv \
  --foldseek structure_prediction/foldseek_results.tsv \
  --output results_figure/integrated_results.tsv
```

### Visualization

```bash
# Generate all publication figures
python3 prot_homology.py plot --all \
  --subject protein_cluster/mmseqs_c0.9_cov0.8.renamed.faa
```

---

## Output Files

| File | Description |
|------|-------------|
| `protein_cluster/*.faa` | Clustered pan-proteome sequences |
| `protein_cluster/panproteome_db.dmnd` | DIAMOND database index |
| `query_vs_pan.tsv` | Full DIAMOND blastp results |
| `query_vs_pan_id*.tsv` | Identity-filtered BLAST hits |
| `domain/pfam_results.tsv` | Pfam domain annotations (GA-filtered) |
| `domain/function_domain_filter.tsv` | BLAST hits filtered by domain overlap |
| `structure_prediction/*.pdb` | ColabFold predicted structures |
| `structure_prediction/foldseek_results.tsv` | Foldseek structural alignments |
| `results_figure/integrated_results.tsv` | **Final results table** (all data merged) |
| `results_figure/*.pdf/.png` | Publication-ready figures |

---

## Key Parameters

### Clustering
- `-c, --clust-identity` — Sequence identity (0–1, default: 0.9)
- `--clust-qcov, --clust-scov` — Query/subject coverage (default: 0.8)

### Sequence Search
- `--identity` — Min % identity for filtered output (default: 20%)
- `--seq-qcov, --seq-scov` — Coverage thresholds passed to DIAMOND
- `--max-hits` — Max targets per query (default: 0 = unlimited)

### Domain Filtering
- `--domains` — Comma-separated Pfam names or accessions (e.g., `Peptidase_M64,PF03009`)

### Foldseek Visualization
- `--blast-threshold` — BLAST identity threshold line (default: 30%)
- `--align-threshold` — Foldseek TM-score for highlighting (default: 0.5)
- `--lddt-threshold` — pLDDT score for highlighting (default: 0.5)

---

## Troubleshooting

### Issue: "Pfam-A.hmm not found"
```bash
# Ensure Pfam database is downloaded and compressed
cd /path/to/pfamdb
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### Issue: "No hits found after domain filter"
- Check that domain names match exactly: `python3 prot_homology.py domain` and inspect `domain/pfam_results.tsv`
- Try a different domain or relax identity thresholds

### Issue: ColabFold job stuck or OOM
- Reduce `--cpus` or `--gpus` allocation
- Check job logs: `cat /scratch/.../logs/*.out`
- Monitor memory: `squeue -j <job_id> -O jobid,jobname,memory,state`

---

## Citation

If you use this pipeline, please cite:
- **DIAMOND**: Buchfink et al. (2021) Nature Methods 18:366–368
- **hmmscan**: Mistry et al. (2021) Nucleic Acids Res. 49:D212–D220
- **ColabFold**: Mirdita et al. (2022) Nat. Methods 19:679–682
- **Foldseek**: van Kempen et al. (2023) bioRxiv 2023.07.04.547681

---

## License

[Your License Here]

---

## Contact

For questions or issues, please open a GitHub issue.
