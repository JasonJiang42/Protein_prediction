#!/usr/bin/env python3
"""
prot_homology.py
----------------
Unified IgA protease homolog prediction pipeline.

Modules
-------
  clustering      Cluster pan-proteome (MMseqs2 or CD-HIT) + build DIAMOND database
  sequence        DIAMOND blastp search + identity/coverage filter
  domain          Pfam hmmscan → domain filter → alignment & species plots
  structure       --prep        Prepare ColabFold input FASTAs (before HPC)
                  --postprocess Integrate + visualize (after HPC results downloaded)
  run-all         clustering → sequence → domain → structure --prep

"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR    = Path(__file__).parent
CLUSTER_DIR = DATA_DIR / "protein_cluster"
DOMAIN_DIR  = DATA_DIR / "domain"
STRUCT_DIR  = DATA_DIR / "structure_prediction"
OUT_DIR     = DATA_DIR / "results_figure"

# Key I/O files — fixed paths (not affected by user flags)
QUERY_FAA      = DATA_DIR / "query.faa"
DEFAULT_INPUT  = DATA_DIR / "proteins_by_species" / "all_species.faa"
DOM_FILTER_TSV = DOMAIN_DIR / "function_domain_filter.tsv"
PFAM_RESULTS   = DOMAIN_DIR / "pfam_results.tsv"
STRUCT_SUM     = STRUCT_DIR / "structure_summary.tsv"
INTEGRATED_TSV = OUT_DIR / "integrated_results.tsv"

# Clustering defaults
DEFAULT_CLUST_OUTDIR  = CLUSTER_DIR
DEFAULT_CLUST_TOOL    = "mmseqs"
DEFAULT_CLUST_ID      = 0.9
DEFAULT_CLUST_QCOV    = 0.8
DEFAULT_CLUST_SCOV    = 0.8

# Sequence defaults
DEFAULT_DB            = CLUSTER_DIR / "panproteome_db"
DEFAULT_BLAST_OUT     = DATA_DIR / "query_vs_pan.tsv"
DEFAULT_SEQ_IDENTITY  = 20.0    # pident % filter for the filtered output file
DEFAULT_SEQ_QCOV      = 0.0     # DIAMOND --query-cover  (0 = no filter)
DEFAULT_SEQ_SCOV      = 0.0     # DIAMOND --subject-cover (0 = no filter)


def _renamed_faa(outdir: Path, tool: str, identity: float,
                 qcov: float, scov: float) -> Path:
    """Reconstruct the renamed FASTA path from clustering parameters."""
    if tool == "mmseqs":
        stem = f"mmseqs_c{identity}_cov{qcov}"
    else:
        stem = f"cdhit_c{identity}_aS{qcov}_aL{scov}"
    return outdir / f"{stem}.renamed.faa"


def _filtered_blast(blast_out: Path, identity: float) -> Path:
    """Derive the identity-filtered BLAST filename from the full output path."""
    return blast_out.parent / f"{blast_out.stem}_id{int(identity)}{blast_out.suffix}"

# External tools / DBs
PFAM_DB = Path("/Users/jason/Desktop/Tools/pfamdb/Pfam-A.hmm")
PYTHON  = sys.executable

# DIAMOND blastp output columns
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "qcovhsp", "scovhsp",
]

# ── Display helpers ────────────────────────────────────────────────────────────

def section(title: str) -> None:
    print(f"\n{'='*62}")
    print(f"  {title}")
    print(f"{'='*62}")


def step(title: str) -> None:
    print(f"\n{'─'*62}")
    print(f"  {title}")
    print(f"{'─'*62}")


def skipped(path: Path, label: str) -> bool:
    """Return True only if path exists AND is non-empty."""
    if path.exists() and path.stat().st_size > 0:
        print(f"  [SKIP] {label} — output exists: {path.name}")
        return True
    return False


def run(cmd: list, label: str) -> None:
    """Run a command, print it, and exit on failure."""
    step(label)
    print("  $", " ".join(str(c) for c in cmd))
    print()
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        sys.exit(f"\n  ERROR: '{label}' failed (exit code {result.returncode})")


def require(path: Path, hint: str) -> None:
    """Exit with a helpful message if a prerequisite file is missing."""
    if not path.exists():
        sys.exit(f"\n  ERROR: Required file not found:\n"
                 f"    {path}\n"
                 f"  Hint: {hint}")


# ── Module: clustering ─────────────────────────────────────────────────────────

def run_clustering(
    input_faa: Path,
    outdir: Path,
    tool: str,
    identity: float,
    qcov: float,
    scov: float,
    threads: int,
) -> tuple[Path, Path]:
    """
    Returns (renamed_faa, diamond_db) paths for downstream modules.
    """
    section("MODULE 1 — clustering")

    require(input_faa, "Provide combined species FASTA with -i / --input")
    outdir.mkdir(parents=True, exist_ok=True)

    tool_flag  = "--mmseqs" if tool == "mmseqs" else None
    tool_label = "MMseqs2" if tool == "mmseqs" else "CD-HIT"

    cmd = [
        PYTHON, str(DATA_DIR / "build_panproteome.py"),
        "-i",    str(input_faa),
        "-o",    str(outdir),
        "-c",    str(identity),
        "-aS",   str(qcov),
        "-aL",   str(scov),
        "-T",    str(threads),
    ]
    if tool_flag:
        cmd.insert(4, tool_flag)   # insert after -o value

    # build_panproteome.py has its own skip logic; always call it
    run(cmd, f"Build pan-proteome: {tool_label} cluster + DIAMOND database")

    renamed_faa = _renamed_faa(outdir, tool, identity, qcov, scov)
    diamond_db  = outdir / "panproteome_db"

    print(f"\n  Outputs:")
    print(f"    {renamed_faa}")
    print(f"    {diamond_db}.dmnd")

    return renamed_faa, diamond_db


# ── Module: sequence ───────────────────────────────────────────────────────────

def run_sequence(
    query_faa: Path,
    diamond_db: Path,
    blast_out: Path,
    identity: float,
    qcov: float,
    scov: float,
    threads: int,
    max_hits: int = 0,
) -> tuple[Path, Path]:
    """
    Returns (blast_full, blast_filtered) paths for downstream modules.
    identity / qcov / scov applied both in DIAMOND and in the post-filter step.
    """
    section("MODULE 2 — sequence (DIAMOND blastp)")

    require(query_faa, f"Query FASTA not found. Provide with -q / --query.")
    # Normalise: strip .dmnd suffix if the user included it, then re-add for the check
    diamond_db = Path(str(diamond_db).removesuffix(".dmnd"))
    if not Path(str(diamond_db) + ".dmnd").exists():
        sys.exit(f"\n  ERROR: DIAMOND database not found: {diamond_db}.dmnd\n"
                 "  Hint: Run 'clustering' module first, or pass --db <path>.")

    blast_out.parent.mkdir(parents=True, exist_ok=True)
    blast_filtered = _filtered_blast(blast_out, identity)

    # Step 2a — DIAMOND blastp
    if not skipped(blast_out, "DIAMOND blastp"):
        cmd = [
            "diamond", "blastp",
            "-q",        str(query_faa),
            "-d",        str(diamond_db),
            "--outfmt", "6", *BLAST_COLS,
            "--header",
            "-o",        str(blast_out),
            "-p",        str(threads),
            "--more-sensitive",
            "--evalue",          "1e-5",
            "--max-target-seqs", str(max_hits),
        ]
        if qcov > 0:
            cmd += ["--query-cover",   str(qcov * 100)]
        if scov > 0:
            cmd += ["--subject-cover", str(scov * 100)]
        run(cmd, f"DIAMOND blastp → {blast_out.name}")

    # Step 2b — filter to pident ≥ identity
    if not skipped(blast_filtered, f"Identity filter (pident ≥ {identity}%)"):
        step(f"Filter hits: pident ≥ {identity}% → {blast_filtered.name}")
        df    = pd.read_csv(blast_out, sep="\t", comment="#", header=None, names=BLAST_COLS)
        n_in  = len(df)
        df_f  = df[df["pident"] >= identity]
        df_f.to_csv(blast_filtered, sep="\t", index=False)
        print(f"  {n_in:,} hits → {len(df_f):,} hits (pident ≥ {identity}%)")
        print(f"  Saved: {blast_filtered.name}")

    print(f"\n  Outputs:")
    print(f"    {blast_out}")
    print(f"    {blast_filtered}")

    return blast_out, blast_filtered


# ── Pfam helpers (inlined from pfam_search.py) ────────────────────────────────

def _parse_blast_hits(blast_tsv: Path, min_pident: float) -> list[str]:
    """Return sorted list of unique sseqid where pident >= min_pident."""
    hits = set()
    with open(blast_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("qseqid"):
                continue
            cols = line.split("\t")
            try:
                if float(cols[2]) >= min_pident:
                    hits.add(cols[1])
            except (IndexError, ValueError):
                continue
    return sorted(hits)


def _extract_seqs(fasta: Path, ids: set, out_faa: Path) -> int:
    """Extract sequences whose header ID is in ids. Returns count written."""
    found, capture, count = set(), False, 0
    with open(fasta) as fh, open(out_faa, "w") as out:
        for line in fh:
            if line.startswith(">"):
                capture = line[1:].split()[0] in ids
                if capture:
                    found.add(line[1:].split()[0])
                    count += 1
            if capture:
                out.write(line)
    missing = ids - found
    if missing:
        print(f"  WARNING: {len(missing)} sseqid(s) not found in FASTA")
    return count


def _run_hmmscan(pfam_db: Path, query_faa: Path,
                 domtblout: Path, out_txt: Path, cpu: int) -> None:
    """Run hmmscan and write raw domtblout (no GA pre-filter; GA applied in post-processing)."""
    cmd = [
        "hmmscan",
        "--domtblout", str(domtblout),
        "-o",          str(out_txt),
        "--cpu",       str(cpu),
        "--noali",
        str(pfam_db), str(query_faa),
    ]
    print("  $", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-2000:], file=sys.stderr)
        sys.exit(f"\n  ERROR: hmmscan failed (exit code {result.returncode})")


DOMTBL_COLS = [
    "query_name", "query_len",
    "hmm_name", "hmm_acc", "hmm_len",
    "seq_score",
    "dom_evalue", "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "domain_coverage_pct", "description",
]


def _parse_domtblout(domtblout: Path) -> list:
    """Parse raw hmmscan domtblout, return list of row dicts (no GA filter, no source column)."""
    rows = []
    with open(domtblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 23:
                continue
            hmm_from, hmm_to, hmm_len = cols[15], cols[16], cols[2]
            try:
                cov = round((int(hmm_to) - int(hmm_from) + 1) / int(hmm_len) * 100, 1)
            except (ValueError, ZeroDivisionError):
                cov = ""
            rows.append(dict(zip(DOMTBL_COLS, [
                cols[3], cols[5],
                cols[0], cols[1], hmm_len,
                cols[7],                                          # seq_score (full-seq)
                cols[11], cols[13], cols[14],
                hmm_from, hmm_to, cols[17], cols[18], cols[19], cols[20],
                cov, " ".join(cols[22:]),
            ])))
    return rows


def _parse_ga_thresholds(pfam_db: Path) -> dict:
    """Return {hmm_name: (seq_ga, dom_ga)} parsed from GA lines in Pfam-A.hmm."""
    thresholds: dict = {}
    name = None
    with open(pfam_db) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("GA ") and name:
                parts = line.split()
                thresholds[name] = (float(parts[1]), float(parts[2].rstrip(";")))
                name = None
    return thresholds


def _ga_filter_rows(rows: list, ga_thresholds: dict) -> list:
    """Keep rows where seq_score >= seq_GA and dom_score >= dom_GA."""
    kept = []
    for row in rows:
        thresh = ga_thresholds.get(row["hmm_name"])
        if thresh is None:
            continue
        seq_ga, dom_ga = thresh
        try:
            if float(row["seq_score"]) >= seq_ga and float(row["dom_score"]) >= dom_ga:
                kept.append(row)
        except (ValueError, TypeError):
            continue
    return kept


# ── Module: domain-filter ──────────────────────────────────────────────────────

def _overlap(qstart: int, qend: int, dom_from: int, dom_to: int) -> int:
    """Return number of overlapping residues (0 = no overlap)."""
    return max(0, min(qend, dom_to) - max(qstart, dom_from) + 1)


def run_domain_filter(
    domains: list,
    pfam_tsv: Path,
    blast_tsv: Path,
    out_tsv: Path,
) -> None:
    """
    Filter blast hits by requiring the query alignment region to overlap a
    user-specified Pfam domain found in the query sequences.

    domains  : list of Pfam hmm_name (e.g. 'Peptidase_M64') or
               hmm_acc (e.g. 'PF03009') — both columns are checked.
    pfam_tsv : combined pfam_results.tsv with a 'source' column (query/subject).
               Query rows are already GA-filtered from the 'domain' module.
    blast_tsv: full or identity-filtered DIAMOND blast output.
    out_tsv  : output path for filtered hits TSV.
    """
    section("MODULE: domain-filter")

    require(pfam_tsv,  "Run 'domain' module first to generate pfam_results.tsv.")
    require(blast_tsv, "Provide BLAST TSV with --blast.")

    df = pd.read_csv(pfam_tsv, sep="\t")
    df_q = df[df["source"] == "query"]
    df_q = df_q[df_q["hmm_name"].isin(domains) | df_q["hmm_acc"].isin(domains)]

    if df_q.empty:
        sys.exit(
            f"\n  ERROR: no query domain regions found for: {domains}\n"
            f"  Check domain names/accessions against pfam_results.tsv columns "
            f"'hmm_name' and 'hmm_acc'."
        )

    # Build {query_name: [(hmm_name, env_from, env_to), ...]}
    query_regions: dict = {}
    for _, row in df_q.iterrows():
        query_regions.setdefault(row["query_name"], []).append(
            (row["hmm_name"], int(row["env_from"]), int(row["env_to"]))
        )

    print(f"  Domains matched in {len(query_regions)} query sequence(s):")
    for qid, regions in query_regions.items():
        for dom, ef, et in regions:
            print(f"    {qid}  {dom}  {ef}–{et} aa")

    blast = pd.read_csv(blast_tsv, sep="\t")
    kept_rows = []
    for _, row in blast.iterrows():
        qid = row["qseqid"]
        if qid not in query_regions:
            continue
        best_dom, best_ov = None, 0
        for dom_name, env_from, env_to in query_regions[qid]:
            ov = _overlap(int(row["qstart"]), int(row["qend"]), env_from, env_to)
            if ov > best_ov:
                best_ov, best_dom = ov, dom_name
        if best_ov > 0:
            r = row.to_dict()
            r["matched_domain"] = best_dom
            r["overlap_aa"]     = best_ov
            kept_rows.append(r)

    extra_cols = ["matched_domain", "overlap_aa"]
    if kept_rows:
        result = pd.DataFrame(kept_rows)[list(blast.columns) + extra_cols]
    else:
        result = pd.DataFrame(columns=list(blast.columns) + extra_cols)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_tsv, sep="\t", index=False)

    print(f"\n  Blast input : {len(blast)} hits")
    print(f"  After filter: {len(result)} hits  ({len(blast) - len(result)} removed)")
    for qid in query_regions:
        sub = result[result["qseqid"] == qid]
        if len(sub):
            print(f"    {qid}: {len(sub)} hits — domains: {sub['matched_domain'].unique().tolist()}")
    print(f"\n  Saved: {out_tsv}")


# ── Module: domain ─────────────────────────────────────────────────────────────

def run_domain(
    subject_faa: Path,
    threads: int,
    pfam_db: Path,
    dom_outdir: Path,
    query_faa: Path,
) -> None:
    section("MODULE 3 — domain (Pfam annotation)")

    require(subject_faa, "Subject FASTA not found. Provide with --subject.")
    require(query_faa,   "Query FASTA not found. Provide with --query.")
    if not pfam_db.exists():
        sys.exit(f"\n  ERROR: Pfam-A.hmm not found at:\n    {pfam_db}\n"
                 "  Download from https://pfam.xfam.org/ and hmmpress it.")

    dom_outdir.mkdir(parents=True, exist_ok=True)

    pfam_results    = dom_outdir / "pfam_results.tsv"
    subj_domtblout  = dom_outdir / "subject.domtblout"
    subj_hmm_txt    = dom_outdir / "subject.hmm.out"
    query_domtblout = dom_outdir / "query.domtblout"
    query_hmm_txt   = dom_outdir / "query.hmm.out"

    # Step 3a — hmmscan on subject sequences
    if not skipped(subj_domtblout, "hmmscan: subject sequences"):
        step("Pfam hmmscan → subject.domtblout")
        print(f"  Running hmmscan on {subject_faa.name} ({threads} CPUs) ...")
        _run_hmmscan(pfam_db, subject_faa, subj_domtblout, subj_hmm_txt, threads)

    # Step 3b — hmmscan on query sequences
    if not skipped(query_domtblout, "hmmscan: query sequences"):
        step("Pfam hmmscan → query.domtblout")
        print(f"  Running hmmscan ({threads} CPUs) ...")
        _run_hmmscan(pfam_db, query_faa, query_domtblout, query_hmm_txt, threads)

    # Step 3c — parse raw domtblouts, apply GA filter, combine → pfam_results.tsv
    require(subj_domtblout,  "subject.domtblout missing — delete and rerun 'domain'.")
    require(query_domtblout, "query.domtblout missing — delete and rerun 'domain'.")
    step(f"GA filter + combine → {pfam_results.name}")
    print(f"  Loading GA thresholds from {pfam_db.name} ...")
    ga = _parse_ga_thresholds(pfam_db)
    print(f"  {len(ga)} HMM profiles with GA thresholds")

    rows_subj_raw  = _parse_domtblout(subj_domtblout)
    rows_query_raw = _parse_domtblout(query_domtblout)
    rows_subj  = _ga_filter_rows(rows_subj_raw,  ga)
    rows_query = _ga_filter_rows(rows_query_raw, ga)
    print(f"  subject : {len(rows_subj_raw)} raw → {len(rows_subj)} GA-passed")
    print(f"  query   : {len(rows_query_raw)} raw → {len(rows_query)} GA-passed")

    df_subj  = pd.DataFrame(rows_subj,  columns=DOMTBL_COLS) if rows_subj  else pd.DataFrame(columns=DOMTBL_COLS)
    df_query = pd.DataFrame(rows_query, columns=DOMTBL_COLS) if rows_query else pd.DataFrame(columns=DOMTBL_COLS)
    df_subj["source"]  = "subject"
    df_query["source"] = "query"
    df_all = pd.concat([df_query, df_subj], ignore_index=True)
    df_all.to_csv(pfam_results, sep="\t", index=False)

    if df_subj.empty:
        print("  WARNING: no subject domain hits passed GA thresholds.")
    if df_query.empty:
        print("  WARNING: no query domain hits passed GA thresholds.")

    print(f"\n  Outputs:")
    print(f"    {subj_domtblout.name} / {query_domtblout.name}  (raw hmmscan domtblout)")
    print(f"    {subj_hmm_txt.name} / {query_hmm_txt.name}    (raw hmmscan text output)")
    print(f"    {pfam_results.name}  (GA-filtered, source=query|subject)")
    print(f"\n  Next: python3 prot_homology.py domain-filter --domains <NAME,...>")


# ── Module: structure — prep ───────────────────────────────────────────────────


def run_structure_prep(fasta: Path, output: Path) -> None:
    section("MODULE 4a — structure: ColabFold API MSA")

    require(fasta, "Provide input FASTA with --fasta")
    output.mkdir(parents=True, exist_ok=True)

    run(
        ["colabfold_batch", "--msa-only", str(fasta), str(output)],
        f"colabfold_batch API MSA only: {fasta.name} → {output}/",
    )

    print(f"\n{'─'*62}")
    print("  HPC submission (agis-hpc):")
    print(f"{'─'*62}")
    print("  1. Transfer input FASTAs to server:")
    print("       scp -r structure_prediction/Q*/input.faa \\")
    print("           jiangshuo@ln01:/scratch/jiangshuo/<project>/input/")
    print()
    print("  2. Submit ColabFold array job:")
    print("       sbatch structure_prediction/submit_colabfold.sh")
    print()
    print("  3. Run Foldseek (after ColabFold completes):")
    print("       sbatch structure_prediction/submit_foldseek.sh")
    print()
    print("  4. Collect pLDDT results on server:")
    print("       python3 structure_prediction/collect_results.py")
    print()
    print("  5. Download results to local machine:")
    print("       scp -r jiangshuo@ln01:/scratch/jiangshuo/.../output/ \\")
    print("           structure_prediction/")
    print("       scp jiangshuo@ln01:.../structure_prediction/structure_summary.tsv \\")
    print("           structure_prediction/")
    print("       scp jiangshuo@ln01:.../foldseek/*_results.tsv \\")
    print("           structure_prediction/foldseek/")
    print()
    print("  6. Then run:")
    print("       python3 prot_homology.py structure --postprocess")
    print(f"{'─'*62}")


def run_structure_predict(
    msa_dir: Path,
    output: Path,
    weights: Path,
    gpus: int = 1,
    cpus: int = 8,
    partition: str = "NV_4090D",
    job_name: str = "colabfold",
    array_jobs: int | None = None,
    array_concurrent: int = 0,
) -> None:
    """Generate and submit ColabFold structure prediction SLURM job."""
    section("MODULE 4b — structure: SUBMIT ColabFold prediction job")

    require(msa_dir, "Provide MSA input folder with --input")
    require(weights, "Provide weights DB path with --weights")
    output.mkdir(parents=True, exist_ok=True)

    # Create logs directory
    logs_dir = Path.cwd() / "logs"
    logs_dir.mkdir(exist_ok=True)

    # Generate SLURM script
    script_path = Path.cwd() / f"submit_{job_name}.sh"
    mem_per_cpu = 7000 if partition == "NV_4090D" else 3500

    # Build array directive if specified
    array_directive = ""
    if array_jobs:
        array_directive = f"#SBATCH --array=0-{array_jobs - 1}"
        if array_concurrent > 0:
            array_directive += f"%{array_concurrent}"
        array_directive += "\n"

    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --gres=gpu:{gpus}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem-per-cpu={mem_per_cpu}M
{array_directive}#SBATCH --output={logs_dir.absolute()}/%j.out
#SBATCH --error={logs_dir.absolute()}/%j.err

source $(conda info --base)/etc/profile.d/conda.sh
conda activate colabfold

echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Start: $(date)"

mkdir -p {output.absolute()}

colabfold_batch \\
  --num-recycle 3 \\
  --data {weights.absolute()} \\
  {msa_dir.absolute()} \\
  {output.absolute()}

echo "Done: $(date)"
"""

    script_path.write_text(slurm_script)
    print(f"\n  Generated SLURM script: {script_path.name}")
    print(f"    - Partition: {partition}")
    print(f"    - GPUs: {gpus}")
    print(f"    - CPUs: {cpus}")
    print(f"    - Memory/CPU: {mem_per_cpu}M")
    if array_jobs:
        max_concurrent = f" (max {array_concurrent} concurrent)" if array_concurrent > 0 else ""
        print(f"    - Array jobs: 0-{array_jobs - 1}{max_concurrent}")

    # Submit job
    step("Submitting SLURM job")
    result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: sbatch submission failed")
        print(result.stderr)
        sys.exit(1)

    job_id = result.stdout.strip().split()[-1]
    print(f"\n  ✓ Job submitted: {job_id}")
    print(f"    Monitor with: squeue -j {job_id}")
    print(f"    Output dir: {output}/")


def run_structure_align(
    query: Path | None,
    query_dir: Path | None,
    targets: Path,
    output: Path,
    threshold: float = 0.5,
    threads: int = 8,
) -> None:
    """Run Foldseek structural alignment."""
    section("Foldseek: Structural Alignment")

    require(targets, "Provide target directory with --targets")

    # Collect query PDB files
    queries = {}
    if query:
        require(query, "Query PDB file not found")
        name = query.stem
        queries[name] = query

    if query_dir:
        require(query_dir, "Query directory not found")
        for pdb in sorted(query_dir.glob("*.pdb")):
            queries[pdb.stem] = pdb

    if not queries:
        sys.exit("ERROR: Provide either --query or --query-dir")

    # Collect target PDB files
    target_pdbs = list(targets.glob("*.pdb"))
    if not target_pdbs:
        sys.exit(f"ERROR: No PDB files found in {targets}")

    print(f"  Queries   : {len(queries)}")
    print(f"  Targets   : {len(target_pdbs)}")
    print(f"  Threshold : TM-score >= {threshold}")

    # Create temporary directory for foldseek
    tmp_dir = output.parent / ".foldseek_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # TSV header
    header = "query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqtmscore\tttmscore\talntmscore\tlddt"

    # Run foldseek for each query
    all_results = []
    for qname, qpdb in sorted(queries.items()):
        step(f"Aligning {qname} against targets")

        # Create temp dir for this query
        q_tmp = tmp_dir / qname
        q_tmp.mkdir(exist_ok=True)

        # Run foldseek easy-search with conda activation
        conda_init = "source $(conda info --base)/etc/profile.d/conda.sh && conda activate foldseek &&"
        cmd_str = (
            f"{conda_init} foldseek easy-search "
            f"{str(qpdb)} {str(targets)} {str(q_tmp / 'results.tsv')} {str(q_tmp / 'search')} "
            f"--format-output {header.replace(chr(9), ',')} "
            f"--threads {threads} -e 10 --exhaustive-search 1"
        )

        result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: foldseek failed for {qname}")
            print(result.stderr)
            continue

        # Read results and filter by threshold
        results_file = q_tmp / "results.tsv"
        if results_file.exists():
            with open(results_file) as f:
                for line in f:
                    fields = line.strip().split("\t")
                    if len(fields) >= 14:
                        try:
                            tm_score = float(fields[14]) if fields[14] else 0.0
                            if tm_score >= threshold:
                                all_results.append(line.rstrip())
                        except (ValueError, IndexError):
                            pass

            # Count hits
            with open(results_file) as f:
                total = len(f.readlines())
            filtered = sum(1 for r in all_results if r.split("\t")[0] == qname)
            print(f"    Hits: {total} total, {filtered} >= {threshold}")

    # Write combined results
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as f:
        f.write(header + "\n")
        for result in all_results:
            f.write(result + "\n")

    # Summary
    print(f"\n  ✓ Results written: {output}")
    print(f"    Total hits >= {threshold}: {len(all_results)}")

    # Cleanup
    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)


# ── Module: integrate ──────────────────────────────────────────────────────────

def run_integrate(
    blast_tsv: Path,
    foldseek_tsv: Path,
    output: Path,
) -> None:
    """
    Merge BLAST (qseqid/sseqid) with Foldseek alignment on matching query+subject pair.

    Join keys:
      blast qseqid  ↔  foldseek query  (strip _unrelaxed/_relaxed suffix)
      blast sseqid  ↔  foldseek target (extract UniRef accession via regex)
    """
    section("MODULE: integrate")

    require(blast_tsv,    "Provide BLAST results with --blast")
    require(foldseek_tsv, "Provide Foldseek results with --foldseek")

    # ── 1. Load BLAST ─────────────────────────────────────────────────────────
    step("Load BLAST results")
    blast = pd.read_csv(blast_tsv, sep="\t", comment="#")
    blast = blast.loc[:, ~blast.columns.str.startswith("Unnamed")]
    print(f"  BLAST hits: {len(blast)}")

    # ── 2. Load Foldseek ──────────────────────────────────────────────────────
    step("Load Foldseek alignment results")
    fs = pd.read_csv(foldseek_tsv, sep="\t", header=0)

    # Extract qseqid from foldseek query column:
    #   Hungatella_hathewayi_DSM_13479_unrelaxed_rank_001_... → Hungatella_hathewayi_DSM_13479
    fs["qseqid"] = fs["query"].str.split("_unrelaxed|_relaxed").str[0]

    # Extract sseqid from foldseek target column:
    #   UniRef90_A0A174XAW7_Alpha_beta_..._unrelaxed_... → UniRef90_A0A174XAW7
    # ColabFold encodes the full FASTA description into the filename; the bare
    # accession (blast sseqid format) is just the UniRef prefix + alphanumeric ID.
    fs["sseqid"] = fs["target"].str.extract(
        r'^((?:UniRef\d+|UniProt|tr|sp)_[A-Z0-9]+)', expand=False
    )
    # Fallback for non-UniRef targets: strip model suffix only
    fallback = fs["sseqid"].isna()
    fs.loc[fallback, "sseqid"] = (
        fs.loc[fallback, "target"].str.split("_unrelaxed|_relaxed").str[0]
    )

    # Rename overlapping columns to avoid clash with BLAST
    fs = fs.rename(columns={
        "query":    "fs_query",
        "target":   "fs_target",
        "mismatch": "fs_mismatch",
        "gapopen":  "fs_gapopen",
        "qstart":   "fs_qstart",
        "qend":     "fs_qend",
        "evalue":   "fs_evalue",
    })
    print(f"  Foldseek hits: {len(fs)}")
    print(f"  Sample qseqid (fs): {fs['qseqid'].head(3).tolist()}")
    print(f"  Sample sseqid (fs): {fs['sseqid'].head(3).tolist()}")

    # ── 3. Merge on (qseqid, sseqid) ─────────────────────────────────────────
    step("Merge BLAST + Foldseek on (qseqid, sseqid)")
    merged = blast.merge(fs, on=["qseqid", "sseqid"], how="left")
    n_matched = merged["alntmscore"].notna().sum()
    print(f"  After merge: {len(merged)} rows  ({n_matched} with Foldseek data)")

    # ── 4. Save ───────────────────────────────────────────────────────────────
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, sep="\t", index=False, float_format="%.4g")
    print(f"\n  ✓ Integrated table: {len(merged)} rows × {len(merged.columns)} columns")
    print(f"  ✓ Saved: {output}")


# ── Module: structure — postprocess ───────────────────────────────────────────

def run_structure_postprocess() -> None:
    section("MODULE 4b — structure: POSTPROCESS (integrate + visualize)")

    require(STRUCT_SUM,
            "Run collect_results.py on the server, then download structure_summary.tsv.")
    require(DOM_FILTER_TSV,
            "Run 'domain' module first.")

    OUT_DIR.mkdir(exist_ok=True)

    # Step 4b-1 — Integrate BLAST + Pfam + Foldseek
    run([PYTHON, str(DATA_DIR / "integrate_results.py")],
        "Integrate BLAST + domain + Foldseek → integrated_results.tsv")

    print(f"\n  Outputs:")
    print(f"    {INTEGRATED_TSV.relative_to(DATA_DIR)}")
    print(f"\n  Next: python3 prot_homology.py plot --foldseek  (generate Foldseek plots)")


# ── Plot helpers ──────────────────────────────────────────────────────────────

_ABBREV = {
    "Escherichia_coli":              "E. coli",
    "Klebsiella_pneumoniae":         "K. pneumoniae",
    "Lactobacillus_johnsonii":       "L. johnsonii",
    "Lacticaseibacillus_rhamnosus":  "L. rhamnosus",
    "Lactiplantibacillus_plantarum": "L. plantarum",
    "Segatella_copri":               "S. copri",
    "Bacteroides_fragilis":          "B. fragilis",
    "Bacteroides_thetaiotaomicron":  "B. thetaiotaomicron",
    "Bacteroides_ovatus":            "B. ovatus",
    "Alistipes_putredinis":          "A. putredinis",
    "Alistipes_finegoldii":          "A. finegoldii",
    "Alistipes_shahii":              "A. shahii",
    "Bifidobacterium_bifidum":       "B. bifidum",
    "Bifidobacterium_breve":         "B. breve",
    "Bifidobacterium_longum":        "B. longum",
    "Enterococcus_faecalis":         "E. faecalis",
    "Enterococcus_faecium":          "E. faecium",
}


def _fasta_lengths(fasta: Path) -> dict:
    """Return {seq_id: length} for all sequences in a FASTA file."""
    lengths: dict = {}
    seq_id, length = None, 0
    with open(fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id:
                    lengths[seq_id] = length
                seq_id = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if seq_id:
        lengths[seq_id] = length
    return lengths


def _hit_label(sseqid: str) -> str:
    """Format sseqid as 'accession  |  Sp. name'."""
    parts   = sseqid.split("|")
    acc     = parts[1] if len(parts) > 1 else sseqid
    species = parts[2] if len(parts) > 2 else ""
    sp_abbr = _ABBREV.get(species, species.replace("_", " "))
    return f"{acc}  |  {sp_abbr}"


def _plot_blast(blast_tsv: Path, query_faa: Path, out_dir: Path) -> None:
    """
    BLAST alignment coverage bars, one figure per query.
    Adapted from plot_blast.py — query lengths read from query FASTA.
    """
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
    except ImportError:
        sys.exit("  ERROR: matplotlib / numpy required.  pip install matplotlib numpy")

    require(blast_tsv, "Run 'sequence' module first.")
    require(query_faa, f"Query FASTA not found: {query_faa}")

    query_lens = _fasta_lengths(query_faa)
    df = pd.read_csv(blast_tsv, sep="\t")
    df["label"] = df["sseqid"].apply(_hit_label)

    vmin, vmax = 20, 65
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "blues_clipped", cm.Blues(np.linspace(0.3, 1.0, 256))
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    for qid, qlen in query_lens.items():
        sub = df[df["qseqid"] == qid].copy()
        sub = sub.sort_values(["sseqid", "pident"], ascending=[True, False])
        if sub.empty:
            print(f"  [skip] no hits for {qid}")
            continue

        n  = len(sub)
        ys = list(range(n))
        fig_h = max(6, n * 0.18 + 2)
        fig, ax = plt.subplots(figsize=(10, fig_h))

        # grey full-length background bars
        ax.barh(ys, [qlen] * n, left=0, height=0.65,
                color="#e8e8e8", zorder=1, linewidth=0)
        # coloured alignment bars
        for y, (_, row) in zip(ys, sub.iterrows()):
            width = row["qend"] - row["qstart"] + 1
            ax.barh(y, width, left=row["qstart"] - 1, height=0.65,
                    color=cmap(norm(row["pident"])), zorder=2, linewidth=0)

        ylabel_size = 5.5 if n > 50 else 7.5
        ax.set_yticks(ys)
        ax.set_yticklabels(sub["label"].tolist(), fontsize=ylabel_size)
        ax.set_xlim(0, qlen * 1.02)
        ax.set_ylim(-0.8, n - 0.2)
        ax.set_xlabel("Query position (aa)", fontsize=10)
        ax.set_title(
            f"DIAMOND blastp — {qid}  ({qlen} aa)\n"
            f"hits ≥ 20% identity  |  n = {n}",
            fontsize=11, fontweight="bold", pad=10,
        )
        ax.axvline(qlen, color="#aaaaaa", linewidth=0.8, linestyle="--", zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", labelsize=9)

        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                            fraction=0.015, pad=0.02, aspect=30)
        cbar.set_label("Sequence identity (%)", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

        plt.tight_layout()
        for ext in ("pdf", "png"):
            out = out_dir / f"blast_alignment_plot_{qid}.{ext}"
            fig.savefig(out, dpi=200, bbox_inches="tight")
            print(f"    Saved: {out.name}")
        plt.close(fig)


def _plot_domain(
    dom_filter_tsv: Path,
    pfam_tsv: Path,
    query_domtblout: Path,
    subject_faa: Path,
    query_faa: Path,
    out_dir: Path,
) -> None:
    """
    Domain architecture plots via gggenes_alignment.R.
    Requires R with packages: gggenes, ggplot2, dplyr, readr.
    """
    r_script = DATA_DIR / "gggenes_alignment.R"
    require(r_script,        "gggenes_alignment.R not found in the project root directory.")
    require(dom_filter_tsv,  "Run 'domain-filter' module first.")
    require(pfam_tsv,        "Run 'domain' module first.")
    require(query_domtblout, "Run 'domain' module first.")
    require(subject_faa,     f"Subject FASTA not found: {subject_faa}")
    require(query_faa,       f"Query FASTA not found: {query_faa}")

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "Rscript", "--vanilla", str(r_script),
        str(dom_filter_tsv),
        str(pfam_tsv),
        str(query_domtblout),
        str(subject_faa),
        str(query_faa),
        str(out_dir),
    ]
    step("Domain architecture plots (gggenes)")
    print("  $", " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        sys.exit(f"\n  ERROR: Rscript failed (exit code {result.returncode})")


def _plot_foldseek(
    integrated_tsv: Path,
    out_dir: Path,
    blast_threshold: float = 30.0,
    align_threshold: float = 0.5,
    lddt_threshold: float = 0.5,
) -> None:
    """
    Integrated scatter: x=BLAST pident, y=Foldseek alntmscore, size=lDDT.
    Star markers where all three criteria are met: pident, alntmscore, lddt.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("  ERROR: matplotlib required.")

    require(integrated_tsv, "Run 'integrate' module first.")
    df = pd.read_csv(integrated_tsv, sep="\t")
    df = df[df["alntmscore"].notna()].copy()

    for col in ("pident", "alntmscore"):
        if col not in df.columns:
            sys.exit(f"\n  ERROR: column '{col}' missing in {integrated_tsv.name}.\n"
                     f"  Available: {list(df.columns)}")

    has_plddt = "lddt" in df.columns
    if not has_plddt:
        print("  [INFO] lddt not found — using fixed dot size")

    out_dir.mkdir(parents=True, exist_ok=True)
    queries  = df["qseqid"].unique().tolist() if "qseqid" in df.columns else [None]
    palette  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                "#9467bd", "#8c564b", "#17becf", "#e377c2"]
    col_map  = {q: palette[i % len(palette)] for i, q in enumerate(queries)}

    # axis ranges from data
    x_pad = max((df["pident"].max() - df["pident"].min()) * 0.05, 2)
    y_pad = max((df["alntmscore"].max() - df["alntmscore"].min()) * 0.05, 0.02)
    x_min = df["pident"].min()     - x_pad
    x_max = df["pident"].max()     + x_pad
    y_min = max(0.0, df["alntmscore"].min() - y_pad)
    y_max = min(1.02, df["alntmscore"].max() + y_pad)

    fig, ax = plt.subplots(figsize=(9, 6.5))

    for qid in queries:
        sub = df[df["qseqid"] == qid] if qid is not None else df
        mask_hi = (
            (sub["pident"]     >= blast_threshold) &
            (sub["alntmscore"] >= align_threshold) &
            (sub["lddt"]       >= lddt_threshold   if has_plddt else True)
        )
        sub_reg = sub[~mask_hi]
        sub_hi  = sub[mask_hi]
        c = col_map.get(qid, "#555555")
        label = f"{qid}  (n={len(sub)})" if qid is not None else f"All  (n={len(sub)})"

        s_reg = (sub_reg["lddt"] * 300).clip(lower=20) if has_plddt else 40
        s_hi  = (sub_hi["lddt"]  * 300).clip(lower=20) if has_plddt else 40

        ax.scatter(sub_reg["pident"], sub_reg["alntmscore"],
                   c=c, s=s_reg, alpha=0.65,
                   edgecolors="white", linewidths=0.4,
                   label=label, zorder=3)
        if len(sub_hi):
            ax.scatter(sub_hi["pident"], sub_hi["alntmscore"],
                       c=c, s=s_hi * 2.5 if has_plddt else 200,
                       alpha=0.95, edgecolors="black", linewidths=1.2,
                       marker="*", zorder=5)

    # Threshold lines
    ax.axvline(blast_threshold, color="#555555", linewidth=1.2, linestyle="--", zorder=2)
    ax.axhline(align_threshold, color="#888888", linewidth=1.2, linestyle=":",  zorder=2)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    lddt_clause = f" & lDDT≥{lddt_threshold}" if has_plddt else ""
    hi_mask = (
        (df["pident"]     >= blast_threshold) &
        (df["alntmscore"] >= align_threshold) &
        (df["lddt"]       >= lddt_threshold   if has_plddt else True)
    )
    n_hi = int(hi_mask.sum())
    ax.scatter([], [], marker="*", c="grey", s=140, edgecolors="black",
               linewidths=1.0,
               label=f"★ pident≥{blast_threshold}% & alntmscore≥{align_threshold}"
                     f"{lddt_clause}  (n={n_hi})")
    handles, labels = ax.get_legend_handles_labels()
    size_note = "  |  size ∝ lDDT" if has_plddt else ""
    ax.legend(handles=handles, labels=labels, fontsize=8.5, frameon=False,
              loc="lower right",
              title=f"Query group{size_note}",
              title_fontsize=7.5)

    ax.set_xlabel("BLAST sequence identity (%)", fontsize=11)
    ax.set_ylabel("Foldseek alntmscore", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = out_dir / f"foldseek_scatter.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"    Saved: {out.name}")
    plt.close(fig)


# ── Module: plot ───────────────────────────────────────────────────────────────

def run_plot(
    blast = None,        # Path to identity-filtered BLAST TSV, or None to skip
    domain = None,       # Path to domain output directory, or None to skip
    foldseek = None,     # Path to integrated results TSV, or None to skip
    subject_faa = None,  # Subject FASTA — required for --domain plot
    query_faa: Path = QUERY_FAA,
    out_dir: Path = OUT_DIR,
    blast_threshold: float = 30.0,
    align_threshold: float = 0.5,
    lddt_threshold: float = 0.5,
) -> None:
    section("MODULE: plot")

    if not any([blast, domain, foldseek]):
        sys.exit(
            "\n  ERROR: no input provided.\n"
            "  Provide one or more of: --blast TSV  --domain DIR  --foldseek TSV"
        )

    out_dir.mkdir(parents=True, exist_ok=True)

    if blast:
        step("BLAST alignment coverage bars")
        _plot_blast(
            blast_tsv = Path(blast),
            query_faa = query_faa,
            out_dir   = out_dir,
        )

    if domain:
        if not subject_faa:
            sys.exit(
                "\n  ERROR: --domain plot requires --subject <FAA>\n"
                "  Provide the same subject FASTA used in the 'domain' module."
            )
        dom_dir = Path(domain)
        _plot_domain(
            dom_filter_tsv  = dom_dir / "function_domain_filter.tsv",
            pfam_tsv        = dom_dir / "pfam_results.tsv",
            query_domtblout = dom_dir / "query.domtblout",
            subject_faa     = Path(subject_faa),
            query_faa       = query_faa,
            out_dir         = out_dir,
        )

    if foldseek:
        step("Foldseek scatter plot")
        _plot_foldseek(
            integrated_tsv  = Path(foldseek),
            out_dir         = out_dir,
            blast_threshold = blast_threshold,
            align_threshold = align_threshold,
            lddt_threshold  = lddt_threshold,
        )

    print(f"\n  Output directory: {out_dir}/")


# ── Module: run-all ────────────────────────────────────────────────────────────

def run_all(args) -> None:
    section("MODULE: run-all")

    # clustering
    renamed_faa, diamond_db = run_clustering(
        input_faa = Path(args.input),
        outdir    = Path(args.clust_output),
        tool      = args.tool,
        identity  = args.clust_identity,
        qcov      = args.clust_qcov,
        scov      = args.clust_scov,
        threads   = args.threads,
    )

    # sequence
    blast_full, blast_filtered = run_sequence(
        query_faa  = Path(args.query),
        diamond_db = Path(args.db) if args.db else diamond_db,
        blast_out  = Path(args.seq_output),
        identity   = args.seq_identity,
        qcov       = args.seq_qcov,
        scov       = args.seq_scov,
        threads    = args.threads,
        max_hits   = args.max_hits,
    )

    # domain
    run_domain(
        subject_faa = renamed_faa,
        threads     = args.threads,
        pfam_db     = Path(args.pfam),
        dom_outdir  = Path(args.dom_outdir),
        query_faa   = Path(args.query),
    )

    # blast plot only — domain needs domain-filter first; foldseek needs HPC results
    run_plot(blast=blast_filtered, query_faa=Path(args.query), out_dir=OUT_DIR)

    # structure prep — skipped in run-all (no fasta/output known at this stage)

    print(f"\n{'='*62}")
    print("  run-all complete (local steps done).")
    print("  Follow the HPC instructions above, then run:")
    print("    python3 prot_homology.py structure --postprocess")
    print("    python3 prot_homology.py plot --foldseek")
    print(f"{'='*62}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def _add_clustering_args(p) -> None:
    p.add_argument(
        "-i", "--input", default=str(DEFAULT_INPUT), metavar="FAA",
        help=f"Combined all-species FASTA (default: {DEFAULT_INPUT.name})",
    )
    p.add_argument(
        "--clust-output", default=str(DEFAULT_CLUST_OUTDIR), metavar="DIR",
        help=f"Output directory (default: {DEFAULT_CLUST_OUTDIR.name}/)",
    )
    p.add_argument(
        "--tool", choices=["mmseqs", "cdhit"], default=DEFAULT_CLUST_TOOL,
        help="Clustering tool: mmseqs (default) or cdhit",
    )
    p.add_argument(
        "-c", "--clust-identity", type=float, default=DEFAULT_CLUST_ID,
        metavar="FLOAT",
        help=f"Sequence identity threshold 0–1 (default: {DEFAULT_CLUST_ID})",
    )
    p.add_argument(
        "--clust-qcov", dest="clust_qcov",
        type=float, default=DEFAULT_CLUST_QCOV, metavar="FLOAT",
        help=f"Query coverage threshold 0–1 (default: {DEFAULT_CLUST_QCOV})",
    )
    p.add_argument(
        "--clust-scov", dest="clust_scov",
        type=float, default=DEFAULT_CLUST_SCOV, metavar="FLOAT",
        help=f"Subject coverage threshold 0–1 (default: {DEFAULT_CLUST_SCOV})",
    )


def _add_thread_arg(p) -> None:
    p.add_argument(
        "-T", "--threads", type=int, default=8,
        help="CPU threads (default: 8)",
    )


def _add_sequence_args(p) -> None:
    p.add_argument(
        "-q", "--query", default=str(QUERY_FAA), metavar="FAA",
        help=f"Query FASTA file (default: {QUERY_FAA.name})",
    )
    p.add_argument(
        "--db", default=None, metavar="PATH",
        help=f"DIAMOND database path, without .dmnd suffix "
             f"(default: {DEFAULT_DB.name} inside clustering output dir)",
    )
    p.add_argument(
        "--seq-output", default=str(DEFAULT_BLAST_OUT), metavar="TSV",
        help=f"BLAST output TSV (default: {DEFAULT_BLAST_OUT.name})",
    )
    p.add_argument(
        "--identity", dest="seq_identity",
        type=float, default=DEFAULT_SEQ_IDENTITY, metavar="PCT",
        help=f"Min %% identity for the filtered output file (default: {DEFAULT_SEQ_IDENTITY})",
    )
    p.add_argument(
        "--seq-qcov", dest="seq_qcov",
        type=float, default=DEFAULT_SEQ_QCOV, metavar="FLOAT",
        help=f"Min query coverage passed to DIAMOND 0–1 (default: {DEFAULT_SEQ_QCOV}, off)",
    )
    p.add_argument(
        "--seq-scov", dest="seq_scov",
        type=float, default=DEFAULT_SEQ_SCOV, metavar="FLOAT",
        help=f"Min subject coverage passed to DIAMOND 0–1 (default: {DEFAULT_SEQ_SCOV}, off)",
    )
    p.add_argument(
        "--max-hits", dest="max_hits", type=int, default=0, metavar="N",
        help="Max hits per query passed to DIAMOND --max-target-seqs "
             "(default: 0 = unlimited; DIAMOND's built-in default is 25)",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="prot_homology.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    sub = parser.add_subparsers(
        dest="module", required=True,
        metavar="{clustering,sequence,domain,domain-filter,plot,structure,run-all}",
    )

    # ── clustering ──
    p_clust = sub.add_parser(
        "clustering",
        help="Cluster pan-proteome (MMseqs2 or CD-HIT) + build DIAMOND database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_clustering_args(p_clust)
    _add_thread_arg(p_clust)

    # ── sequence ──
    p_seq = sub.add_parser(
        "sequence",
        help="DIAMOND blastp + identity/coverage filter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_sequence_args(p_seq)
    _add_thread_arg(p_seq)

    # ── domain ──
    p_func = sub.add_parser(
        "domain",
        help="Pfam domain annotation (hmmscan → pfam_results.tsv)",
    )
    p_func.add_argument(
        "--subject", default=str(_renamed_faa(DEFAULT_CLUST_OUTDIR,
                                              DEFAULT_CLUST_TOOL,
                                              DEFAULT_CLUST_ID,
                                              DEFAULT_CLUST_QCOV,
                                              DEFAULT_CLUST_SCOV)),
        metavar="FAA",
        help="Subject (pan-proteome) FASTA file from clustering step",
    )
    p_func.add_argument(
        "-q", "--query", default=str(QUERY_FAA), metavar="FAA",
        help=f"Query FASTA file (default: {QUERY_FAA.name})",
    )
    p_func.add_argument(
        "--pfam", default=str(PFAM_DB), metavar="HMM",
        help=f"Path to Pfam-A.hmm (must be hmmpress-ed) "
             f"(default: {PFAM_DB})",
    )
    p_func.add_argument(
        "--dom-outdir", default=str(DOMAIN_DIR), metavar="DIR",
        help=f"Output directory for domain results (default: {DOMAIN_DIR.name}/)",
    )
    _add_thread_arg(p_func)

    # ── domain-filter ──
    p_dfilt = sub.add_parser(
        "domain-filter",
        help="Filter blast hits by overlap with a Pfam domain in query sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_dfilt.add_argument(
        "--domains", required=True, metavar="NAME[,NAME,...]",
        help="Comma-separated Pfam domain names or accessions to filter by "
             "(e.g. 'Peptidase_M64' or 'PF03009' or 'Peptidase_M26_N,Peptidase_M26_C')",
    )
    p_dfilt.add_argument(
        "--pfam-tsv", default=str(PFAM_RESULTS), metavar="TSV",
        help=f"pfam_results.tsv from 'domain' module "
             f"(default: {PFAM_RESULTS.relative_to(DATA_DIR)})",
    )
    p_dfilt.add_argument(
        "--blast", default=str(DEFAULT_BLAST_OUT), metavar="TSV",
        help=f"BLAST TSV to filter (full or identity-filtered) "
             f"(default: {DEFAULT_BLAST_OUT.name})",
    )
    p_dfilt.add_argument(
        "--out", default=str(DOM_FILTER_TSV), metavar="TSV",
        help=f"Output filtered hits TSV "
             f"(default: {DOM_FILTER_TSV.relative_to(DATA_DIR)})",
    )

    # ── plot ──
    p_plot = sub.add_parser(
        "plot",
        help="Generate figures: BLAST coverage, domain architecture, Foldseek scatter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_plot.add_argument(
        "--blast", default=None, metavar="TSV",
        help="Identity-filtered BLAST TSV (e.g. query_vs_pan_id20.tsv) → alignment coverage plot",
    )
    p_plot.add_argument(
        "--domain", default=None, metavar="DIR",
        help="Domain output directory (pfam_results.tsv, function_domain_filter.tsv, "
             "query.domtblout) → gggenes plot; requires --subject",
    )
    p_plot.add_argument(
        "--foldseek", default=None, metavar="TSV",
        help="Integrated results TSV (from 'integrate' module) → scatter plot",
    )
    p_plot.add_argument(
        "-q", "--query", default=str(QUERY_FAA), metavar="FAA",
        help=f"Query FASTA — needed for --blast and --domain plots "
             f"(default: {QUERY_FAA.name})",
    )
    p_plot.add_argument(
        "--subject", default=None, metavar="FAA",
        help="Subject FASTA — required when using --domain",
    )
    p_plot.add_argument(
        "--blast-threshold", type=float, default=30.0, metavar="FLOAT",
        help="BLAST identity threshold line in Foldseek scatter (default: 30.0)",
    )
    p_plot.add_argument(
        "--align-threshold", type=float, default=0.5, metavar="FLOAT",
        help="Foldseek alntmscore threshold for star highlight (default: 0.5)",
    )
    p_plot.add_argument(
        "--lddt-threshold", type=float, default=0.5, metavar="FLOAT",
        help="Foldseek lDDT threshold for star highlight (default: 0.5)",
    )
    p_plot.add_argument(
        "--all", dest="all_plots", action="store_true",
        help="Run all three plots using default input paths",
    )
    p_plot.add_argument(
        "--out", default=str(OUT_DIR), metavar="DIR",
        help=f"Output directory for all plots (default: {OUT_DIR.name}/)",
    )

    # ── integrate ──
    p_int = sub.add_parser(
        "integrate",
        help="Merge BLAST + Foldseek (+ optional domain/pLDDT) into one TSV",
    )
    p_int.add_argument(
        "--blast", required=True, metavar="TSV",
        help="BLAST results TSV (filtered or full)",
    )
    p_int.add_argument(
        "--foldseek", required=True, metavar="TSV",
        help="Foldseek alignment results TSV (from structure-align)",
    )
    p_int.add_argument(
        "--output", default=str(INTEGRATED_TSV), metavar="TSV",
        help=f"Output integrated TSV (default: {INTEGRATED_TSV.relative_to(DATA_DIR)})",
    )

    # ── structure-msa ──
    p_msa = sub.add_parser(
        "structure-msa",
        help="ColabFold: generate MSAs using EBI API (run on login node)",
    )
    p_msa.add_argument(
        "--fasta", required=True, metavar="FASTA",
        help="Input FASTA file with query sequences",
    )
    p_msa.add_argument(
        "--output", required=True, metavar="DIR",
        help="Output directory for MSA alignments",
    )
    _add_thread_arg(p_msa)

    # ── structure-predict ──
    p_pred = sub.add_parser(
        "structure-predict",
        help="ColabFold: submit structure prediction job to GPU cluster",
    )
    p_pred.add_argument(
        "--input", required=True, metavar="DIR",
        help="Input directory with MSA alignments",
    )
    p_pred.add_argument(
        "--output", required=True, metavar="DIR",
        help="Output directory for predicted structures",
    )
    p_pred.add_argument(
        "--weights", required=True, metavar="DIR",
        help="Path to ColabFold weights database",
    )
    p_pred.add_argument(
        "--gpus", type=int, default=1, metavar="N",
        help="Number of GPUs (default: 1)",
    )
    p_pred.add_argument(
        "--cpus", type=int, default=8, metavar="N",
        help="Number of CPUs per task (default: 8)",
    )
    p_pred.add_argument(
        "--partition", default="NV_4090D", metavar="PARTITION",
        help="SLURM partition (default: NV_4090D)",
    )
    p_pred.add_argument(
        "--job-name", default="colabfold", metavar="NAME",
        help="SLURM job name (default: colabfold)",
    )
    p_pred.add_argument(
        "--array-jobs", type=int, metavar="N",
        help="Number of array jobs (0 to N-1). Leave empty for single job",
    )
    p_pred.add_argument(
        "--array-concurrent", type=int, default=0, metavar="N",
        help="Max concurrent array jobs (0 = unlimited, default: 0)",
    )
    _add_thread_arg(p_pred)

    # ── structure-align ──
    p_align = sub.add_parser(
        "structure-align",
        help="Foldseek: structural alignment of predicted structures",
    )
    p_align.add_argument(
        "--query", metavar="PDB",
        help="Query PDB file (single structure)",
    )
    p_align.add_argument(
        "--query-dir", metavar="DIR",
        help="Query directory with multiple PDB files",
    )
    p_align.add_argument(
        "--targets", required=True, metavar="DIR",
        help="Target directory with PDB files to align against",
    )
    p_align.add_argument(
        "--output", required=True, metavar="TSV",
        help="Output TSV file with alignment results",
    )
    p_align.add_argument(
        "--threshold", type=float, default=0.5, metavar="FLOAT",
        help="TM-score threshold for filtering (default: 0.5)",
    )
    _add_thread_arg(p_align)

    # ── run-all ──
    p_all = sub.add_parser(
        "run-all",
        help="Run clustering → sequence → domain → structure --prep",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _add_clustering_args(p_all)
    _add_sequence_args(p_all)
    p_all.add_argument(
        "--pfam", default=str(PFAM_DB), metavar="HMM",
        help=f"Path to Pfam-A.hmm (default: {PFAM_DB})",
    )
    p_all.add_argument(
        "--dom-outdir", default=str(DOMAIN_DIR), metavar="DIR",
        help=f"Output directory for domain results (default: {DOMAIN_DIR.name}/)",
    )
    _add_thread_arg(p_all)

    # ── dispatch ──
    args = parser.parse_args()

    if args.module == "clustering":
        run_clustering(
            input_faa = Path(args.input),
            outdir    = Path(args.clust_output),
            tool      = args.tool,
            identity  = args.clust_identity,
            qcov      = args.clust_qcov,
            scov      = args.clust_scov,
            threads   = args.threads,
        )

    elif args.module == "sequence":
        db = Path(args.db) if args.db else DEFAULT_DB
        run_sequence(
            query_faa  = Path(args.query),
            diamond_db = db,
            blast_out  = Path(args.seq_output),
            identity   = args.seq_identity,
            qcov       = args.seq_qcov,
            scov       = args.seq_scov,
            threads    = args.threads,
            max_hits   = args.max_hits,
        )

    elif args.module == "domain":
        run_domain(
            subject_faa = Path(args.subject),
            threads     = args.threads,
            pfam_db     = Path(args.pfam),
            dom_outdir  = Path(args.dom_outdir),
            query_faa   = Path(args.query),
        )

    elif args.module == "domain-filter":
        domains = [d.strip() for d in args.domains.split(",") if d.strip()]
        run_domain_filter(
            domains   = domains,
            pfam_tsv  = Path(args.pfam_tsv),
            blast_tsv = Path(args.blast),
            out_tsv   = Path(args.out),
        )

    elif args.module == "plot":
        _blast_default = str(_filtered_blast(DEFAULT_BLAST_OUT, DEFAULT_SEQ_IDENTITY))
        run_plot(
            blast           = args.blast    or (_blast_default  if args.all_plots else None),
            domain          = args.domain   or (str(DOMAIN_DIR) if args.all_plots else None),
            foldseek        = args.foldseek or (str(INTEGRATED_TSV) if args.all_plots else None),
            subject_faa     = args.subject,
            query_faa       = Path(args.query),
            out_dir         = Path(args.out),
            blast_threshold = args.blast_threshold,
            align_threshold = args.align_threshold,
            lddt_threshold  = args.lddt_threshold,
        )

    elif args.module == "integrate":
        run_integrate(
            blast_tsv    = Path(args.blast),
            foldseek_tsv = Path(args.foldseek),
            output       = Path(args.output),
        )

    elif args.module == "structure-msa":
        run_structure_prep(
            fasta  = Path(args.fasta),
            output = Path(args.output),
        )

    elif args.module == "structure-predict":
        run_structure_predict(
            msa_dir         = Path(args.input),
            output          = Path(args.output),
            weights         = Path(args.weights),
            gpus            = args.gpus,
            cpus            = args.cpus,
            partition       = args.partition,
            job_name        = args.job_name,
            array_jobs      = args.array_jobs,
            array_concurrent = args.array_concurrent,
        )

    elif args.module == "structure-align":
        run_structure_align(
            query     = Path(args.query) if args.query else None,
            query_dir = Path(args.query_dir) if args.query_dir else None,
            targets   = Path(args.targets),
            output    = Path(args.output),
            threshold = args.threshold,
            threads   = args.threads,
        )

    elif args.module == "run-all":
        run_all(args)


if __name__ == "__main__":
    main()
