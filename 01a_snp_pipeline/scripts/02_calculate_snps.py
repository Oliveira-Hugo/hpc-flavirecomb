#!/usr/bin/env python3

import csv
from pathlib import Path

# =========================
# Paths
# =========================

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent.parent

INPUT_FILE = PIPELINE_ROOT / "00_input" / "recomb_and_parents.csv"
OUTPUT_FILE = PIPELINE_ROOT / "01a_snp_pipeline" / "results" / "recombinant_snps.csv"

SNIPIT_OUTPUTS_DIR = Path("/scratch/cenapadrjsd/hugo.oliveira2/hpc-flavirecomb/01a_snp_pipeline/results/snipit_outputs")
FRAGMENTS_DIR = PIPELINE_ROOT / "01a_snp_pipeline" / "results" / "fragments"

# =========================
# Helpers
# =========================

def normalize_parent(value):
    v = (value or "").strip()
    if v == "" or v.lower() == "unknown":
        return "NA"
    return v

def read_csv_dict(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def read_snps(snps_csv, parent):
    if parent == "NA":
        return 0
    with open(snps_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("record") == parent:
                try:
                    return int(row.get("num_snps", 0))
                except ValueError:
                    return 0
    return 0

def fasta_length(fasta_path):
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.strip().replace("-", "").replace(".", ""))
    return length

# =========================
# Load input
# =========================

rows = read_csv_dict(INPUT_FILE)

for row in rows:
    row["SNPs (in recombinant region)"] = ""
    row["SNPs (in non-recombinant region)"] = ""
    row["Recombinant Length (bp)"] = "0"
    row["Non-Recombinant Length (bp)"] = "0"

# =========================
# Main loop
# =========================

for row in rows:
    recombinant = row["Recombinant"].strip()
    minor = normalize_parent(row["Minor parent"])
    major = normalize_parent(row["Major parent"])

    prefix = f"snipit_{recombinant}_{minor}_{major}"

    minor_snps = 0
    major_snps = 0
    rec_len = 0
    nonrec_len = 0

    for subdir in SNIPIT_OUTPUTS_DIR.iterdir():
        if not subdir.is_dir():
            continue
        name = subdir.name

        if not name.startswith(prefix):
            continue

        if "_frag1_" in name:
            snps_csv = subdir / "snps.csv"
            fasta = FRAGMENTS_DIR / f"{name}.fasta"

            if snps_csv.exists():
                minor_snps += read_snps(snps_csv, minor)
            if fasta.exists():
                rec_len += fasta_length(fasta)

        elif "_frag2_outside_" in name:
            snps_csv = subdir / "snps.csv"
            fasta = FRAGMENTS_DIR / f"{name}.fasta"

            if snps_csv.exists():
                major_snps += read_snps(snps_csv, major)
            if fasta.exists():
                nonrec_len += fasta_length(fasta)

    rec_pct = (minor_snps / rec_len * 100) if rec_len > 0 else 0
    nonrec_pct = (major_snps / nonrec_len * 100) if nonrec_len > 0 else 0

    row["SNPs (in recombinant region)"] = f"{minor_snps} ({rec_pct:.2f}%)"
    row["SNPs (in non-recombinant region)"] = f"{major_snps} ({nonrec_pct:.2f}%)"
    row["Recombinant Length (bp)"] = str(rec_len)
    row["Non-Recombinant Length (bp)"] = str(nonrec_len)

# =========================
# Write output
# =========================

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"SNP results written to: {OUTPUT_FILE}")
