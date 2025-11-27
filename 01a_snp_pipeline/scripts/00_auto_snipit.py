#!/usr/bin/env python3
from Bio import SeqIO
import csv
import os
from pathlib import Path

root = Path(__file__).resolve().parent.parent
script_dir = Path(__file__).resolve().parent
project_root = script_dir.parent.parent

fasta_file = project_root / "00_input" / "annotated_denv_genomes_nm.fasta"
table_file = project_root / "00_input" / "recomb_and_parents.csv"
output_dir = root / "results" / "fragments"
alignment_length = 15213  # fixed alignment length

# ==== CREATE OUTPUT DIR ====
os.makedirs(output_dir, exist_ok=True)

# ==== LOAD FASTA ====
seq_dict = {rec.id: rec for rec in SeqIO.parse(fasta_file, "fasta")}

# ==== READ TABLE AND PROCESS ====
with open(table_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=",")
    for row in reader:
        recombinant = row["Recombinant"].strip()
        minor = row["Minor parent"].strip()
        major = row["Major parent"].strip()
        begin = int(row["Begin"])
        end = int(row["End"])

        taxa = [t for t in [recombinant, minor, major] if t != "Unknown"]

        if not taxa:
            continue

        missing = [t for t in taxa if t not in seq_dict]
        if missing:
            print(f"Warning: Missing sequences for {missing}, skipping line.")
            continue

        base_name = f"snipit_{recombinant}_{minor}_{major}".replace("Unknown", "NA")
        out_fasta_path = os.path.join(output_dir, f"{base_name}.fasta")

        SeqIO.write([seq_dict[t] for t in taxa], out_fasta_path, "fasta")

        frag1_records = []
        for t in taxa:
            seq = seq_dict[t].seq
            frag1_seq = seq[begin - 1:end]
            frag1_records.append(seq_dict[t][:0])
            frag1_records[-1].seq = frag1_seq

        frag1_path = os.path.join(output_dir, f"{base_name}_frag1_{begin}_{end}.fasta")
        SeqIO.write(frag1_records, frag1_path, "fasta")

        frag2_records = []
        for t in taxa:
            seq = seq_dict[t].seq
            frag2_seq = seq[0:begin - 1] + seq[end:alignment_length]
            frag2_records.append(seq_dict[t][:0])
            frag2_records[-1].seq = frag2_seq

        frag2_path = os.path.join(output_dir, f"{base_name}_frag2_outside_{begin}_{end}.fasta")
        SeqIO.write(frag2_records, frag2_path, "fasta")

print("Done.")

