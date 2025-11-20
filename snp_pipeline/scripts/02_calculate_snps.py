import os
import pandas as pd
from Bio import SeqIO

# -------------------------------------------------------------------
# Correct paths based on your directory tree
# -------------------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
FRAG_DIR = os.path.join(DATA_DIR, "fragments")

input_file = os.path.join(DATA_DIR, "recomb_and_parents.csv")
output_file = os.path.join(DATA_DIR, "recombinant_snps.csv")

genome_size = 15213

# -------------------------------------------------------------------
# Function to compute ungapped length
# -------------------------------------------------------------------
def get_ungapped_length(fasta_file):
    if not os.path.exists(fasta_file):
        print(f"Warning: FASTA file not found: {fasta_file}")
        return 0
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return len(str(record.seq).replace('-', '').replace('.', ''))
    return 0

# -------------------------------------------------------------------
# Load *input* file but DO NOT overwrite it
# -------------------------------------------------------------------
df = pd.read_csv(input_file)

# -------------------------------------------------------------------
# Add output columns to the new dataframe
# -------------------------------------------------------------------
required_cols = [
    "SNPs (in recombinant region)",
    "SNPs (in non-recombinant region)",
    "Recombinant Length (bp)",
    "Non-Recombinant Length (bp)"
]

for col in required_cols:
    if col not in df.columns:
        df[col] = "" if "SNPs" in col else 0

# -------------------------------------------------------------------
# Main processing
# -------------------------------------------------------------------
for idx, row in df.iterrows():
    recombinant = str(row["Recombinant"]).strip()
    minor_parent = str(row["Minor parent"]).strip() if pd.notna(row["Minor parent"]) else "NA"
    major_parent = str(row["Major parent"]).strip() if pd.notna(row["Major parent"]) else "NA"

    base_name = f"snipit_{recombinant}_{minor_parent}_{major_parent}"

    # -----------------------------------------------------
    # Recombinant region (frag1)
    # -----------------------------------------------------
    minor_snps_total = 0
    recombinant_length_total = 0

    all_frag1_files = [
        f for f in os.listdir(FRAG_DIR)
        if f.startswith(base_name) and "_frag1_" in f and f.endswith(".csv")
    ]

    for frag1_csv in all_frag1_files:
        frag1_csv_path = os.path.join(FRAG_DIR, frag1_csv)
        frag1_fasta_path = frag1_csv_path.replace(".csv", ".fasta")

        recombinant_length_total += get_ungapped_length(frag1_fasta_path)

        frag1_df = pd.read_csv(frag1_csv_path)
        if minor_parent != "NA":
            row_minor = frag1_df[frag1_df["record"] == minor_parent]
            if not row_minor.empty:
                minor_snps_total += int(row_minor["num_snps"].iloc[0])

    rec_percentage = (minor_snps_total / recombinant_length_total) * 100 if recombinant_length_total > 0 else 0

    # -----------------------------------------------------
    # Non-recombinant region (frag2)
    # -----------------------------------------------------
    major_snps_total = 0
    non_recombinant_length_total = 0

    all_frag2_files = [
        f for f in os.listdir(FRAG_DIR)
        if f.startswith(base_name) and "_frag2_outside_" in f and f.endswith(".csv")
    ]

    for frag2_csv in all_frag2_files:
        frag2_csv_path = os.path.join(FRAG_DIR, frag2_csv)
        frag2_fasta_path = frag2_csv_path.replace(".csv", ".fasta")

        non_recombinant_length_total += get_ungapped_length(frag2_fasta_path)

        frag2_df = pd.read_csv(frag2_csv_path)
        if major_parent != "NA":
            row_major = frag2_df[frag2_df["record"] == major_parent]
            if not row_major.empty:
                major_snps_total += int(row_major["num_snps"].iloc[0])

    non_rec_percentage = (major_snps_total / non_recombinant_length_total) * 100 if non_recombinant_length_total > 0 else 0

    # Store results in the NEW OUTPUT TABLE
    df.at[idx, "SNPs (in recombinant region)"] = f"{minor_snps_total} ({rec_percentage:.2f}%)"
    df.at[idx, "SNPs (in non-recombinant region)"] = f"{major_snps_total} ({non_rec_percentage:.2f}%)"
    df.at[idx, "Recombinant Length (bp)"] = recombinant_length_total
    df.at[idx, "Non-Recombinant Length (bp)"] = non_recombinant_length_total

# -------------------------------------------------------------------
# SAVE NEW TABLE — DO NOT MODIFY THE ORIGINAL INPUT FILE
# -------------------------------------------------------------------
df.to_csv(output_file, index=False)

print(f"✔ SNP results written to: {output_file}")

