#!/usr/bin/env python3
import os
import re

def convert_tnt_to_newick(content):
    """Apply TNT → Newick conversion rules."""

    content = content.replace('*', ';')
    content = content.replace(' ', ',')
    content = content.replace(')', '),')
    content = re.sub(r',\)', ')', content)
    content = re.sub(r',;', ';', content)

    return content.strip()


def process_tnt_file(filename):
    """Process a single TNT file."""
    with open(filename, "r") as infile:
        lines = infile.readlines()

    if len(lines) < 3:
        print(f"⚠️ File too short, skipping: {filename}")
        return

    inner = ''.join(lines[1:-1])
    newick = convert_tnt_to_newick(inner)
    base = os.path.splitext(filename)[0]
    clean_base = re.sub(r"^consensus_", "", base)
    output_filename = f"consensus_{clean_base}.tre"

    with open(output_filename, "w") as outfile:
        outfile.write(newick + "\n")

    print(f"Converted: {filename} → {output_filename}")


def main():
    print("=== Converting TNT Trees to Newick Format ===")

    files = [f for f in os.listdir(".") if f.startswith("consensus") and f.endswith(".tnt")]

    if not files:
        print("No .tnt files found in this directory.")
        return

    for f in files:
        process_tnt_file(f)

    print("=== Conversion Completed ===")


if __name__ == "__main__":
    main()

