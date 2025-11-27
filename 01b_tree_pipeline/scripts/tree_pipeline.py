#!/usr/bin/env python3
"""
tree_pipeline_with_tnt_script.py

Pipeline:
  1) run IQ-TREE2 (model selection MFP + ML + SH-aLRT)
  2) write a TNT-friendly NEXUS for TNT
  3) create a TNT run-file (embedded in this script, modeled on your original)
  4) run TNT (non-interactively) using the generated run-file
  5) locate TNT output tree (auto-detected or provided)
  6) merge TNT branch lengths into IQ-TREE tree (SH-aLRT supports)
  7) write final combined tree

Usage example:
  python3 tree_pipeline_with_tnt_script.py \
    --alignment ../data/annotated_denv_genomes_nm.fasta \
    --outdir ./results \
    --prefix denv_tree \
    --iqtree-bin iqtree2 \
    --tnt-bin tnt \
    --tnt-max-ram 6000

Notes:
  - IQ-TREE2 must be in PATH (or specify the path).
  - TNT must be in PATH (or specify the path).
  - This script generates a TNT run-file (tnt_script_<prefix>.run) in the working directory and runs TNT < tnt_script.run.
"""
import argparse
import subprocess
import time
from pathlib import Path
from Bio import SeqIO, Phylo

def run(cmd, cwd=None, env=None):
    print(f"RUN: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=cwd, env=env)

def write_tnt_nexus_from_fasta(fasta_path: Path, nexus_out: Path):
    seqs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not seqs:
        raise SystemExit("No sequences found in alignment.")
    seq_len = len(seqs[0].seq)
    num_taxa = len(seqs)
    with nexus_out.open("w") as f:
        f.write("xread\n")
        f.write(f"{seq_len} {num_taxa}\n")
        for rec in seqs:
            # Avoid spaces in taxon names — user should ensure names match across tools
            f.write(f"{rec.id} {str(rec.seq)}\n")
        f.write(";\n")
    print(f"Wrote TNT NEXUS: {nexus_out}")

def find_recent_tree_candidate(directory: Path, since_ts: float = 0.0):
    candidates = []
    for ext in (".nex", ".nwk", ".tre", ".treefile", ".tree"):
        for p in directory.rglob(f"*{ext}"):
            try:
                mtime = p.stat().st_mtime
            except Exception:
                continue
            if mtime >= since_ts:
                candidates.append((mtime, p))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]

def leafset_map(tree):
    mapping = {}
    for clade in tree.find_clades(order="preorder"):
        leaves = [t.name for t in clade.get_terminals()]
        if leaves:
            mapping[frozenset(leaves)] = clade
    return mapping

def merge_branch_lengths(branch_tree, support_tree):
    src_map = leafset_map(branch_tree)
    dst_map = leafset_map(support_tree)
    copied = 0
    for leafset, dst_clade in dst_map.items():
        src_clade = src_map.get(leafset)
        if src_clade is not None:
            dst_clade.branch_length = src_clade.branch_length
            copied += 1
    print(f"Copied branch lengths for {copied} clades (matched by leaf sets).")
    return support_tree

def create_tnt_input_tree(iqtree_treefile: Path, tnt_tree_path: Path):
    """
    Save the IQ-TREE topology into a file that TNT can read.
    We write the Newick tree verbatim to tnt_tree_path.
    """
    text = iqtree_treefile.read_text()
    # write as-is; user may adjust if TNT requires special formatting
    tnt_tree_path.write_text(text)
    print(f"Wrote TNT input tree: {tnt_tree_path}")

def write_tnt_runfile(outdir: Path, prefix: str, tnt_max_ram: int):
    """
    Compose a TNT runfile similar to your original script.
    This file will be written to outdir/tnt_script_<prefix>.run
    """
    runfile_path = outdir / f"tnt_script_{prefix}.run"
    # The runfile uses commands inspired by your original.
    # Note: TNT syntax must be compatible with your local TNT version.
    content = f"""mxram {tnt_max_ram};
nstates dna;
proc {prefix}.nex;
proc TNT_input_tree_{prefix}.nwk;
# optimize branch lengths (parsimony-based)
blength *;
# export the tree with branch lengths
export > TNT_original_output_tree_{prefix}.nex;
quit;
"""
    runfile_path.write_text(content)
    print(f"Wrote TNT run-file: {runfile_path}")
    return runfile_path

def main():
    p = argparse.ArgumentParser(description="IQ-TREE2 -> TNT (embedded runfile) -> merge")
    p.add_argument("--alignment", "-s", required=True, type=Path, help="Input MSA FASTA")
    p.add_argument("--outdir", "-od", type=Path, default=Path("."), help="Output directory")
    p.add_argument("--prefix", "-p", type=str, default="tree", help="Output prefix")
    p.add_argument("--iqtree-bin", type=str, default="iqtree2", help="IQ-TREE binary (default: iqtree2)")
    p.add_argument("--tnt-bin", type=str, default="tnt", help="TNT binary (default: tnt)")
    p.add_argument("--tnt-max-ram", type=int, default=6000, help="mxram setting for TNT (MB)")
    p.add_argument("--tnt-output", type=Path, default=None, help="Optional: explicit TNT output tree file")
    p.add_argument("--threads", "-nt", type=int, default=0, help="Threads for IQ-TREE (0 = AUTO)")
    args = p.parse_args()

    alignment = args.alignment.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix
    iqtree = args.iqtree_bin
    tnt = args.tnt_bin
    tnt_max_ram = args.tnt_max_ram
    tnt_output = args.tnt_output.resolve() if args.tnt_output else None
    threads = args.threads

    if not alignment.exists():
        raise SystemExit(f"Alignment not found: {alignment}")

    # 1) Run IQ-TREE2 with model selection (MFP) + SH-aLRT support
    iq_pre = outdir / f"{prefix}_iqtree"
    m_arg = "MFP"
    nt_arg = "-nt AUTO" if threads == 0 else f"-nt {threads}"
    # -redo ensures reruns overwrite previous files
    iq_cmd = f"{iqtree} -s {alignment} -m {m_arg} -alrt 1000 {nt_arg} -pre {iq_pre} -redo"
    run(iq_cmd, cwd=outdir)

    # detect IQ-TREE output treefile
    iq_treefile = iq_pre.with_suffix(".treefile")
    if not iq_treefile.exists():
        cand = iq_pre.with_suffix(".bestTree")
        if cand.exists():
            iq_treefile = cand
        else:
            raise SystemExit("IQ-TREE did not produce expected treefile.")
    print(f"IQ-TREE tree-with-supports: {iq_treefile}")

    # 2) produce TNT-friendly NEXUS
    tnt_nexus = outdir / f"{prefix}.nex"
    write_tnt_nexus_from_fasta(alignment, tnt_nexus)

    # 3) create TNT input tree file (from IQ-TREE topology)
    tnt_input_tree = outdir / f"TNT_input_tree_{prefix}.nwk"
    create_tnt_input_tree(iq_treefile, tnt_input_tree)

    # 4) write the TNT run-file (embedded)
    tnt_runfile = write_tnt_runfile(outdir, prefix, tnt_max_ram)

    # 5) run TNT non-interactively using generated run-file
    print("Running TNT (this will execute the generated run-file).")
    tnt_start = time.time()
    # run TNT by redirecting the runfile into TNT; outputs will be created in outdir if TNT writes relative paths
    tnt_cmd = f"{tnt} < {tnt_runfile}"
    run(tnt_cmd, cwd=outdir)
    tnt_end = time.time()
    print("TNT run completed.")

    # 6) locate TNT output tree
    if tnt_output:
        tnt_tree_path = tnt_output
        if not tnt_tree_path.exists():
            raise SystemExit(f"Provided --tnt-output does not exist: {tnt_tree_path}")
    else:
        tnt_tree_path = find_recent_tree_candidate(outdir, since_ts=tnt_start - 1.0)
        if tnt_tree_path is None:
            raise SystemExit("Could not find TNT output tree automatically. Use --tnt-output to specify it.")
    print(f"Using TNT output tree: {tnt_tree_path}")

    # 7) Read trees (TNT branch-length tree and IQ-TREE support tree)
    print("Reading trees...")
    try:
        branch_tree = Phylo.read(str(tnt_tree_path), format="newick")
    except Exception:
        branch_tree = Phylo.read(str(tnt_tree_path), format="nexus")

    support_tree = Phylo.read(str(iq_treefile), format="newick")

    # 8) Merge branch lengths into the IQ-TREE support tree
    merged = merge_branch_lengths(branch_tree, support_tree)

    final_out = outdir / f"final_tree_{prefix}.nwk"
    Phylo.write(merged, str(final_out), "newick")
    print(f"Final merged tree written to: {final_out}")

if __name__ == "__main__":
    main()
