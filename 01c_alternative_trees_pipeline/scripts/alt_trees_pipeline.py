#!/usr/bin/env python3
import subprocess

BASE = "/home/hugo/hpc_flavirecomb"

def run(cmd):
    print(f"\n=== Running: {cmd} ===")
    subprocess.run(cmd, shell=True, check=True)

def main():
    print("=== 02 Alternative Trees Pipeline ===")

    run(f"python3 {BASE}/01c_alternative_trees_pipeline/scripts/00_prepare_alt_alignments.py")
    run(f"python3 {BASE}/01c_alternative_trees_pipeline/scripts/01_prepare_tnt_scripts.py")
    run(f"bash {BASE}/01c_alternative_trees_pipeline/scripts/02_run_tnt_scripts.sh")
    run(f"python3 {BASE}/01c_alternative_trees_pipeline/scripts/03_convert_trees.py")
    
    print("\n=== 02 Alternative Trees Pipeline Completed ===")

if __name__ == "__main__":
    main()

