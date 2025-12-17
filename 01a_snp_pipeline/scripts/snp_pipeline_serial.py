#!/usr/bin/env python3
import subprocess
import os

BASE = "/home/hugo/hpc_flavirecomb"

def run(cmd):
    print(f"\n=== Running: {cmd} ===")
    subprocess.run(cmd, shell=True, check=True)

def main():
    print("=== 00 SNP Pipeline ===")

    run(f"python3 {BASE}/01a_snp_pipeline/scripts/00_auto_snipit.py")
    run(f"bash {BASE}/01a_snp_pipeline/scripts/01_run_snipit.sh")
    run(f"python3 {BASE}/01a_snp_pipeline/scripts/02_calculate_snps.py")

    print("\n=== 00 SNP Pipeline Completed ===")


if __name__ == "__main__":
    main()

