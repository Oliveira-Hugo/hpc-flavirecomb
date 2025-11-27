#!/usr/bin/env python3
import os
import subprocess

BASE = "/home/hugo/hpc_flavirecomb"

FINAL_TREE = f"{BASE}/01b_tree_pipeline/scripts/topology_final.nwk"
ALT_TREES_DIR = f"{BASE}/01c_alternative_trees_pipeline/scripts"
OUTPUT_DIR = f"{BASE}/02_comp_trees_pipeline/results"
YBYRA = "/home/hugo/ybyra/ybyra_sa.py"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def run(cmd):
    print(f"\n=== Running: {cmd} ===")
    subprocess.run(cmd, shell=True, check=True)

def main():
    print("=== 03 Compare Trees Pipeline ===")

    consensus_files = [
        f for f in os.listdir(ALT_TREES_DIR)
        if f.startswith("consensus_") and f.endswith(".tre")
    ]

    if not consensus_files:
        print("ERROR: No consensus*.tre files found in alternative trees directory.")
        return

    config_files = []

    for idx, treefile in enumerate(consensus_files, start=1):

        alt_tree_path = os.path.join(ALT_TREES_DIR, treefile)
        config_path = os.path.join(OUTPUT_DIR, f"config_{idx}.txt")

        config_content = f""">id = calculating_topological_distances
<begin files
    {FINAL_TREE} ;
    {alt_tree_path} ;
end files>

>n = 1 {FINAL_TREE} ] 
>opt = 3
>compare = 0
>verbose
"""

        with open(config_path, "w") as cfg:
            cfg.write(config_content)

        config_files.append(config_path)

    for cfg in config_files:
        print(f"\nRunning comparison for: {cfg}")
        run(f"python3 {YBYRA} -d -f {cfg}")

    print("\n=== 03 Compare Trees Pipeline Completed ===")

if __name__ == "__main__":
    main()

