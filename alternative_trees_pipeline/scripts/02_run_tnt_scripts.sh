#!/bin/bash
TNT_SCRIPTS_DIR="/home/hugo/hpc_flavirecomb/alternative_trees_pipeline/results"

TNT_EXEC="/home/hugo/Downloads/TNT-bin/tnt"

for tnt_script in "$TNT_SCRIPTS_DIR"/script_*.RUN; do
    echo "Running TNT script: $tnt_script"
    "$TNT_EXEC" < "$tnt_script"
done

echo "All TNT runs completed."
