#!/bin/bash

mamba activate brocspin

MAIN_OUTPUT_DIR="/home/allenchen/brocspin/trajectories/temp_batch"
mkdir "$MAIN_OUTPUT_DIR"

RUN_SCRIPT="/home/allenchen/brocspin/script/run_batch.py"
PDB_FILE="/home/allenchen/brocspin/structures/CG_pdb/8k7w_CG.pdb"
PSF_FILE="/home/allenchen/brocspin/structures/psf/8k7w.psf"

SALT=150.0
TEMP=("270" "280" "290" "300")

for ((i=0; i<${#TEMP[@]}; i++)); do
  TEM="${TEMP[$i]}"

  OUTPUT_DIR="$MAIN_OUTPUT_DIR/b${TEM}"
  mkdir "$OUTPUT_DIR"

  rm -f system*

  echo "Simulating with TEMP=$TEM, nonperiodic with SALT=$SALT"
  python "$RUN_SCRIPT" \
    --pdb "$PDB_FILE" \
    --psf "$PSF_FILE" \
    --temp "$TEM" \
    --ens non \
    --salt "$SALT" \
    > "$OUTPUT_DIR/run.log" 2>&1 

  if [ $? -eq 0 ]; then
    mv system.pdb system.log system.chk system.xtc system.xml "$OUTPUT_DIR/" 2>/dev/null || true
    echo "Simulation for $PDB_FILE completed. Results are in $OUTPUT_DIR."
  else
    echo "Simulation for $PDB_FILE failed. Check $OUTPUT_DIR/run.log for details."
  fi
done

echo "All simulations completed."
