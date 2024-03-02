#!/bin/bash

# This script runs the CAT (Comparative Annotation Toolkit) pipeline for fly genome annotation.

# Check if the required command line arguments are provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 <clade_file> <hints_file>"
  exit 1
fi

# Assign command line arguments to variables
clade_file=$1
hints_file=$2

# Extract necessary information from the clade file
root_anc=$(awk 'NR==3' "$clade_file" | sed 's/#Root: //')
hal_name="${clade_file%_tree.nw}.hal"
refsp=$(awk -F " " '/#Reference_species:/ {print $2}' "$clade_file")

# Create the working directory
work_dir="${clade_file%_tree.nw}_${hints_file}/${clade_file%_tree.nw}_work_dir"
mkdir -p "$work_dir"

# Run the CAT pipeline using specified parameters
luigi --module cat RunCat \
  --hal="CAT_run/$hal_name" \
  --ref-genome="$refsp" \
  --config="CAT_run/${clade_file%_tree.nw}.cat.config" \
  --work-dir "$work_dir" \
  --out-dir "${clade_file%_tree.nw}_${hints_file}/${clade_file%_tree.nw}_out_dir" \
  --augustus \
  --local-scheduler \
  --augustus-species fly \
  --augustus-cgp \
  --cgp-param augustus_cfgs/cgp_parameters.cfg \
  --augustus-cgp-cfg-template augustus_cfgs/cgp_extrinsic_template.cfg \
  --tm-cfg augustus_cfgs/extrinsic.ETM1.cfg \
  --tmr-cfg augustus_cfgs/extrinsic.ETM2.cfg \
  --assembly-hub \
  --binary-mode local \
  --workers=40 \
  --maxCores=8


