#!/bin/bash
#BSUB -J separate_cnv_table
#BSUB -o separate_cnv_table.out
#BSUB -e separate_cnv_table.err
#BSUB -n 10
#BSUB -M 10000
#BSUB -R "select[mem>10000] rusage[mem=10000]"
#BSUB -q "normal"

# Define suffix
suffix="v1.0_"
# Define path to mapper file
MAPPER_FILE="/nfs/users/nfs_j/jn11/my_lustre/analysis/sample_set_mapper.csv"
# Define input file name
INPUT_FILE="/nfs/users/nfs_j/jn11/my_lustre/analysis/v1.0_afun/coverage_archive/full_raw_CNV_table.csv"
# Define output directory
DESTINATION_DIR="/nfs/users/nfs_j/jn11/my_lustre/analysis" #/$suffix$GROUP
# Define working directory
WORKING_DIR="/nfs/users/nfs_j/jn11/my_lustre/analysis/v1.0_afun/coverage_archive"

# Store the header of the INPUT_FILE
HEADER=$(head -n 1 "$INPUT_FILE")

# Create a temporary file
TMP_FILE=$(mktemp)

# Create a CSV file for each group in the working directory
echo "Creating directories for groups:"
cut -d ',' -f1 "$MAPPER_FILE" | tail -n +2 | sort | uniq | while read -r GROUP; do
  # Add the header to each CSV file
  echo "$HEADER" > "$WORKING_DIR/$GROUP.csv"
done

# Process the input file and distribute lines to appropriate CSV files
echo "Processing the input file:"
while IFS=$'\t' read -r LINE; do
  # Extract the ID from the first column
  ID=$(echo "$LINE" | cut -d ',' -f1 | cut -d ':' -f1)

  # Look for a matching group in the mapper file
  GROUP=$(grep ",$ID$" "$MAPPER_FILE" | cut -d ',' -f1)

  # If a matching group was found, append the line to the corresponding CSV file
  if [ -n "$GROUP" ]; then
    echo "${LINE%$'\t'*}" >> "$WORKING_DIR/$GROUP.csv"
  else
    # If no matching group was found, append the line to all CSV files (including those without a group)
    echo "${LINE%$'\t'*}" >> "$WORKING_DIR/.csv"
  fi
done < "$INPUT_FILE" > "$TMP_FILE"

# Move the CSV files to the appropriate directory in DESTINATION_DIR
echo "Moving CSV files to the destination directory:"
cut -d ',' -f1 "$MAPPER_FILE" | sort | uniq | while read -r GROUP; do
  mkdir -p "$DESTINATION_DIR/$suffix$GROUP/coverage_CNV"
  if [ -f "$WORKING_DIR/$GROUP.csv" ]; then
    mv "$WORKING_DIR/$GROUP.csv" "$DESTINATION_DIR/$suffix$GROUP/coverage_CNV/full_raw_CNV_table.csv"
  fi
done

# Cleanup temporary file
rm "$TMP_FILE"
