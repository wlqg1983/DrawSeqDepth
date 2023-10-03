#!/bin/bash

# Check if all required arguments are provided
if [[ $# -lt 4 ]]; then
    echo "Missing required arguments."
    echo "Usage: bash hisat2.paired.sh reference.fasta input1.fastq input2.fastq output_dir [flag]"
    exit 1
fi

# Parse command line arguments
reference="$1"
input_fastq1="$2"
input_fastq2="$3"
output_dir="$4"
flag="${5:-4}" # Set default value of flag to 4 if not provided

# Create output directory if it does not exist
mkdir -p "$output_dir"

# Change working directory to the output directory
cd "$output_dir"

# Copy reference genome to output directory
cp "../$reference" .

# Get the filename of the reference genome
reference_name=$(basename "$reference")

# Index reference genome
hisat2-build "$reference_name" "$reference_name"

# Align paired-end reads to the reference genome and save intermediate result in the current working directory
hisat2 -x "$reference_name" -1 "../$input_fastq1" -2 "../$input_fastq2" -S "output.sam"

# Convert SAM file to BAM file and filter unmapped reads using specified flag (or default value)
samtools view -b -F "$flag" -o "output.bam" "output.sam"

# Sort the BAM file
samtools sort -o "sorted.bam" "output.bam"

# Generate BAM index
samtools index "sorted.bam"

# Calculate coverage
samtools depth "sorted.bam" > "coverage.txt"

# Remove intermediate result files in the current working directory
rm -f "$reference_name".*".ht2" output.sam output.bam

echo "Processing complete."