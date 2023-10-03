#!/bin/bash

# Read command line arguments
reference="$1"
long_reads="$2"
output_prefix="$3"
flag="${4:-4}"  # Set default value of flag to 4 if not provided

# Check if output directory already exists
if [ -d "$output_prefix" ]; then
    echo "Output directory '$output_prefix' already exists. Storing results in existing directory."
else
    # Create the output directory
    mkdir "$output_prefix"
    echo "Output directory '$output_prefix' created."
fi

# Copy reference genome to output directory
cp "$reference" "$output_prefix/"

# Run minimap2 to perform mapping
minimap2 -ax map-ont "$reference" "$long_reads" > "$output_prefix/output.sam"

# Convert SAM file to BAM file and filter unmapped reads
echo $flag
samtools view -b -F "$flag" -o "$output_prefix/output.bam" "$output_prefix/output.sam"

# Sort the BAM file
samtools sort -o "$output_prefix/sorted.bam" "$output_prefix/output.bam"

# Generate BAM index
samtools index "$output_prefix/sorted.bam"

# Calculate coverage
samtools depth "$output_prefix/sorted.bam" > "$output_prefix/coverage.txt"

# Remove intermediate result files
rm "$output_prefix/output.sam"
rm "$output_prefix/output.bam"