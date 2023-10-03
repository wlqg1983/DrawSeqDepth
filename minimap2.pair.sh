#!/bin/bash

# Read command line arguments
reference="$1"
reads1="$2"
reads2="$3"
output_prefix="$4"
flag="${5:-4}"  # Set default value of flag to 4 if not provided

# Append "_results" to the output directory name
output_dir="$output_prefix"

# Check if output directory already exists
if [ -d "$output_dir" ]; then
    echo "Output directory '$output_dir' already exists. Storing results in existing directory."
else
    # Create the output directory
    mkdir "$output_dir"
    echo "Output directory '$output_dir' created."
fi

# Copy reference genome to output directory
cp "$reference" "$output_dir/"

# Run minimap2 to perform mapping
minimap2 -ax sr "$reference" "$reads1" "$reads2" > "$output_dir/output.sam"

# Convert SAM file to BAM file and filter unmapped reads
echo $flag
samtools view -b -F $flag -o "$output_dir/output.bam" "$output_dir/output.sam"

# Sort the BAM file
samtools sort -o "$output_dir/sorted.bam" "$output_dir/output.bam"

# Generate BAM index
samtools index "$output_dir/sorted.bam"

# Calculate coverage
samtools depth "$output_dir/sorted.bam" > "$output_dir/coverage.txt"

# Remove intermediate result files
rm "$output_dir/output.sam"
rm "$output_dir/output.bam"