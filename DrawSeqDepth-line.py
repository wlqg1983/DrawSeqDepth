#!/usr/bin/env python
import argparse
import os
import subprocess
import sys
import shutil
import matplotlib
import matplotlib.pyplot as plt
import time

matplotlib.use('Agg')


def draw_seq_depth(alignment, reference, sequencing_type, genome_type, output_format, output_prefix, input_fastq, flag):
    # Get the absolute path of the current working directory
    abs_work_dir = os.path.dirname(os.path.realpath(__file__))

    # Call the alignment script based on sequencing type and alignment type
    if sequencing_type == "single":
        if alignment == "minimap2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "minimap2.single.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "bowtie2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "bowtie2.single.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "bwa":
            subprocess.call(["bash", os.path.join(abs_work_dir, "bwa.single.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "hisat2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "hisat2.single.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
    elif sequencing_type == "pair":
        if alignment == "minimap2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "minimap2.pair.sh"), reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "bowtie2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "bowtie2.pair.sh"), reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "bwa":
            subprocess.call(["bash", os.path.join(abs_work_dir, "bwa.pair.sh"), reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "hisat2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "hisat2.pair.sh"), reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
    elif sequencing_type == "third":
        if alignment == "minimap2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "minimap2.third.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "bowtie2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "bowtie2.third.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "bwa":
            subprocess.call(["bash", os.path.join(abs_work_dir, "bwa.third.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
        elif alignment == "hisat2":
            subprocess.call(["bash", os.path.join(abs_work_dir, "hisat2.third.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), str(flag)])
    
    # Set the absolute paths for coverage_file and output_file using os.getcwd()
    coverage_file = os.path.join(os.getcwd(), output_prefix + "_results", "coverage.txt")
    print(coverage_file)
    time.sleep(10)
    output_file = os.path.join(os.getcwd(), output_prefix + "_results", output_prefix)

    # Loop over each output format
    for format in output_format:
        # Call plot.py script and pass command line arguments
        #subprocess.call([str(sys.executable), "plot.py", coverage_file, output_file, format])
        subprocess.call(['plot.py', coverage_file, output_file, format])
        
    # Get the current path
    current_path = os.path.join(os.getcwd(), output_prefix + "_results")

    # Get all files and directories in the current path
    files = os.listdir(current_path)

    # Iterate through and print the absolute paths of files and directories
    for file in files:
        file_path = os.path.join(current_path, file)
        print(file_path)



output_prefix = "default_output_results"

if __name__ == "__main__":
    # Create argparse object and set up command line arguments
    parser = argparse.ArgumentParser(prog="DrawSeqDepth-line.py", description="Draw sequence depth.")

    parser.add_argument("-alignment",
                        choices=["minimap2", "hisat2", "bowtie2", "bwa"],default="minimap2",
                        help="Specifies the alignment software to be used. Default is minimap2.")
    parser.add_argument("-reference", type=str,
                        help="Specifies the reference genome or a specific segment of the genome.")
    parser.add_argument("-single", metavar="FASTQ", help="Specifies single-end sequencing data.")
    parser.add_argument("-pair", metavar=("FASTQ1", "FASTQ2"), nargs=2, help="Specifies paired-end sequencing data.")
    parser.add_argument("-third", metavar="FASTQ", help="Specifies third-generation sequencing data.")
    parser.add_argument("-format", default=["pdf"], nargs='+', choices=["png", "pdf", "jpeg", "bmp", "tiff", "gif", "eps", "svg"],
                        help="Specify the format of the output image. Default is pdf.")
    parser.add_argument("-output", type=str, help="Specify the output prefix for the image file.")
    parser.add_argument("-flag", type=int, default=4, help="Specifies the value of the flag parameter. Default is 4.")

    args = parser.parse_args()

    alignment = args.alignment
    reference = args.reference

    sequencing_type = ""
    if args.single:
        sequencing_type = "single"
        single_fastq = args.single
    elif args.pair:
        sequencing_type = "pair"
        pair_fastq1, pair_fastq2 = args.pair
    elif args.third:
        sequencing_type = "third"
        third_fastq = args.third

    genome_type = "line"  # deals with circular genomes

    output_prefix = "default_output_prefix"
    if args.output:
        output_prefix = args.output

    output_format = args.format

    input_fastq = ""
    if sequencing_type == "single":
        input_fastq = single_fastq
    elif sequencing_type == "pair":
        input_fastq = (pair_fastq1, pair_fastq2)
    elif sequencing_type == "third":
        input_fastq = third_fastq

    # Create the output folder if it doesn't exist
    output_folder = output_prefix + "_results"
    os.makedirs(output_folder, exist_ok=True)

    # Copy the reference file to the output folder
    shutil.copy(reference, output_folder)

    # Pass the 'flag' parameter to the draw_seq_depth function
    draw_seq_depth(alignment, reference, sequencing_type, genome_type, output_format, output_prefix, input_fastq, args.flag)
    #print(alignment, reference, sequencing_type, genome_type, output_format, output_prefix, input_fastq, args.flag)
