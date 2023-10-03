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


def cut_and_join(reference_file, output_folder):
    # Call the cut_join.py script to perform cut-and-join on the reference genome
    subprocess.call(["cut_join.py", reference_file])
    
    # Generate the path to the resulting circular sequence in the current directory
    circ_sequence_src = "junction.seq.fasta"
    
    # Generate the path to the target location in the _junction_results folder
    circ_sequence_dest = os.path.join(output_folder, "junction.seq.fasta")
    
    # move the generated circular sequence to the _junction_results folder
    shutil.move(circ_sequence_src, circ_sequence_dest)
    
    return circ_sequence_dest


def draw_seq_depth(alignment, reference, sequencing_type, genome_type, output_format, output_prefix, input_fastq, flag):
    # Get the directory of the script file
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Construct absolute paths for different scripts based on input parameters
    if sequencing_type == "single":
        script_path = os.path.join(script_dir, "minimap2.single.sh")
        if alignment == "minimap2":
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "bowtie2":
            script_path = os.path.join(script_dir, "bowtie2.single.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "bwa":
            script_path = os.path.join(script_dir, "bwa.single.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "hisat2":
            script_path = os.path.join(script_dir, "hisat2.single.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
    elif sequencing_type == "pair":
        script_path = ""
        if alignment == "minimap2":
            script_path = os.path.join(script_dir, "minimap2.pair.sh")
            subprocess.call(["bash", script_path, reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "bowtie2":
            script_path = os.path.join(script_dir, "bowtie2.pair.sh")
            subprocess.call(["bash", script_path, reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "bwa":
            script_path = os.path.join(script_dir, "bwa.pair.sh")
            subprocess.call(["bash", script_path, reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "hisat2":
            script_path = os.path.join(script_dir, "hisat2.pair.sh")
            subprocess.call(["bash", script_path, reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
    elif sequencing_type == "third":
        script_path = ""
        if alignment == "minimap2":
            script_path = os.path.join(script_dir, "minimap2.third.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "bowtie2":
            script_path = os.path.join(script_dir, "bowtie2.third.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "bwa":
            script_path = os.path.join(script_dir, "bwa.third.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])
        elif alignment == "hisat2":
            script_path = os.path.join(script_dir, "hisat2.third.sh")
            subprocess.call(["bash", script_path, reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_junction_results"), str(flag)])

    # Extract coverage.txt from the output folder and convert to absolute path
    coverage_file = os.path.abspath(os.path.join(os.getcwd(), output_prefix + "_junction_results", "coverage.txt"))
    print(coverage_file)

    # Loop over each output format
    for format in output_format:
        output_file = os.path.abspath(os.path.join(os.getcwd(), output_prefix + "_junction_results", output_prefix))
    
        # Call plot.py script and pass command line arguments
        plot_script_path = os.path.join(script_dir, "plot.py")
        subprocess.call([str(sys.executable), plot_script_path, coverage_file, output_file, format])
        print([str(sys.executable), plot_script_path, coverage_file, output_file, format])
        
    # Get the current path
    current_path = os.path.join(os.getcwd(), output_prefix + "_results")

    # Get all files and directories in the current path
    files = os.listdir(current_path)

    # Iterate through and print the absolute paths of files and directories
    for file in files:
        file_path = os.path.join(current_path, file)
        print(file_path)
    

def main():
    # Create argparse object and set up command line arguments
    parser = argparse.ArgumentParser(prog="DrawSeqDepth-circle.py", description="Draw sequence depth for circular sequences.")

    parser.add_argument("-alignment",
                        choices=["minimap2", "hisat2", "bowtie2", "bwa"],default="minimap2",
                        help="Specifies the alignment software to be used. Default is minimap2.")    
    parser.add_argument("-reference", type=str, help="Specifies the reference genome.")    
    parser.add_argument("-single", metavar="FASTQ", help="Specifies single-end sequencing data.")
    parser.add_argument("-pair", metavar=("FASTQ1", "FASTQ2"), nargs=2, help="Specifies paired-end sequencing data.")
    parser.add_argument("-third", metavar="FASTQ", help="Specifies third-generation sequencing data.")
    parser.add_argument("-format", default=["pdf"], nargs='+', choices=["png", "pdf","jpeg","bmp","tiff","gif","eps","svg"],
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

    genome_type = "circle"  # deals with circular genomes

    output_prefix = "default_output_results"
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
    output_folder = os.path.join(os.getcwd(), output_prefix + "_junction_results")
    os.makedirs(output_folder, exist_ok=True)

    # Call the cut_and_join function to get the junction_sequence
    junction_sequence = cut_and_join(reference, output_folder)
    
    # Call the draw_seq_depth function with modified parameters, including 'flag'
    draw_seq_depth(alignment, junction_sequence, sequencing_type, genome_type, output_format, output_prefix, input_fastq, args.flag)
    print(alignment, junction_sequence, sequencing_type, genome_type, output_format, output_prefix, input_fastq, args.flag)

if __name__ == "__main__":
    main()
