#!/usr/bin/env python
import argparse
import subprocess
import time
import sys


def main():
    # Create an argparse object and set up command line arguments
    parser = argparse.ArgumentParser(prog="DrawSeqDepth.py", description="DrawSeqDepth Main Program.")
    
    # Add -line and -circle options
    parser.add_argument("-line", action="store_true", help='Process linear genomes')
    parser.add_argument("-circle", action="store_true", help='Process circular genomes')

    # Add other parameters
    parser.add_argument("-alignment",
                        choices=["minimap2", "hisat2", "bowtie2", "bwa"], default="minimap2",
                        help="Specify the alignment software to be used. Default is minimap2.")
    parser.add_argument("-reference", type=str, help="Specify the reference genome.")
    parser.add_argument("-single", metavar="FASTQ", help="Specify single-end sequencing data.")
    parser.add_argument("-pair", metavar=("FASTQ1", "FASTQ2"), nargs=2, help="Specify paired-end sequencing data.")
    parser.add_argument("-third", metavar="FASTQ", help="Specify third-generation sequencing data.")
    parser.add_argument("-format", default="pdf", choices=["png", "pdf", "jpeg", "bmp", "tiff", "gif", "eps", "svg"],
                        help="Specify the output image format. Default is pdf.")
    parser.add_argument("-output", type=str, help="Specify the output prefix for image files.")
    parser.add_argument("-flag", type=int, default=4, help="Specify the value of the flag parameter. Default is 4.")

    # Parse command line arguments
    args = parser.parse_args()

    # Check if no arguments are given, print help and exit
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(0)

    # Check if any required parameters are missing
    if args.reference is None or args.output is None or (args.single is None and args.pair is None and args.third is None):
        parser.print_help()
        sys.exit(0)

    # Check if both -line and -circle options are selected
    if args.line and args.circle:
        warning = "\033[93mWarning: Please use one of -line/-circle one time.\033[0m"
        print(warning)
        return

    # Check if third-generation sequencing data is selected without minimap2 alignment software
    if args.third and args.alignment != "minimap2":
        warning = "\033[93mWarning: Using this software for aligning the third-generation sequencing data to the reference genome/nucleotide/gene will be computationally time-consuming. Recommended to use minimap2 for alignment.I will do all the work still.\033[0m"
        print(warning)
        time.sleep(10)

    # Build the argument list
    arg_list = []

    # Add other parameters
    for arg_name, arg_value in vars(args).items():
        if arg_name == 'pair' and arg_value is not None:
            arg_list.extend([f'-pair', arg_value[0], arg_value[1]])
        elif arg_name not in ["line", "circle"] and arg_value is not None:
            arg_list.extend([f'-{arg_name}', str(arg_value)])

    if args.line:
        # When '-line' is chosen, call DrawSeqDepth-line.py
        subprocess.call(["DrawSeqDepth-line.py"] + arg_list)
    else:
        # If neither -line nor -circle is specified, or '-line' is explicitly selected
        # Call DrawSeqDepth-line.py
        subprocess.call(["DrawSeqDepth-line.py"] + arg_list)

    if args.circle:
        # When '-circle' is chosen, call DrawSeqDepth-circle.py
        subprocess.call(["DrawSeqDepth-circle.py"] + arg_list)

if __name__ == "__main__":
    main()
    print("\033[93m\nThe program has finished! Please check the results!!\n\033[0m")