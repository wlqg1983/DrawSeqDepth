#!/usr/bin/env python
import sys

def print_help():
    print("Usage: cut_join.py <reference_file> [result_file_name]\n")
    print("Arguments:")
    print("  <reference_file>       Path to the reference file containing sequence data.")
    print("  [result_file_name]     Name for the result file (optional, defaults to 'junction.seq.fasta').")
    print("Options:")
    print("  -h, --help             Show this help message and exit.\n")
    print("Description:")
    print("  This program reads sequence data from the specified reference file and performs the following steps:")
    print("  1. Extracts the first 1000 bases from the 5' end and the last 1000 bases from the 3' end of the reference sequence.")
    print("  2. Joins the extracted sequences to create a junction sequence.")
    print("  3. Saves the junction sequence to a FASTA file.\n")
    print("Example:")
    print("cut_join.py reference.fasta")
    print("OR:")
    print("cut_join.py reference.fasta my_junction.fasta")

def cut_and_join(reference_file, output_file_name="junction.seq.fasta"):
    with open(reference_file, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip()
        sequence = ''.join(lines[1:]).replace(' ', '').replace('\n', '')

    # Extract the 5' and 3' ends of the reference sequence
    five_prime_end = sequence[:1000]
    three_prime_end = sequence[-1000:]

    # Join the two ends to create the junction sequence
    junction_sequence = five_prime_end + three_prime_end

    # Save the junction sequence to a FASTA file
    with open(output_file_name, 'w') as file:
        file.write(f'>{header}\n{junction_sequence}\n')

    print(f'Junction sequence saved to {output_file_name}')

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        if sys.argv[1] in ['-h', '--help']:
            print_help()
            sys.exit(0)

        reference_file = sys.argv[1]
        output_file_name = "junction.seq.fasta" if len(sys.argv) < 3 else sys.argv[2]
        cut_and_join(reference_file, output_file_name)
    else:
        print_help()
