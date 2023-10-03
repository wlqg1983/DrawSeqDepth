#!/usr/bin/env python
import os
import requests
import sys
import argparse

def download_fasta_sequence(accession):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        sys.stderr.write(f"Failed to download FASTA sequence for accession: {accession}\n")
        return None

def download_fasta_sequences(accessions, output_dir):
    for accession in accessions:
        fasta_sequence = download_fasta_sequence(accession)
        if fasta_sequence:
            filename = f"{accession}.fasta"
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                filepath = os.path.join(output_dir, filename)
            else:
                filepath = filename
            with open(filepath, 'w') as file:
                file.write(fasta_sequence)
            print(f"Downloaded and saved {accession} to {filepath}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download FASTA sequences from NCBI')
    parser.add_argument('accessions', nargs='+', help='Accession numbers of the sequences')
    parser.add_argument('-o', '--output', help='Output directory for saving FASTA files')

    # Check if -h or --help is provided
    if '-h' in sys.argv or '--help' in sys.argv:
        parser.print_help()
    else:
        try:
            args = parser.parse_args()
            download_fasta_sequences(args.accessions, args.output)
        except argparse.ArgumentError:
            parser.print_help()