#!/usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import sys
import time

matplotlib.use('Agg')

# Read data
if len(sys.argv) != 4 or sys.argv[1] in ['-h', '--help']:
    print("Usage: plot_coverage.py <input_data> <output_data> <output_format>")
    print("Arguments:")
    print("  <input_data>      Path to the input data file.")
    print("  <output_data>     Path to save the output plot file.")
    print("  <output_format>   Output plot format (e.g., png, pdf).")
    sys.exit()

input_data = sys.argv[1]
output_data = sys.argv[2]
output_format = sys.argv[3]

print(input_data)
time.sleep(10)
print(output_data)
time.sleep(10)
print(output_format)
time.sleep(10)

column_names = ['Genome Name', 'Genomic Position', 'Coverage Depth']
data = pd.read_csv(input_data, sep='\t', header=None, names=column_names)

# Get data
x = data['Genomic Position'].values
y = data['Coverage Depth'].values

# Calculate average coverage depth, maximum depth, and minimum depth
average_coverage_depth = y.mean()
maximal_depth = y.max()
minimal_depth = y.min()
total_genome_length = x[-1]

# Count coverage with zero reads
zero_coverage_count = sum(y == 0)

# Create the plot
plt.figure(figsize=(12, 6))
plt.plot(x, y, color='blue', linewidth=1, linestyle='-')

# Customize the plot
plt.title('Sequencing Depth and Coverage Map', fontsize=18, pad=16)
plt.xlabel('Genomic Position (bp)', fontsize=16)
plt.ylabel('Sequencing Depth (×)', fontsize=16)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, linestyle='--', alpha=0.5)

# Add text annotations
total_genome_length_str = f'{total_genome_length:,} bp'
text_str_1 = f'(1) Total sequence length = {total_genome_length_str}'
text_str_3 = f'(3) Maximum depth = {maximal_depth} ×'
text_str_2 = f'(2) Average depth = {average_coverage_depth:.2f} ×'
text_str_4 = f'(4) Minimum depth = {minimal_depth} ×'
text_str_5 = f'(5) Number of bases not covered by any reads: {zero_coverage_count:,} bp'

plt.gca().annotate(text_str_1, fontsize=16, xy=(0, -0.15), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_3, fontsize=16, xy=(0, -0.22), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_5, fontsize=16, xy=(0, -0.29), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_2, fontsize=16, xy=(0.5, -0.15), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_4, fontsize=16, xy=(0.5, -0.22), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')

# Adjust layout to fit the text
plt.subplots_adjust(bottom=0.2)

# Save the plot
plt.savefig(output_data + "." + output_format, bbox_inches='tight', dpi=600)

print("Plot saved as", output_data + "." + output_format)