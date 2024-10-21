import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor

def generate_read_position(input_file, output_file, sample_name):
    # Load the data
    df = pd.read_csv(input_file, sep=",", header=0)
    df = df[['query_id', 'most_common_q_start', 'most_common_q_end', 'part']]

    abr_data = df[df['part'] == 'ABR']
    sixteen_s_data = df[df['part'] == '16S']

    # Merge ABR and 16S data based on query_id
    merged_data = pd.merge(abr_data, sixteen_s_data, on='query_id', suffixes=('_abr', '_16s'))
    merged_data['x_enum'] = np.arange(1, len(merged_data) + 1)

    # Plotting
    plt.figure(figsize=(10, 6))
    for i, row in merged_data.iterrows():
        x_val = row['x_enum']
        
        # Use the most_common_q_start and most_common_q_end directly from the dataframe
        abr_start = row['most_common_q_start_abr']
        abr_end = row['most_common_q_end_abr']
        abr_length = abr_end - abr_start
        
        sixteen_s_start = row['most_common_q_start_16s']
        sixteen_s_end = row['most_common_q_end_16s']
        sixteen_s_length = sixteen_s_end - sixteen_s_start
        
        # Determine the overlap (if any)
        overlap_start = max(abr_start, sixteen_s_start)
        overlap_end = min(abr_end, sixteen_s_end)
        
        # Plot ABR
        plt.bar(x_val, abr_length, bottom=abr_start, label='ABR' if i == 0 else "", color='lightblue')
        
        # Plot 16S
        plt.bar(x_val, sixteen_s_length, bottom=sixteen_s_start, label='16S' if i == 0 else "", color='orange')
        
        # Plot overlap if it exists
        if overlap_start < overlap_end:
            overlap_length = overlap_end - overlap_start
            plt.bar(x_val, overlap_length, bottom=overlap_start, color='red', label='Overlap' if i == 0 else "")
    
    plt.xlabel('Query ID Enumeration')
    plt.ylabel('Position')
    plt.title(f'{sample_name} - Stacked Bar Plot: most_common_q_start and most_common_q_end for ABR and 16S Parts with Overlap Highlighted')
    plt.xticks(np.linspace(1, len(merged_data), 5).astype(int))
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    input_file = snakemake.input.filtered_data
    output_file = snakemake.output[0]
    sample_name = snakemake.params.sample_name    
    generate_read_position(input_file, output_file, sample_name)
