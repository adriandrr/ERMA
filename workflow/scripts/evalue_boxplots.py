import os
import sys
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_and_process_data(csv_file, min_similarity):
    """Read and process the integrated_filtered_results.csv.gz file."""
    
    if os.path.exists(csv_file):
        fields = ['evalue','part']
        df = pd.read_csv(csv_file, header=0, sep=',',usecols=fields)
        sample_name = csv_file.split("/")[-2]
        df['sample'] = sample_name
        return df
    
    print(f"File {csv_file} does not exist.")
    return None

def plot_boxplots(data, output_file):
    """Plot boxplots based on the e-values for ABR and 16S parts across samples."""
    plt.figure(figsize=(15, 10))
    flierprops = dict(markerfacecolor='0.75', markersize=2, linestyle='none')
    sns.boxplot(x='sample', y='evalue', hue='part', data=data, flierprops=flierprops)
    plt.yscale('log')
    plt.title('Boxplot of e-values for ABR and 16S parts across samples -Filtered-')
    plt.xlabel('Sample')
    plt.ylabel('e-value (log scale)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main(csv_files, min_similarity, output_dir):
    """Main function to process the CSV files and generate the plot."""
    output_plot_path = output_dir
    all_data = []

    # Loop over each CSV file, process it, and append the data
    for csv_file in csv_files:
        data = read_and_process_data(csv_file, min_similarity)
        if data is not None:
            all_data.append(data)

    if all_data:
        combined_data = pd.concat(all_data)
        print("Creating combined boxplot")
        plot_boxplots(combined_data, output_plot_path)
        print(f"Combined boxplot has been saved to {output_plot_path}")
    else:
        print("No data found.")

if __name__ == "__main__":
    # Assuming snakemake is used to pass inputs
    csv_files = snakemake.input.csv_files  # List of integrated_filtered_results.csv.gz files
    output_dir = snakemake.output[0]  # Directory where output will be saved
    min_similarity = snakemake.params.min_similarity  # Minimum similarity filter

    main(csv_files, min_similarity, output_dir)