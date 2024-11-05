import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

necessary_columns = [
    "query_id",
    "part",
    "path",
    "AMR Gene Family",
    "perc_identity",
]

dtype_dict = {
    "query_id": "string",
    "part": "string",
    "path": "string",
    "AMR Gene Family": "string",
    "perc_identity": "float"
}

def generate_percentage_idt_per_genus(input_files, output_file):
    all_data = []  # List to hold DataFrames from all input files

    for input_file in input_files:
        df = pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, dtype=dtype_dict)
        df['species'] = df['path'].apply(lambda x: x.split(';')[-2] if pd.notna(x) else None)
        all_data.append(df)
        
    # Combine all partitions into a single DataFrame
    combined_data = pd.concat(all_data)
    
    # Calculate species query counts and species order
    species_query_counts = combined_data.groupby('species')['query_id'].nunique().reset_index()
    species_query_counts.columns = ['species', 'unique_query_count']
    combined_data = pd.merge(combined_data, species_query_counts, on='species')
    species_order = combined_data.groupby('species')['perc_identity'].median().sort_values(ascending=False).index
    
    # Plotting
    fig, ax1 = plt.subplots(figsize=(15, 8))
    sns.boxplot(x='species', y='perc_identity', data=combined_data, ax=ax1, order=species_order)
    ax1.set_xlabel("Bacterial Species")
    ax1.set_ylabel("Percentage Identity")
    ax1.set_title("Boxplot of Percentage Identity and Read Counts for Each Bacterial Species")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

    # Add a second y-axis for unique query counts
    ax2 = ax1.twinx()
    sns.barplot(x='species', y='unique_query_count', data=species_query_counts, ax=ax2, alpha=0.3, color='blue', order=species_order)
    ax2.set_ylabel("Number of Unique Reads (query_id)", color='blue')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    input_files = snakemake.input.filtered_data  # List of partitioned files
    output_file = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    generate_percentage_idt_per_genus(input_files, output_file)
