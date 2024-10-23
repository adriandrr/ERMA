import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

dtype_dict = {
    "taxid": "float",
    "ARO Accession": "string",
    "CVTERM ID": "float",
    "Model Sequence ID": "float",
    "Model ID": "float",
    "Model Name": "string",
    "ARO Name": "string",
    "Protein Accession": "string",
    "DNA Accession": "string",
    "AMR Gene Family": "string",
    "Drug Class": "string",
    "Resistance Mechanism": "string",
    "CARD Short Name": "string",
}

def generate_percentage_idt_per_genus(input_file, output_file):
    df = pd.read_csv(input_file,dtype=dtype_dict)
    df['species'] = df['path'].apply(lambda x: x.split(';')[-2] if pd.notna(x) else None)

    species_query_counts = df.groupby('species')['query_id'].nunique().reset_index()
    species_query_counts.columns = ['species', 'unique_query_count']
    df = pd.merge(df, species_query_counts, on='species')
    species_order = df.groupby('species')['perc_identity'].median().sort_values(ascending=False).index
    fig, ax1 = plt.subplots(figsize=(15, 8))

    sns.boxplot(x='species', y='perc_identity', data=df, ax=ax1, order=species_order)
    ax1.set_xlabel("Bacterial Species")
    ax1.set_ylabel("Percentage Identity")
    ax1.set_title("Boxplot of Percentage Identity and Read Counts for Each Bacterial Species")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
    ax2 = ax1.twinx()
    sns.barplot(x='species', y='unique_query_count', data=species_query_counts, ax=ax2, alpha=0.3, color='blue', order=species_order)
    ax2.set_ylabel("Number of Unique Reads (query_id)", color='blue')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    input_file = snakemake.input.filtered_data
    output_file = snakemake.output[0]
    sys.stderr = open(snakemake.log[0], "w")
    generate_percentage_idt_per_genus(input_file, output_file)