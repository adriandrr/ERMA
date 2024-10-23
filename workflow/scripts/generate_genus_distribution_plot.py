import pandas as pd
import altair as alt
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

def generate_genus_distribution_plot(input_file, output_file, sample_name):
    df_full = pd.read_csv(input_file, sep=",", header=0,dtype=dtype_dict)    
    abr = df_full[df_full["part"] == "ABR"]
    sixteen_s = df_full[df_full["part"] == "16S"]
    df = pd.merge(abr,sixteen_s, on='query_id', suffixes=('_abr','_16S'))
    df['genus'] = df['path_16S'].apply(lambda x: x.split(';')[-2] if pd.notna(x) else None)
    
    # Group by AMR gene family and genus
    grouped_data = df.groupby(['AMR Gene Family_abr', 'genus']).size().reset_index(name='count')
    
    # Create the Altair bar chart
    chart = alt.Chart(grouped_data).mark_bar().encode(
        x='AMR Gene Family_abr',
        y='count',
        color='genus',
        tooltip=['AMR Gene Family_abr', 'genus', 'count']
    ).properties(
        title=f'Genus Distribution for {sample_name}',
        width=800,
        height=400
    )
    
    # Save the plot as an HTML file
    chart.save(output_file)

if __name__ == "__main__":
    input_file = snakemake.input.filtered_data
    output_file = snakemake.output[0]
    sample_name = snakemake.params.sample_name
    sys.stderr = open(snakemake.log[0], "w")
    generate_genus_distribution_plot(input_file, output_file, sample_name)
