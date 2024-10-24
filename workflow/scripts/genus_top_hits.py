import dask.dataframe as pd
import altair as alt
import sys
import gzip

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

# Define a function to create plots for a given approach and top N
def create_plots(merged_data, title_suffix):
    # Sort merged_data by 'query_id' and 'perc_identity_16s' in descending order
    sorted_data = merged_data.sort_values(['query_id', 'perc_identity_16s'], ascending=[True, False])
    plots = []
    
    # Loop through top 1, 3, 5, and 10 hits
    for top_n in [1, 3, 5, 10]:
        # Select the top N hits for each query_id
        top_data = sorted_data.groupby('query_id').head(top_n)
        # Count unique query_ids per genus and AMR Gene Family
        final_grouped_data = top_data.groupby(['AMR Gene Family_abr', 'genus']).agg({'query_id': 'nunique'}).reset_index()
        final_grouped_data.rename(columns={'query_id': 'unique_query_count'}, inplace=True)
        # Total counts for ordering
        total_counts = final_grouped_data.groupby('genus')['unique_query_count'].sum().reset_index()
        total_counts = total_counts.sort_values(by='unique_query_count', ascending=False)
        # Set ordered categories for the genus column
        ordered_genus = pd.Categorical(final_grouped_data['genus'], categories=total_counts['genus'], ordered=True)
        final_grouped_data['genus'] = ordered_genus

        # Create the bar plot
        bars = alt.Chart(final_grouped_data).mark_bar().encode(
            x=alt.X('AMR Gene Family_abr', title='AMR Gene Family'),
            y=alt.Y('unique_query_count', title='Unique Query Count', axis=alt.Axis(labelAngle=0)),
            color=alt.Color('genus', title='Genus'),
            tooltip=['AMR Gene Family_abr', 'genus', 'unique_query_count']
        ).properties(
            title=f'Top {top_n} Hits {title_suffix}',
            width=200,
            height=400
        )
        
        plots.append(bars)
    
    # Combine the plots vertically
    combined_chart = alt.hconcat(*plots).configure_axis(labelAngle=45)
    
    return combined_chart

def process_blast_data(input_file,output_file,sample_name):
    filtered_df = pd.read_csv(input_file, sep=",", header=0,dtype=dtype_dict)

    # Filter ABR and 16S parts with respective conditions
    abr_data = filtered_df[filtered_df['part'] == 'ABR']
    sixteen_s_data = filtered_df[filtered_df['part'] == '16S']
    
    # Merge ABR and 16S data on 'query_id'
    merged_data_path_16s = pd.merge(abr_data, sixteen_s_data, on='query_id', suffixes=('_abr', '_16s'))
    
    # Extract genus from the path
    merged_data_path_16s['genus'] = merged_data_path_16s['path_16s'].apply(lambda x: x.split(';')[-2] if pd.notna(x) else None)
    
    # Create the plot for genus comparison
    plot_path_16s = create_plots(merged_data_path_16s, 'Genus')

    # Configure and display the combined charts
    combined_charts = plot_path_16s.resolve_scale(
        color='independent'
    ).configure_legend(
        titleFontSize=12,
        labelFontSize=10,
        orient='right'
    ).configure_title(
        fontSize=20
    ).configure_axis(
        labelAngle=45
    ).properties(
        title=alt.TitleParams(
            f'{sample_name}: Comparison of Genus for AMR Gene Families (Top 1, 3, 5, 10)', 
            anchor='middle'
            )
    )

    combined_charts.configure_title(
        fontSize=20,
        font='Courier',
        anchor='start',
        color='gray'
    )

    combined_charts.save(output_file)

if __name__ == "__main__":
    input_file = snakemake.input.filtered_data
    output_file = snakemake.output[0]
    sample_name = snakemake.params.sample_name
    sys.stderr = open(snakemake.log[0], "w")  
    process_blast_data(input_file,output_file,sample_name)
