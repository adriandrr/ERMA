import pandas as pd
import altair as alt
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
    "AMR Gene Family_abr": "string",
}

def create_plots(merged_data, title_suffix):
    plots = []
    
    for top_n in [1, 3, 5, 10]:
        top_data = merged_data.groupby('query_id').head(top_n)
        final_grouped_data = top_data.groupby(['AMR Gene Family_abr', 'genus']).agg({'query_id': 'nunique'}).reset_index()
        final_grouped_data.rename(columns={'query_id': 'unique_query_count'}, inplace=True)
        total_counts = final_grouped_data.groupby('genus')['unique_query_count'].sum().reset_index()
        total_counts = total_counts.sort_values(by='unique_query_count', ascending=False)
        ordered_genus = pd.Categorical(final_grouped_data['genus'], categories=total_counts['genus'], ordered=True)
        final_grouped_data['genus'] = ordered_genus

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
    
    combined_chart = alt.hconcat(*plots).configure_axis(labelAngle=45)
    
    return combined_chart

def process_blast_data(input_files, output_file, sample_name,top_hits):
    merged_data_list = []
    query_id_set = set()  # Set to track query_ids
    overlapping_ids = set()  # Set to track overlapping query_ids

    for input_file in input_files:
        df = pd.read_csv(input_file, sep=",", usecols=necessary_columns, header=0, dtype=dtype_dict)
        
        # Check for overlapping query IDs
        current_query_ids = set(df['query_id'])
        overlap = query_id_set.intersection(current_query_ids)
        if overlap:
            overlapping_ids.update(overlap)
        query_id_set.update(current_query_ids)
        
        abr_data = df[df['part'] == 'ABR']
        sixteen_s_data = df[df['part'] == '16S']
        
        # Merge based on 'query_id'
        merged_df = pd.merge(abr_data, sixteen_s_data, on='query_id', suffixes=('_abr', '_16s'))
        merged_df['genus'] = merged_df['path_16s'].apply(lambda x: x.split(';')[-2] if pd.notna(x) else None)
        sorted_data = merged_df.sort_values(['query_id', 'perc_identity_16s'], ascending=[True, False])
        top_df = sorted_data.groupby('query_id').head(top_hits[-1])
        merged_data_list.append(top_df)

    # Concatenate all merged data
    merged_data_path_16s = pd.concat(merged_data_list, ignore_index=True)

    plot_path_16s = create_plots(merged_data_path_16s, 'Genus')

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

    # Log overlapping query_ids if any
    if overlapping_ids:
        print(f"Overlapping query_ids detected: {overlapping_ids}", file=sys.stderr)

if __name__ == "__main__":
    input_files = snakemake.input.filtered_data
    output_file = snakemake.output[0]
    sample_name = snakemake.params.sample_name
    show_top_hits = [1,3,5,10] # largest number needs to be last
    sys.stderr = open(snakemake.log[0], "w")  
    process_blast_data(input_files, output_file, sample_name,show_top_hits)
