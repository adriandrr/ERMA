import pandas as pd

def filter_group(group):
    most_common_q_start = group["most_common_q_start"].dropna().iloc[0]
    abr_part = group[group["part"] == "ABR"]
    filtered_16s_part = group[
        (group["part"] == "16S") & 
        (group["q_start"].between(most_common_q_start - 30, most_common_q_start + 30))
    ]
    if filtered_16s_part.empty:
        return pd.DataFrame()     
    return pd.concat([abr_part, filtered_16s_part])

def process_orientation_and_counts(group):
    most_common_q_start = group["q_start"].mode().iloc[0]
    most_common_q_end = group["q_end"].mode().iloc[0]
    return pd.Series({"most_common_q_start": most_common_q_start, "most_common_q_end": most_common_q_end})

def filter_blast_results(input_file, output_file):
    df = pd.read_csv(input_file, compression='gzip', header=0, sep=',')
    
    # Filter results based on percentage identity and alignment length for ABR
    abr_data = df[
        (df['part'] == 'ABR') &
        (df["perc_identity"] > 93.0) &
        (df['orientation'] == 'forward')
    ]
    orientation_counts = abr_data.groupby("query_id").apply(process_orientation_and_counts).reset_index()
    abr_data = abr_data.merge(orientation_counts, on="query_id")  # Assign the merged data back to abr_data

    # Filter results based on percentage identity and alignment length for 16S
    sixteen_s_data = df[
        (df['part'] == '16S') &
        (df['align_length'].between(200, 300)) &
        (df["perc_identity"] > 87.0) &
        (df['orientation'] == 'reverse')
    ]
    
    orientation_counts = sixteen_s_data.groupby("query_id").apply(process_orientation_and_counts).reset_index()
    sixteen_s_data = sixteen_s_data.merge(orientation_counts, on="query_id")  # Assign the merged data back to sixteen_s_data
    
    # Filter for common query IDs between abr_data and sixteen_s_data
    common_query_ids = pd.Index(abr_data['query_id']).intersection(sixteen_s_data['query_id'])
    abr_data_filtered = abr_data[abr_data['query_id'].isin(common_query_ids)]
    sixteen_s_data_filtered = sixteen_s_data[sixteen_s_data['query_id'].isin(common_query_ids)]
    
    # Concatenate the filtered ABR and 16S data
    merged_data = pd.concat([abr_data_filtered, sixteen_s_data_filtered])
    
    # Write the filtered data to output
    merged_data.to_csv(output_file, index=False)

if __name__ == "__main__":
    input_file = snakemake.input.integrated_data
    output_file = snakemake.output.filtered_data
    filter_blast_results(input_file, output_file)
