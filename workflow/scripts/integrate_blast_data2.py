import pandas as pd
import gzip
import sys

def integrate_blast_data(card_results_path, silva_results_path, aro_mapping_path, taxa_mapping_path, output_path):
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    
    # Read in chunks to avoid memory issues
    chunksize = 1000  # Adjust based on memory availability

    aro_df = pd.read_csv(aro_mapping_path, sep="\t")
    taxa_df = pd.read_csv(taxa_mapping_path, sep="\t")

    # Flag for writing the header once
    first_chunk = True

    # Process CARD results in chunks
    with gzip.open(card_results_path, 'rt') as f:
        for chunk in pd.read_csv(f, sep="\t", names=blast_columns, chunksize=chunksize):
            chunk["part"] = "ABR"
            chunk['ARO Accession'] = chunk['subject_id'].str.split("|", expand=True)[2]
            chunk["orientation"] = chunk["q_start"] - chunk["q_end"]
            merged_chunk = chunk.merge(aro_df, on='ARO Accession', how='left')
            merged_chunk = process_orientation_and_counts(merged_chunk)
            merged_chunk.to_csv(output_path, mode='a', header=first_chunk, index=False)
            first_chunk = False  # After the first chunk, set to False

    # Process SILVA results in chunks
    with gzip.open(silva_results_path, 'rt') as f:
        for chunk in pd.read_csv(f, sep="\t", names=blast_columns, chunksize=chunksize):
            chunk["part"] = "16S"
            chunk["primaryAccession"] = chunk['subject_id'].str.split(".", expand=True)[0]
            chunk["orientation"] = chunk["q_start"] - chunk["q_end"]
            merged_chunk = chunk.merge(taxa_df, on="primaryAccession", how="left")
            merged_chunk = process_orientation_and_counts(merged_chunk)
            merged_chunk.to_csv(output_path, mode='a', header=False, index=False)  # No header after first write

def process_orientation_and_counts(group):
    orientation = "mixed"
    if (group["orientation"] < 0).all():
        orientation = "forward"
    elif (group["orientation"] >= 0).all():
        orientation = "reverse"
    most_common_q_start = group["q_start"].mode().iloc[0]
    most_common_q_end = group["q_end"].mode().iloc[0]
    
    # Add orientation as a column
    group["orientation"] = orientation
    group["most_common_q_start"] = most_common_q_start
    group["most_common_q_end"] = most_common_q_end
    return group

# Snakemake inputs/outputs
card_results = snakemake.input.card_results
silva_results = snakemake.input.silva_results
aro_mapping = snakemake.params.aro_mapping
taxa_mapping = snakemake.params.taxa_mapping
output = snakemake.output.integrated_data
integrate_blast_data(card_results, silva_results, aro_mapping, taxa_mapping, output)
