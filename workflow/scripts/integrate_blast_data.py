import dask.dataframe as dd
import pandas as pd
import gzip
import sys
import os

def process_orientation_and_counts(group):
    orientation = "mixed"
    if (group["distance"] < 0).all():
        orientation = "reverse"
    if (group["distance"] >= 0).all():
        orientation = "forward"
    return pd.Series({"orientation": orientation})

def process_card_results(card_results_path, aro_mapping_path, output_path):
    """Process CARD results using Dask and save them to an intermediate output file."""
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    
    aro_df = pd.read_csv(aro_mapping_path, sep="\t")  # Read ARO mapping with pandas (small dataset)

    # Read large CARD data with Dask
    card_df = dd.read_csv(card_results_path, sep="\t", names=blast_columns, compression='gzip')
    card_df["part"] = "ABR"
    card_df['ARO Accession'] = card_df['subject_id'].str.split("|", expand=True)[2]
    card_df["distance"] = card_df["q_start"] - card_df["q_end"]

    # Apply orientation counts using map_partitions (works on each partition of data)
    orientation_counts = card_df.groupby("query_id").apply(process_orientation_and_counts, meta=pd.DataFrame()).reset_index()

    # Merge with orientation counts and ARO mapping
    merged_df = card_df.merge(orientation_counts, on="query_id")
    merged_df = merged_df.merge(aro_df, on='ARO Accession', how='left')

    # Save results to CSV
    merged_df.compute().to_csv(output_path, index=False)

def process_silva_results(silva_results_path, taxa_mapping_path, output_path):
    """Process SILVA results using Dask and save them to an intermediate output file."""
    blast_columns = ["query_id", "subject_id", "perc_identity", "align_length", "mismatches",
                     "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    
    taxa_df = pd.read_csv(taxa_mapping_path, sep="\t")  # Read Taxa mapping with pandas (small dataset)

    # Read large SILVA data with Dask
    silva_df = dd.read_csv(silva_results_path, sep="\t", names=blast_columns, compression='gzip')
    silva_df["part"] = "16S"
    silva_df[["primaryAccession", "acc_start", "acc_stop"]] = silva_df['subject_id'].str.split(".", expand=True)
    silva_df["distance"] = silva_df["q_start"] - silva_df["q_end"]

    # Apply orientation counts using map_partitions (works on each partition of data)
    orientation_counts = silva_df.groupby("query_id").apply(process_orientation_and_counts, meta=pd.DataFrame()).reset_index()

    # Merge with orientation counts and Taxa mapping
    merged_df = silva_df.merge(orientation_counts, on="query_id")
    merged_df = merged_df.merge(taxa_df, on="primaryAccession", how="left")

    # Save results to CSV
    merged_df.compute().to_csv(output_path, index=False)

def merge_results(card_output, silva_output, final_output):
    """Merge processed CARD and SILVA results using pandas into one final output file."""
    card_df = pd.read_csv(card_output)
    silva_df = pd.read_csv(silva_output)
    
    combined_df = pd.concat([silva_df, card_df])
    
    combined_df.to_csv(final_output, index=False)

if __name__ == "__main__":
    card_results = snakemake.input.card_results
    silva_results = snakemake.input.silva_results
    aro_mapping = snakemake.input.aro_mapping
    taxa_mapping = snakemake.input.taxa_mapping
    card_output = snakemake.output.intermed_card_results
    silva_output = snakemake.output.intermed_silva_results
    final_output = snakemake.output.integrated_data
    sys.stderr = open(snakemake.log[0], "w")

    # Process CARD and SILVA results using Dask
    process_card_results(card_results, aro_mapping, card_output)
    process_silva_results(silva_results, taxa_mapping, silva_output)

    # Merge results using pandas
    merge_results(card_output, silva_output, final_output)
    
    print(f"Processing complete. Final merged output saved to {final_output}")
