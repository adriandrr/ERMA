# EpicPCR Sequencing Analysis Pipeline

This pipeline is designed to process sequencing reads obtained from epicPCR experiments, linking antimicrobial resistance (AMR) genes with 16S rRNA genes. The pipeline uses **Snakemake** to manage the workflow and integrates tools for downloading necessary databases, running sequence alignments, filtering, and generating visual reports.

## Features
- Downloads and prepares SILVA and CARD databases for use in the analysis.
- Converts raw FASTQ sequencing files to FASTA format.
- Performs BLAST sequence alignments against both CARD (AMR genes) and SILVA (16S rRNA).
- Integrates BLAST results to identify linked AMR and microbial markers.
- Generates filtered results based on alignment quality and similarity thresholds.
- Produces graphical outputs such as genus distribution plots, e-value boxplots, and read positions.
- Generates an HTML report summarizing the analysis.

## Prerequisites

### Install Snakemake
The pipeline requires **Snakemake** with support for conda environments. You can install it using conda (via Miniconda or Anaconda):

```bash
	conda install -c bioconda -c conda-forge snakemake
```

Alternatively, follow the official [[!Snakemake installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)] guide for more options.

Install Dependencies
This pipeline uses conda environments to manage its dependencies. Snakemake will automatically create and manage these environments when run with the 
--use-conda 
flag.

## Usage Instructions

Clone the repository: First, clone the pipeline repository to your local machine:
bash
Copy code
git clone https://github.com/your-username/epicPCR-pipeline.git
cd epicPCR-pipeline
Prepare Data Folder: You need to place your raw sequencing files (FASTQ format) in the data/fastq/ directory. This folder must exist before running the pipeline.
Modify the Config File: Open the config/config.yaml file and change the base_dir parameter to the base directory where the pipeline is located. The config file should look like this:
```yaml
Copy code
base_dir: "/path/to/your/project"

min_similarity: "0.8"

silva:
  download-path-seq: "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSUParc_tax_silva.fasta.gz"
  download-path-tax: "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_parc_138.2.txt.gz"

card:
  download-path: "https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2"
```
Replace /path/to/your/project with the actual path to your local pipeline directory.
Run the Pipeline: To start the pipeline, run the following command from the base directory:

```bash
Copy code
snakemake --use-conda --cores N
```
Replace N with the number of cores (threads) you want to use.

The pipeline will automatically:

- Download and prepare the required SILVA and CARD databases.
- Decompress and convert your raw FASTQ files.
- Perform BLAST alignments against the CARD and SILVA databases.
- Generate integrated, filtered results.
- Produce plots and an HTML report summarizing the analysis.
- View Results: After the pipeline finishes running, results will be saved in the results/ directory. You can find the final report as a ZIP file in results/report.zip. This file contains HTML reports and visual summaries of the analysis.

## Additional Notes

The pipeline is designed to handle large sequencing datasets in parallel, so it's recommended to run it on a machine with sufficient computational resources.
If any errors occur during the pipeline run, Snakemake will provide detailed logs, allowing you to debug and troubleshoot any issues.
License

This project is licensed under the MIT License.
