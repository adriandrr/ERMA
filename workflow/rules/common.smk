import os
from Bio import SeqIO

def get_base_dir():
    return config["base_dir"]

def get_fastq_dir():
    return os.path.join(get_base_dir(), "data", "fastq")

def get_blast_db_dir():
    return os.path.join(get_base_dir(), "data", "blast_db")

def get_output_dir():
    return os.path.join(get_base_dir(), "results")

def get_card_db_dir():
    return os.path.join(get_base_dir(), "data", "card_db")

def get_silva_db_dir():
    return os.path.join(get_base_dir(), "data", "silva_db")

def count_fasta_lines(fasta_file,maxlines):
    lines = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
    num_parts = lines // maxlines +1 
    return(num_parts)

def get_num_parts(wildcards):
    with open(f"{wildcards.base_dir}/data/fastq/{wildcards.sample}_num_parts.txt") as f:
        return (f.read().strip())