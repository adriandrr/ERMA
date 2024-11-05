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

def get_numpart_list():
    return(list([f"00{i}" if i < 10 else f"0{i}" for i in range(1,config["num_parts"]+1)]))

def get_num_parts(wildcards):
    with open(f"{wildcards.base_dir}/data/fastq/{wildcards.sample}_num_parts.txt") as f:
        return int(f.read().strip())