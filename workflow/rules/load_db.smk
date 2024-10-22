rule get_16S_db:
    output:
        seq = "{base_dir}/data/silva_db/silva_seq.fasta.gz",
        tax = "{base_dir}/data/silva_db/silva_tax.txt.gz"
    params:
        seq = config["silva"]["download-path-seq"],
        tax = config["silva"]["download-path-tax"],
        path = "{base_dir}/data/silva_db"
    conda:
        "../envs/python.yaml"
    shell:
        """
        mkdir -p {params.path};
        cd {params.path};
        wget -O silva_seq.fasta.gz {params.seq};
        wget -O silva_tax.txt.gz {params.tax};
        """

rule unzip_silva_db:
    input:
        seq = "{base_dir}/data/silva_db/silva_seq.fasta.gz",
        tax = "{base_dir}/data/silva_db/silva_tax.txt.gz"
    output:
        seq = "{base_dir}/data/silva_db/silva_seq.fasta",
        tax = "{base_dir}/data/silva_db/silva_tax.txt"
    shell:
        """
        gzip -dk {input.seq};
        gzip -dk {input.tax};
        """

rule get_card_db:
    output:
        seq = "{base_dir}/data/card_db/card_seq.tar.bz2"
    params:
        seq = config["card"]["download-path"],
        path = "{base_dir}/data/card_db"
    shell:
        """
        mkdir -p {params.path};
        cd {params.path};
        wget -O card_seq.tar.bz2 {params.seq};
        """

rule unzip_card_db:
    input:
        seq = "{base_dir}/data/card_db/card_seq.tar.bz2"
    output:
        seq = "{base_dir}/data/card_db/protein_fasta_protein_homolog_model.fasta"
    params:
        path = "{base_dir}/data/card_db"
    shell:
        """
        tar -xvjf {input.seq} -C {params.path};
        """

rule makeblastdb_card:
    input:
        seq = "{base_dir}/data/card_db/protein_fasta_protein_homolog_model.fasta"
    output:
        db = "{base_dir}/data/blast_db/card_db.pdb"
    params:
        path = "{base_dir}/data/blast_db/card_db"    
    shell:
        """
        makeblastdb -in {input.seq} -dbtype prot -out {params.path};
        """

rule makeblastdb_silva:
    input:
        seq = "{base_dir}/data/silva_db/silva_seq.fasta"
    output:
        db = "{base_dir}/data/blast_db/silva_db.ndb"
    params:
        path = "{base_dir}/data/blast_db/silva_db"           
    shell:
        """
        makeblastdb -in {input.seq} -dbtype nucl -out {params.path};
        """
