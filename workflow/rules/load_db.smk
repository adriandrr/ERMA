rule get_16S_db:
        output:
            seq="data/silva_db/silva_seq.fasta.gz",
            tax="data/silva_db/silva_tax.txt.gz",
        params:
            seq=str(config["silva"]["download-path-seq"]),
            tax=str(config["silva"]["download-path-tax"]),
            path=str(config["blast_db"]["silva"])
        conda:
            "../envs/python.yaml"
        shell:
            "mkdir {params.path}; "
            "cd {params.path}; "
            "wget -O silva_seq.fasta.gz {params.seq}; "
            "wget -O silva_tax.txt.gz {params.tax}; "

rule unzip_silva_db:
    input:
        seq="data/silva_db/silva_seq.fasta.gz",
        tax="data/silva_db/silva_tax.txt.gz",
    output:
        seq=temp("data/silva_db/silva_seq.fasta"),
        tax=temp("data/silva_db/silva_tax.txt"),
    shell:
        """
        gzip -dk {input.seq}
        gzip -dk {input.tax}
        """

rule get_card_db:
        output:
            seq="data/card_db/card_seq.tar.bz2",
        params:
            seq=str(config["card"]["download-path-seq"]),
            path=str(config["blast_db"]["card"])
        shell:
            "mkdir {params.path}; "
            "cd {params.path}; "
            "wget -O card_seq.tar.bz2 {params.seq}; "

rule unzip_card_db:
    input:
        seq="data/card_db/silva_seq.tar.bz2",
    output:
        seq=temp("data/card_db/rotein_fasta_protein_homolog_model.fasta"),
    shell:
        "tar -xvjf {input.seq}"

rule makeblastdb:
    input:
        seq1="data/card_db/rotein_fasta_protein_homolog_model.fasta"
        seq2="data/silva_db/silva_seq.fasta"
    output:
        output1="data/blast_db/card_db.pdb"
        output1="data/blast_db/silva_db.ndb"
    params:
        path=str(config["blast_db"]["blast"])
    shell:
        """
        makeblastdb -in silva_seq.fasta -dbtype nucl -out {params.path}/silva_db
        makeblastdb -in protein_fasta_protein_homolog_model.fasta -dbtype prot -out {params.path}/card_db
        """