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
            seq="data/silva_db/card_seq.tar.bz2",
        params:
            seq=str(config["card"]["download-path-seq"]),
            tax=str(config["silva"]["download-path-tax"]),
            path=str(config["blast_db"]["silva"])
        shell:
            "mkdir {params.path}; "
            "cd {params.path}; "
            "wget -O silva_seq.fasta.gz {params.seq}; "
            "wget -O silva_tax.txt.gz {params.tax}; "

rule unzip_card_db:
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