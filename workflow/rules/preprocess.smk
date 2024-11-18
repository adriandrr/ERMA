if config["seq_tech"] == "Illumina":
    rule merge_fastq:
        input:
            r1="{base_dir}/data/fastq/{sample}_R1_001.fastq.gz",
            r2="{base_dir}/data/fastq/{sample}_R2_001.fastq.gz"
        output:
            out = temp("{base_dir}/data/fastq/{sample}.fastq.gz")
        log:
            "{base_dir}/logs/merge_fastq/{sample}.log"
        conda:
            "../envs/python.yaml"                
        shell:
            "NGmerge -1 {input.r1} -2 {input.r2} -o {output.out}"

rule decompress_fastq:
    input:
        "{base_dir}/data/fastq/{sample}.fastq.gz"
    log:
        "{base_dir}/logs/decompress_fastq/{sample}.log"    
    output:
        temp("{base_dir}/data/fastq/{sample}.fastq")
    shell:
        "gzip -dk {input} 2> {log}"

rule convert_fastq_to_fasta:
    input:
        "{base_dir}/data/fastq/{sample}.fastq"
    output:
        temp("{base_dir}/data/fastq/{sample}.fasta")
    log:
        "{base_dir}/logs/convert_fastq_to_fasta/{sample}.log"    
    conda:
        "../envs/python.yaml"        
    shell:
        "seqtk seq -a {input} > {output} 2> {log}"

rule split_fasta_file:
    input:
        "{base_dir}/data/fastq/{sample}.fasta"
    output:
        temp("{base_dir}/data/fastq/{sample}.part_{part}.fasta")
    params:
        outdir = config["base_dir"],
        num_parts = config["num_parts"],    
    run:
      for sample in samples:
            fasta_file = f"{params.outdir}/data/fastq/{sample}.fasta"
            shell("seqkit split2 --by-part {params.num_parts} {fasta_file} --out-dir {params.outdir}/data/fastq/")

rule gzip_blast_results:
    input:
        silva_results = "{base_dir}/results/{sample}/{part}/SILVA_results.txt",
        card_results = "{base_dir}/results/{sample}/{part}/card_results.txt"
    output:
        silva_zip = "{base_dir}/results/{sample}/{part}/SILVA_results.txt.gz",
        card_zip = "{base_dir}/results/{sample}/{part}/card_results.txt.gz"
    log:
        "{base_dir}/logs/gzip_blast_results/{sample}_{part}.log"        
    shell:
        """
        gzip {input.silva_results} 2> {log}
        gzip {input.card_results} 2>> {log}
        """

rule gzip_filtered_blast_data:
    input:
        int_data = "{base_dir}/results/{sample}/{part}/filtered_results.csv"
    log:
        "{base_dir}/logs/gzip_integrated_blast_data/{sample}_{part}.log"    
    output:
        int_data_zip = "{base_dir}/results/{sample}/{part}/filtered_results.csv.gz"
    shell:
        "gzip -9 {input.int_data} 2> {log}"