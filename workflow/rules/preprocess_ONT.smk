rule decompress_fastq:
    input:
        "{base_dir}/data/input/{sample}.fastq.gz"
    output:
        "{base_dir}/data/input/{sample}.fastq"
    log:
        "{base_dir}/logs/decompress_fastq/{sample}.log"    
    shell:
        "gzip -dk {input} 2> {log}"

rule convert_fastq_to_fasta:
    input:
        "{base_dir}/data/input/{sample}.fastq"
    output:
        "{base_dir}/data/input/temp/{sample}.fasta"
    log:
        "{base_dir}/logs/convert_fastq_to_fasta/{sample}.log"    
    conda:
        "../envs/python.yaml"        
    shell:
        "seqtk seq -a {input} > {output} 2> {log}"

rule split_fasta_file:
    input:
        "{base_dir}/data/input/temp/{sample}.fasta"
    output:
        "{base_dir}/data/input/{sample}.part_{part}.fasta"
    params:
        outdir = config["base_dir"],
        num_parts = config["num_parts"],    
    shell:
        """
        seqkit split2 --by-part {params.num_parts} {input} --out-dir {params.outdir}/data/input/
        """