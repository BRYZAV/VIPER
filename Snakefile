import pandas as pd

#Load the configuration file for loading certain parameters
configfile: "config.yaml"

#Load the sample list into a table
samples = pd.read_table(config["samples"]).set_index("samples", drop=False)


rule all: #the target files
    input:
        expand("viral_matches/{sample}_viral_matches.tsv", sample=samples.index),
        # Aggregated stats produced by viral_search_stats
        "viral_stats/Viral_Species_counts.tsv",
        "viral_stats/Viral_Family_counts.tsv"


rule fastp_filter:
    input:
        config["fastq_dir"]+"/{sample}.fastq.gz"
    output:
        trim="trimmed_reads/{sample}_trimmed.fastq.gz",
        html="fastp_results/{sample}.html",
        json="fastp_results/{sample}.json"
    threads: 12 
    shell: 
        "fastp -i {input} -e 20 -l 50 -o {output.trim} -j {output.json} -h {output.html} -w {threads}"


rule bowtie2_map:
    input: 
        fastq="trimmed_reads/{sample}_trimmed.fastq.gz",
    output:
        gz="host_filtered_data/unmapped_reads/unmapped_{sample}.fastq.gz",
        sam="host_filtered_data/unmapped_reads/sam_files/{sample}_aln.sam"
    log:
        "logs/bowtie2/{sample}.log"
    params:
        bowtie_index=config["bowtie_index"]
    threads: 24 
    shell:
        "bowtie2 -x {params.bowtie_index} -U {input.fastq} -p {threads} --un-gz {output.gz} "
        "-S {output.sam} 2> {log}"


rule megahit_assembly:
    input:
        "host_filtered_data/unmapped_reads/unmapped_{sample}.fastq.gz"
    output:
        directory("megahit_results/{sample}_contigs")
    threads: 24 
    shell:
        "megahit -r {input} -o {output} -t {threads}" 


rule bbwrap_map:
    input:
        unmap="host_filtered_data/unmapped_reads/unmapped_{sample}.fastq.gz",
        contigs="megahit_results/{sample}_contigs"
    output:
        sam="host_filtered_data/singleton/sam_files/{sample}_singleton_aln.sam.gz"
    threads: 24 
    shell:
        "bbwrap.sh ref={input.contigs}/final.contigs.fa in={input.unmap} out={output.sam} path={input.contigs}/bbwrap_index threads={threads}"


rule fastq_convert:
    input:
	    "host_filtered_data/singleton/sam_files/{sample}_singleton_aln.sam.gz"
    output:
	    "host_filtered_data/singleton/{sample}_singleton.fasta"
    shell:
	    "samtools fastq -f 4 {input} | seqkit fq2fa -o {output}"



rule contig_diamond_search:
    input:
        "megahit_results/{sample}_contigs"
    output:
        "diamond_matches/contig/{sample}_dmnd_matches.tsv"
    log:
        "logs/diamond/contig/{sample}_diamond.log"
    threads: 32 
    shell:
        "diamond blastx -q {input}/final.contigs.fa -d /90daydata/zedru/nr_db.dmnd -k 1 -p {threads} "
        "--outfmt 6 qseqid sseqid pident evalue qcovhsp bitscore length mismatch gapopen qstart qend qlen qframe "
        "sstart send slen sscinames skingdoms sskingdoms sphylums staxids stitle slineages "
        "--header -o {output} 2> {log}"



rule singleton_diamond_search:
    input:
        "host_filtered_data/singleton/{sample}_singleton.fasta"
    output:
        "diamond_matches/singleton/{sample}_dmnd_matches.tsv"
    log:
        "logs/diamond/singleton/{sample}_diamond.log"
    threads: 32 
    shell:
        "diamond blastx -q {input} -d /90daydata/zedru/nr_db.dmnd -k 1 -p {threads} "
        "--outfmt 6 qseqid sseqid pident evalue qcovhsp bitscore length mismatch gapopen qstart qend qlen qframe "
        "sstart send slen sscinames skingdoms sskingdoms sphylums staxids stitle slineages "
        "--header -o {output} 2> {log}"


rule viral_search:
    input:
        contig="diamond_matches/contig/{sample}_dmnd_matches.tsv",
        singleton="diamond_matches/singleton/{sample}_dmnd_matches.tsv"
    output:
        "viral_matches/{sample}_viral_matches.tsv"
    shell:
        "cat {input.contig} {input.singleton} | "
        "grep 'Viruses' | "
        "grep -Eiv 'hypothetical|phage|Escherichia|Caudoviricetes|Pseudomonadota|Uroviricota' | "
        "csvtk add-header -t -n qseqid,sseqid,pident,evalue,qcovhsp,bitscore,length,mismatch,gapopen,qstart,qend,qlen,qframe,sstart,"
        "send,slen,sscinames,skingdoms,sskingdoms,sphylums,staxids,stitle,slineages > {output}"

rule viral_search_stats:
    input:
        # all per-sample viral match TSVs
        expand("viral_matches/{sample}_viral_matches.tsv", sample=samples.index)
    output:
        species="viral_stats/Viral_Species_counts.tsv",
        family="viral_stats/Viral_Family_counts.tsv"
    shell:
        "cat {input} | csvtk freq -t -f sscinames -n -r | csvtk rename -t -f frequency -n Sequence_Counts > {output.species} && "
        "cut -f 21 {input}/*tsv | tail -n +2 | taxonkit reformat -I 1 -f {{f}} | csvtk -t add-header -n staxids,Family "
        "| csvtk freq -t -f Family -n -r | csvtk rename -t -f frequency -n Sequence_Counts > {output.family}"
