import pandas as pd

#Load the configuration file for loading certain parameters
configfile: "config.yaml"

#Load the sample list into a table
samples = pd.read_table(config["samples"]).set_index("samples", drop=False)


rule all: #the target files
    input:
        expand("viral_matches/{sample}_viral_matches.tsv", sample=samples.index),
        # Create Fastqc filesa
        expand("fastqc_results/{sample}_fastqc.html", sample=samples.index),
        expand("fastqc_results/{sample}_fastqc.zip", sample=samples.index),
        # Aggregated stats produced by viral_search_stats
        "viral_stats/Viral_Species_counts.tsv",
        "viral_stats/Viral_Family_counts.tsv",
        #MultiQC Report
        "MultiQC_report/MultiQC_report.html"


rule fastp_filter:
    input:
        config["fastq_dir"]+"/{sample}.fastq.gz"
    output:
        trim="trimmed_reads/{sample}_trimmed.fastq.gz",
        html="fastp_results/{sample}.html",
        json="fastp_results/{sample}.json"
    params:
        min_read_length=config["min_read_length"],
        min_mean_qual=config["min_mean_qual"]  
    threads: 
        config["fastp_threads"]
    shell: 
        "fastp -i {input} -e {params.min_mean_qual} -l {params.min_read_length} -o {output.trim} -j {output.json} -h {output.html} -w {threads}"


rule fastqc_results:
    input:
        "trimmed_reads/{sample}_trimmed.fastq.gz"
    output:
        html="fastqc_results/{sample}_fastqc.html",
        zip="fastqc_results/{sample}_fastqc.zip"
    threads:
        config["fastqc_threads"]
    shell:
        "mkdir -p fastqc_results/ && "
        "fastqc {input} -t {threads} -o fastqc_results/ && "
        "mv fastqc_results/{wildcards.sample}_trimmed_fastqc.html fastqc_results/{wildcards.sample}_fastqc.html && "
        "mv fastqc_results/{wildcards.sample}_trimmed_fastqc.zip fastqc_results/{wildcards.sample}_fastqc.zip"


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
    threads: 
        config["bowtie_threads"]
    shell:
        "bowtie2 -x {params.bowtie_index} -U {input.fastq} -p {threads} --un-gz {output.gz} "
        "-S {output.sam} 2> {log}"


rule megahit_assembly:
    input:
        "host_filtered_data/unmapped_reads/unmapped_{sample}.fastq.gz"
    output:
        directory("megahit_results/{sample}_contigs")
    log:
        "logs/megahit/{sample}.log"
    threads:
        config["megahit_threads"]
    shell:
        "megahit -r {input} -o {output} --out-prefix {wildcards.sample} -t {threads} 2> {log}" 


rule bbwrap_map:
    input:
        unmap="host_filtered_data/unmapped_reads/unmapped_{sample}.fastq.gz",
        contigs="megahit_results/{sample}_contigs"
    output:
        sam="host_filtered_data/singleton/sam_files/{sample}_singleton_aln.sam.gz",
        fastq="host_filtered_data/singleton/{sample}_singleton.fastq"
    threads:
        config["bbwrap_threads"]
    shell:
        "bbwrap.sh ref={input.contigs}/{wildcards.sample}.contigs.fa in={input.unmap} out={output.sam} path={input.contigs}/bbwrap_index threads={threads} && "
        "samtools fastq -f 4 {output.sam} > {output.fastq}"


rule contig_diamond_search:
    input:
        "megahit_results/{sample}_contigs"
    output:
        "diamond_matches/contig/{sample}_dmnd_matches.tsv"
    log:
        "logs/diamond/contig/{sample}_diamond.log"
    threads:
        config["diamond_threads"]
    shell:
        "diamond blastx -q {input}/{wildcards.sample}.contigs.fa -d /90daydata/zedru/nr_db.dmnd -k 1 -p {threads} "
        "--outfmt 6 qseqid sseqid pident evalue qcovhsp bitscore length mismatch gapopen qstart qend qlen qframe "
        "sstart send slen sscinames skingdoms sskingdoms sphylums staxids stitle slineages "
        "--header -o {output} 2> {log}"



rule singleton_diamond_search:
    input:
        "host_filtered_data/singleton/{sample}_singleton.fastq"
    output:
        "diamond_matches/singleton/{sample}_dmnd_matches.tsv"
    log:
        "logs/diamond/singleton/{sample}_diamond.log"
    threads:
        config["diamond_threads"]
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
        counts=temp("viral_stats/Viral_counts.tsv"),
        taxids=temp("viral_stats/Viral_Taxonomy.tsv"),
        species="viral_stats/Viral_Species_counts.tsv",
        family="viral_stats/Viral_Family_counts.tsv"
    shell:
        "cat {input} | csvtk freq -t -f sscinames -n -r | csvtk rename -t -f frequency -n Sequence_Counts > {output.counts} && "
        "cat {input} | csvtk freq -t -f staxids -n -r | cut -f 1 | tail -n +2 | taxonkit reformat -I 1 -f {{f}} | csvtk -t add-header -n staxids,Family > {output.taxids} && "
        "paste {output.counts} {output.taxids} > {output.species} && "
        "cut -f 21 {input} | tail -n +2 | taxonkit reformat -I 1 -f {{f}} | csvtk -t add-header -n staxids,Family "
        "| csvtk freq -t -f Family -n -r | csvtk rename -t -f frequency -n Sequence_Counts > {output.family}"

rule multiqc_report:
    input:
        expand("fastqc_results/{sample}_fastqc.zip", sample=samples.index),
        expand("fastp_results/{sample}.json", sample=samples.index),
        expand("logs/bowtie2/{sample}.log", sample=samples.index),
        expand("logs/megahit/{sample}.log", sample=samples.index)
    output:
        html="MultiQC_report/MultiQC_report.html"
    shell:
        "multiqc {input} -n MultiQC_report -o MultiQC_report/"

