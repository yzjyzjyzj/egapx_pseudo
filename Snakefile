from snakemake.io import glob_wildcards

SAMPLES = glob_wildcards("data/{sample}.gff3").sample

rule all:
    input:
        expand("results/{sample}.final.gff3", sample=SAMPLES),
        expand("results/{sample}.stats.txt", sample=SAMPLES)

rule extract_pseudogene_ids:
    input:
        gff="data/{sample}.gff3"
    output:
        ids="work/{sample}.pseudogene.ids"
    log:
        "logs/{sample}.extract_ids.log"
    shell:
        r"""
        grep -P '\t(pseudogene)\t' {input.gff} \
        | cut -f9 \
        | sed -n 's/.*ID=\([^;]*\).*/\1/p' \
        > {output.ids} 2>{log}
        """

rule extract_pseudogene_gff:
    input:
        gff="data/{sample}.gff3",
        ids="work/{sample}.pseudogene.ids"
    output:
        gff_pseudo="work/{sample}.pseudogene.gff3",
        gff_gene="work/{sample}.nonpseudo.gff3"
    log:
        "logs/{sample}.extract_gff.log"
    conda:
        "envs/agat.yaml"
    shell:
        """
        agat_sp_filter_feature_from_keep_list.pl \
        --gff {input.gff} --keep_list {input.ids} \
        --output {output.gff_pseudo} 2>{log}
        
        agat_sp_filter_feature_from_kill_list.pl --gff {input.gff} \
        --kill_list {input.ids} --output {output.gff_gene} 2>> {log}
        """

rule convert_gff_to_gtf:
    input:
        genome="data/{sample}.genome.fasta",
        gff3="work/{sample}.pseudogene.gff3"
    output:
        gtf="work/{sample}.pseudogene.gtf"
    log:
        "logs/{sample}.gff_to_gtf.log"
    conda:
        "envs/td2.yaml"
    shell:
        """
        scripts/gff3_gene_to_gtf_format.pl {input.gff3} {input.genome} \
        > {output.gtf} 2>{log}
        """

rule gtf_to_cdna:
    input:
        gtf="work/{sample}.pseudogene.gtf",
        genome="data/{sample}.genome.fasta"
    output:
        cdna="work/{sample}.cdna.fasta",
        gff3="work/{sample}.pseudogene.gtf.gff3"
    log:
        "logs/{sample}.gtf_to_cdna.log"
    conda:
        "envs/td2.yaml"
    shell:
        """
        scripts/gtf_genome_to_cdna_fasta.pl {input.gtf} {input.genome} \
        > {output.cdna} 2>{log}
        
        scripts/gtf_to_alignment_gff3.pl {input.gtf} > {output.gff3} 2>> {log}
        """


rule extract_orf:
    input:
        cdna="work/{sample}.cdna.fasta"
    output:
        outdir = directory("work/{sample}.cdna"),
        gff3="work/{sample}.cdna.fasta.TD2.gff3"
    log:
        "logs/{sample}.extract_orf.log"
    conda:
        "envs/td2.yaml"
    shell:
        """
        TD2.LongOrfs -t {input.cdna} --output-dir {output.outdir} 2>{log}
        
        TD2.Predict -t {input.cdna} -O {output.outdir} 2>>{log}
        
        mv {wildcards.sample}.cdna.fasta.TD2* work/
        """

rule generate_cds_gff:
    input:
        td2_gff3="work/{sample}.cdna.fasta.TD2.gff3",
        gtf_gff3="work/{sample}.pseudogene.gtf.gff3",
        cdna="work/{sample}.cdna.fasta"
    output:
        gff3="work/{sample}.cdna.fasta.TD2.genome.gff3"
    log:
        "logs/{sample}.generate_cds_gff.log"
    conda:
        "envs/td2.yaml"
    shell:
        """
        scripts/cdna_alignment_orf_to_genome_orf.pl \
        {input.td2_gff3} {input.gtf_gff3} {input.cdna} \
        > {output.gff3} 2>{log}
        """

rule filter_incomplete:
    input:
        gff3="work/{sample}.cdna.fasta.TD2.genome.gff3",
        genome="data/{sample}.genome.fasta",
        gff3_gene="work/{sample}.nonpseudo.gff3"
    output:
        gff3="work/{sample}.cdna.fasta.TD2.genome.filter.gff3",
        gff3_gene_filter="work/{sample}.nonpseudo.filter.gff3",
        gff3_combine="work/{sample}.combined.gff3",
        gff3_final="results/{sample}.final.gff3",
        gff3_stats="results/{sample}.stats.txt"
    log:
        "logs/{sample}.filter_incomplete.log"
    conda:
        "envs/agat.yaml"
    shell:
        """
        agat_sp_filter_incomplete_gene_coding_models.pl --gff {input.gff3} \
        --fasta {input.genome} -o {output.gff3} 2>{log}
        
        agat_sp_filter_incomplete_gene_coding_models.pl --gff {input.gff3_gene} \
        --fasta {input.genome} -o {output.gff3_gene_filter} 2>> {log}
        
        cat {output.gff3} {output.gff3_gene_filter} > {output.gff3_combine}
        
        bedtools sort -i {output.gff3_combine} > {output.gff3_final}
        
        agat_sp_statistics.pl --gff {output.gff3_final} -o {output.gff3_stats}  2>>{log}
        """


