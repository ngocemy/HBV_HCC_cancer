# Snakemake workflow to handle the RNAseq data analysis
# NOTE: Could add QC, trimming and differential expression analysis


R_IN = join(IN, 'rnaseq')
R_TMP = join(TMP, 'rnaseq')
R_OUT = join(OUT, 'rnaseq')

#"PHH_NI_3",
# 00: Generate genome index for HISAT2 RNAseq mapping
# NOTE: A wrapper for hisat2 index will be available when wrappers moves to 0.35
# It was added in 480ef74 but is not yet available (bitbucket.org/snakemake/snakemake-wrappers)
rule hisat2_index:
  input: 
    fasta = GENOME
  output:
    touch(join(R_TMP, 'indexing.done'))
  params:
    prefix = join(R_TMP, 'genome')
  threads: 12
  conda: "../envs/rnaseq_env.yaml"
  shell:
    """
    hisat2-build -p {threads} {input.fasta} {params.prefix}
    """
  #wrapper:
  #  "0.34.0/bio/hisat2/index"


# 01: Map RNAseq read using HISAT2
rule hisat2_align:
    input:
      idx_flag = join(R_TMP, 'indexing.done'),
      end1 = join(TMP, "reads", "{sample}_rnaseq_{rep}.end1.fq.gz")
    output:
      join(R_TMP, '{sample}_{rep}.bam')
    params:                             # idx is required, extra is optional
      idx = join(R_TMP, 'genome'),
      end2 = join(TMP, "reads", "{sample}_rnaseq_{rep}.end2.fq.gz")
    threads: 8
    conda: "../envs/rnaseq_env.yaml"
    shell:
      """
      hisat2 --threads {threads} \
             -1 {input.end1} \
             -2 {params.end2} \
             -x {params.idx} \
             --min-intronlen 1000 |
        samtools view -bh -o {output} -
      """
    #wrapper:
    #  "0.34.0/bio/hisat2/align"


# 02: Download hg19 transcriptome
rule get_transcriptome:
    output: join(R_IN, "hg19_transcriptome.fa.gz")
    shell:
        """
        ensembl="ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/"
        # Download all coding RNA transcripts
        wget -O {output} ${{ensembl}}/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        # Append noncoding RNA transcripts
        wget -O - ${{ensembl}}/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz >> {output}
        """


# 03: Combine BAM files from replicates and convert to BEDgraph of RNAseq coverage
rule bam_to_bedg:
    input: lambda w: expand(join(R_TMP, '{{sample}}_{rep}.bam'), rep=units.loc[units['sample'] == f'{w.sample}', 'replicate'])
    output: join(OUT, 'rnaseq', 'rna_coverage_{sample}.bedgraph')
    threads: 12
    conda: "../envs/rnaseq_env.yaml"
    shell:
        """
        samtools merge - {input} |
        samtools sort -@ {threads} | 
            bedtools genomecov -ibam - -bga > {output}
        """


# 04: Index transcriptome FASTA file to assign reads to it
rule salmon_index:
    input:
        join(R_IN, "hg19_transcriptome.fa.gz")
    output:
        directory(join(R_TMP, 'salmon', 'index'))
    threads: 12
    conda: "../envs/rnaseq_env.yaml"
    shell:
      """
      salmon index -p {threads} -i {output} -t {input}
      """
    #wrapper:
    #    "0.34.0/bio/salmon/index"


# 05: Quantifying transcript abundance. Note libtype=A tells salmon to 
# automatically determine whether the sample is stranded or unstranded
# Validatemappings runs extension alignment at the end to prevent spurious mappings
rule salmon_quant:
  input: 
    r1 = join(TMP, 'reads', '{sample}_rnaseq_{rep}.end1.fq.gz'),
    r2 = join(TMP, 'reads', '{sample}_rnaseq_{rep}.end2.fq.gz'),
    index = join(R_TMP, 'salmon', 'index')
  output: 
    quant = join(R_OUT, 'salmon', '{sample}', '{rep}', 'quant.sf'),
    lib = join(R_OUT, 'salmon', '{sample}', '{rep}', 'lib_format_counts.json')
  threads: 12
  conda: "../envs/rnaseq_env.yaml"
  params:
    libtype = "A",
    extra = '--validateMappings --rangeFactorizationBins 4'
  wrapper:
    "0.34.0/bio/salmon/quant"


rule diff_expr:
  input:
    expand(
      join(R_OUT, 'salmon', '{sample}', '{rep}', 'quant.sf'),
      sample=units.loc[units.libtype == "rnaseq", "sample"],
      rep=units.loc[units.libtype == "rnaseq", "replicate"]
    )
  output: join(R_OUT, 'diff_expr', 'integration_vs_control.tsv')
  params:
    salmon_dir = join(R_OUT, 'salmon'),
    samples = config['samples'],
    units = config['units']
  conda: "../envs/rnaseq_env.yaml"
  shell:
    """
    Rscript scripts/diff_expr.R {params.samples} \
                                {params.units} \
                                {params.salmon_dir} \
                                {output}
    """
