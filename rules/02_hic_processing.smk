!# usr/bin/env python3
conda: "envs/hic_processing.yaml"
rule bt2_index:
  input: GENOME
  output: touch(join(TMP, 'index_genome.done'))
  params:
    idx = join(TMP, 'genome')
  singularity: config['containers']['hicstuff']
  shell: "bowtie2-build {input} {params.idx}"

# Make splits from Hi-C fastq files to speed up mapping. [between 1 and 999]
N_SPLITS = 100
split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)] #base names of split files

rule split_hic_fastq:
  input: join(TMP, 'reads', '{sample}_{libtype}_1.{end}.fq.gz')
  output: 
    expand(
      join(
        TMP,
        'split_reads',
        '{{sample}}_{{libtype}}_1_{{end}}',
        '{{sample}}_{{libtype}}_1.{{end}}.{split}.fq.gz'
      ),
      split=split_names
    )
  params:
    n_splits = N_SPLITS,
    split_dir = lambda w: join(TMP, 'split_reads', f"{w.sample}_{w.libtype}_1_{w.end}")
  message: "Splitting {wildcards.sample}_{wildcards.end} into {params.n_splits} split fastq"
  shell:
     """
     mkdir -p {params.split_dir}
     # 100 split fastqs will be created with name pattern 00000.fq - 000100.fq
     seqkit split2 -p {params.n_splits} \
                   -w 0 \
                   -f \
                   -1 {input} \
                   -O {params.split_dir}
     """

# Iterative alignment of a single fastq split from a Hi-C sample
rule split_iter_align_hic:
  input:
    index_flag = join(TMP, 'index_genome.done'),
    fq = join(
      TMP,
      'split_reads',
      '{sample}_{libtype}_1_{end}',
      '{sample}_{libtype}_1.{end}.{split}.fq.gz'
    ),
  output:
    join(
      TMP,
      'split_reads',
      '{sample}_{libtype}_1_{end}',
      '{sample}_{libtype}_1.{end}.{split}.bam'
    )
  params:
    tmp_dir = lambda w: join(TMP, "split_reads", f"{w.sample}_hic_{w.end}", f"{w.split}"),
    index = join(TMP, 'genome'),
    iteralign_presets = config['params']['hicstuff_iteralign']
  threads: 12
  singularity: config['containers']['hicstuff']
  shell:
    """
    hicstuff iteralign {params.iteralign_presets} \
                       -t {threads} \
                       -T {params.tmp_dir} \
                       -g {params.index} \
                       -o {output} \
                       {input.fq}
    """

# Merge splits from individual mapping jobs into one bam file per library end
rule merge_split_alignments:
  input:
    expand(
      join(
        TMP,
        'split_reads',
        '{{sample}}_{{libtype}}_1_{{end}}',
        '{{sample}}_{{libtype}}_1.{{end}}.{split}.bam'),
      split=split_names
    )
  output: join(TMP, 'bam', '{sample}_{libtype}_1.{end}.bam')
  threads: 12
  singularity: config['containers']['hicstuff']
  shell: "samtools merge -n -O BAM -@ {threads} {output} {input}"

  
## 00 Generate Hi-C pairs files
rule generate_pairs:
  input:
    bam1 = join(TMP, 'bam', "{sample}_{libtype}_1.end1.bam"),
    bam2 = join(TMP, 'bam', "{sample}_{libtype}_1.end2.bam"),
    index_flag = join(TMP, 'index_genome.done')
  output:
    hicdir = temporary(directory(join(TMP, 'hicstuff', '{sample}_{libtype}'))),
    pairs = join(TMP, 'pairs', '{sample}_{libtype}.pairs')
  params:
    res = config['contact_maps']['max_res'],
    idx = join(TMP, 'genome')
  threads: 12
  singularity: config['containers']['hicstuff']
  shell:
    """
    hicstuff pipeline -e {params.res} \
                      -g {params.idx} \
                      -o {output.hicdir} \
                      -S bam \
                      -t {threads} \
                      -nD \
                      {input.bam1} \
                      {input.bam2}

    cp {output.hicdir}/tmp/valid_idx_pcrfree.pairs {output.pairs}
    """

# 00: Generate chrom sizes file
rule chrom_sizes:
  input: GENOME
  output: join(TMP, "chrom.sizes")
  params:
    genome = GENOME,
    res = MAX_RES,
    digest_dir = join(TMP, 'digest')
  singularity: config['containers']['hicstuff']

  shell:
    """
    hicstuff digest -e {params.res} \
                    -o {params.digest_dir} \
                    {input}

    tail -n +2 {params.digest_dir}/info_contigs.txt \
      | awk -vOFS='\t' '{{print $1,$2}}' \
      > {output}
    """


rule compress_pairs:
  input: join(TMP, 'pairs', '{sample}_{libtype}.pairs')
  output: join(TMP, 'pairs', '{sample}_{libtype}.pairs.gz')
  shell: "gzip {input}"

# 01: Convert pairs files to cool format and normalize it.
# Generates a high res and a low res format for each pairs file.
rule pairs_to_cool:
  input: 
    pairs = join(TMP, 'pairs', '{sample}_{libtype}.pairs.gz'),
    chromsizes = join(TMP, 'chrom.sizes')
  output:
    cool_max = join(OUT, 'cool', '{sample}_{libtype}_' + f'{MAX_RES_STR}.cool'),
    cool_comp = join(OUT, 'cool', '{sample}_{libtype}_' + f'{COMP_RES_STR}.cool')
  threads: 3
  params:
    max_res = MAX_RES,
    comp_res = COMP_RES
  threads: 12
  singularity: config['containers']['cooler']
  shell: 
    """
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
                       {input.chromsizes}:{params.max_res} \
                       {input.pairs} {output.cool_max}
    cooler balance -p {threads} --mad-max 10 {output.cool_max}
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
                       {input.chromsizes}:{params.comp_res} \
                       {input.pairs} {output.cool_comp}
    cooler balance -p {threads} --mad-max 10 {output.cool_comp}
    """

# 02: Add chromosome information to cool files in place
rule add_chrom_arms:
  input:
    cool = join(OUT, 'cool', '{sample}_{libtype}_{res}.cool'),
    centro = join(IN, 'centro_hg19.bed')
  output: temp(touch(join(TMP, '{sample}_{libtype}_{res}_chrom_arm.done')))
  run:
    add_chromosome_arms(input['cool'], input['centro'])

# 03a: Compute contacts with virus along the genome to detect insertion events
rule virus_coverage:
  input: join(TMP, 'pairs', '{sample}_captn.pairs')
  output: join(OUT, 'cov_virus_{sample}.bedgraph')
  threads: 3
  run:
    chrom_contact_to_bedgraph(
      input[0], 
      output[0],
      "HBVayw",
      binsize=MAX_RES,
      threads=threads,
      buffer="4G"
    )

# 03c: Compute insulation scores along the matrix (TAD boundaries) the diamond 
# window is set to 50x50 pixels
rule insulation_score:
  input: join(OUT, 'cool', '{sample}_hic_' + f'{COMP_RES_STR}.cool')
  output: join(OUT, 'insulation_{sample}.bedgraph')
  params:
    win_size_bp = 10 * COMP_RES
  singularity: config['containers']['cooler']
  shell: 
    """
    cooltools diamond-insulation {input} {params.win_size_bp} |
      tail -n +2 |
      awk -vOFS="\t" '{{print $1,$2,$3,$5}}' > {output}
    """
# 04: Extracting read names of human-virus contacts
# those reads are not necessary to be chimeric, but based on 3D contacts human-virus
rule hetero_reads:
    input: join(TMP,'pairs', '{sample}_{libtype}.pairs')
    output: join(TMP, 'hetero_reads_{sample}_{libtype}.txt')
    threads: 4
    script:
        "extract_hetero_reads.py"
