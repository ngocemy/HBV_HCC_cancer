
# Download Hint reference file, bwa index and background matrix from their dropbox...
rule hint_download:
  output:
    bwa_index = directory(join(TMP, 'hint', 'bwa_index', 'hg19')),
    bg_matrix = directory(join(TMP, 'hint', 'bg_matrix', 'hg19')),
    ref = directory(join(TMP, 'hint', 'ref', 'hg19'))
  params:
    base_url = "https://www.dropbox.com/sh"
  shell:
    """
    wget -O "{output.bwa_index}.zip" "{params.base_url}/2ufsyu4wvrboxxp/AAB9zR4BlbpwBnVabfYFPfkZa/bwaIndex/hg19.zip?dl=1"
    unzip -d "$(dirname {output.bwa_index})" "{output.bwa_index}.zip"
    wget -O "{output.bg_matrix}.zip" "{params.base_url}/2ufsyu4wvrboxxp/AAD1kpUdm7Zgu6EM8dss6BjRa/backgroundMatrices/hg19.zip?dl=1"
    unzip -d "$(dirname {output.bg_matrix})" "{output.bg_matrix}.zip"
    wget -O "{output.ref}.zip" "{params.base_url}/qas48d7409t2syz/AAABcBBuR4Nhgk0Gpnb_uhdDa/hg19.zip?dl=1"
    unzip -d "$(dirname {output.ref})" "{output.ref}.zip"
    """

rule make_100kb_matrix:
  input: join(OUT, 'cool', '{library}.cool')
  output: join(TMP, 'hint', 'cool100kb', '{library}_100kb.cool')
  threads : 12
  shell:
    """
    cooler coarsen \
      -k 10 \
      -n {threads} \
      -o {output} \
      {input}
    cooler balance -p {threads} {output}
    """

# Detect translocations using hint
rule hin_tl:
  input:
    cool10 = join(OUT, 'cool', '{library}.mcool'),
    cool100 = join(TMP, 'hint', 'cool100kb', '{library}_100kb.cool'),
    bwa_index = join(TMP, 'hint', 'bwa_index', 'hg19'),
    bg_matrix = join(TMP, 'hint', 'bg_matrix', 'hg19'),
    ref = join(TMP, 'hint', 'ref', 'hg19')
  params:
    res = MAX_RES,
    lib = lambda w: f"{w.library}"
  conda: "../envs/hic_processing.yaml"
  output: directory(join(OUT, 'hint', '{library}'))
  shell:
    """
    hint tl \
      -m {input.cool10}::/resolutions/{params.res},{input.cool100} \
      -f cooler \
      --refdir {input.ref} \
      --backdir {input.bg_matrix} \
      -g hg19 \
      -n {params.lib} \
      -c 0.05 \
      --ppath $(which pairix) \
      -p 12\
      -o {output}
    """

