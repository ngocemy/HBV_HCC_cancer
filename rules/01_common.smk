#!/bin/env snakemake -s


def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    targets = units.loc[
      (units['sample'] == wildcards.sample) &
      (units['libtype'] == wildcards.libtype) & 
      (units['replicate'] == wildcards.rep)
      ,
      f'fq{wildcards.end}'
    ]
    return targets.tolist()


# Combine all fastq files from the same sample / library type combination
rule combine_units:
  input: get_fastqs
  message:
    """
	Running combine_units
	Wildcards: {wildcards}
    Merging :
      {input} into:  {output}
    """
  output:
    join(TMP, "reads", "{sample}_{libtype}_{rep}.end{end}.fq.gz")
  threads: 3
  shell: "cat {input} > {output}"

rule merge_reps:
    input:
      lambda w: expand(
        join(TMP, "reads", "{sample}_{libtype}_{rep}.end{end}.fq.gz"),
        sample=w['sample'],
		libtype=w['libtype'],
        rep=units.loc[(units['sample'] == w['sample']) & (units.libtype == w['libtype']), 'replicate'],
        end=w['end']
    )
    output: temporary(join(TMP, '{libtype}_merge','{sample}_{libtype}.end{end}.fq.gz'))
    message:
      """
	  Running merge_reps
      Wildcard is: {wildcards}
      Merging :
        {input} into:  {output}
      """
    shell: "cat {input} > {output}"
