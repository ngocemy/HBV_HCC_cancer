#! usr/bin/env bash
HUMAN_FA='hbv_data/input/fa/bt2_hg19/hg19' # Index of bowtie2 indexing file for human hg19
VIRUS_FA='hbv_data/input/fa/bt2_hbv/hbv' # Index of bowtie2 indexing file for HBV
FASTQ_FILES="fq/hic/PM37_cap1.end1.fastq.gz, fq/hic/PM31_cap2.end1.fastq.gz, fq/hic/PM32_capture.rep1.end1.fastq.gz, fq/hic/PM25a.capture.rep1.end1.fastq.gz" #single end not paired end fastq, one end at a time
OUTDIR="polyidus/data/output/"
for f in "$FASTQ_FILES"; do
    mkdir "$OUTDIR/$(find "$f" | cut -d"/" -f 6)"
    python polyidus/src/polyidus.py "$HUMAN_FA" "$VIRUS_FA" --fastq "$f" --outdir "$OUTDIR/$(find "$f" | cut -d"/" -f 6)" 
done
