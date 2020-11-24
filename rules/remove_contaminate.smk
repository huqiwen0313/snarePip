rule A_remove_ATAC_contaminated_1:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2)
    output:
        r1_filtered = join(FASTQ_DIR, '{sample}_R1_filtered.fastq.gz'),
        r2_filtered = join(FASTQ_DIR, '{sample}_R2_filtered.fastq.gz')
    log:
        "tmp/logs/A_remove_ATAC_contaminated_1.{sample}.log"
    run:
        shell("cutadapt --no-indels --cores=2 --overlap 8 -e 0.2 -a ACGTACTGCAX --discard-trimmed -o {output.r2_filtered} -p {output.r1_filtered} {input.r2} {input.r1} > {log}")
