rule B_FastQC_2:
    input:
        r1 = join(FASTQ_DIR, '{sample}_R1_filtered.fastq.gz')
    output:
        FQzip = 'tmp/FastQC/{sample}_R1_filtered_fastqc.zip'
    params:
        output_prefix = 'tmp/FastQC/',
        fastqc = config["fastqc_dir"]
    run:
        shell("{params.fastqc} -o '{params.output_prefix}' {input.r1}")
