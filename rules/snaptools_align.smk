rule snaptools_align_2:
    input:
        r1 = "tmp/snapTag/{sample}.R1.fastq",
        r2 = "tmp/snapTag/{sample}.R3.fastq"
    output:
        'tmp/bam/{sample}.bam'

    params:
        ref_fa = config['genome_mmi'],
        aligener_path = config['snap_aligner']

    log:
        "tmp/logs/snaptools_align_2.{sample}.log"

    shell:
        "snaptools align-paired-end --input-reference='{params.ref_fa}' "
        "--input-fastq1='{input.r1}' --input-fastq2='{input.r2}' "
        "--output-bam='{output}' --aligner=minimap2 --path-to-aligner='{params.aligener_path}' "
        "--read-fastq-command=cat --min-cov 0 --num-threads 6 --if-sort True --tmp-folder ./ --overwrite True"
        " 2> {log}"
