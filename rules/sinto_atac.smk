rule run_sinto_5:
    input:
        bam = "tmp/bam/{sample}.sorted.bam",
        index = rules.samtools_index_4.output
    output:
        out1 = 'tmp/fragments/{sample}.fragments',
        out2 = 'tmp/fragments/{sample}.fragements.sort.bed'
    log:
        "tmp/logs/run_sinto_5.{sample}.log"
    run:
        shell('sinto fragments -b {input.bam} -f {output.out1} --barcode_regex "[^:]*" 2> {log}')
        shell('sort -k 1,2 -k2,2n {output.out1} > {output.out2}')
