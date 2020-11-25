rule samtools_sort_3:
    input:
        "tmp/bam/{sample}.bam"
    output:
        bam = "tmp/bam/{sample}.sorted.bam"
    log:
        "tmp/logs/samtools_sort_3.{sample}.log"
    params:
        sort = "tmp/bam/{sample}.sorted.bam",
        scripts_path = config['scripts_dir'],
        whitelist = config['white_list'],
        collapsed =  "tmp/bam/{sample}.collapsed.bam",
        cb = "tmp/bam/{sample}.collapsed.bam.bam"
    run:
        shell("perl {params.scripts_path}/collapseBarcodesBam.pl {params.whitelist} {input} {params.collapsed}")
        shell("cp {params.cb} {input}")
        shell("samtools sort -o {output.bam} {input}")
        
rule samtools_index_4:
    input:
        "tmp/bam/{sample}.sorted.bam"
    output:
        "tmp/bam/{sample}.sorted.bam.bai"
    log:
        "tmp/logs/samtools_sort_3.{sample}.index.log"
    run:
        shell("samtools index {input} {output}")

