rule C_dropTag_3:
    input:
        r1 = join(FASTQ_DIR, '{sample}_R1_filtered.fastq.gz'),
        r2 = join(FASTQ_DIR, '{sample}_R2_filtered.fastq.gz'),
        fastqc = 'tmp/FastQC/{sample}_R1_filtered_fastqc.zip'
    output:
        'tmp/tagged/{sample}.fastq.gz'

    params:
        config_file = config['dropest_config'],
        prefix = 'tmp/tagged/{sample}',
        log_prefix = 'tmp/logs/B_dropTag_2.{sample}'
    log:
        'tmp/logs/B_dropTag_2.{sample}_tag_main.log'

    shell:
        "droptag -p 2 -r 0 -c {params.config_file} -n {params.prefix} -l {params.log_prefix} "
        "{input.r2} {input.r1} > {log}"

rule C_alignment_4:
    input:
        'tmp/tagged/{sample}.fastq.gz'
    output:
        bam = 'tmp/alignment/{sample}.Aligned.out.bam',
        log_STAR = 'tmp/alignment/{sample}.Log.final.out'
    log:
        'tmp/alignment/{sample}.Log.out'
    params:
        star_index = config['index_file'],
        prefix = 'tmp/alignment/{sample}'
    shell:
        "STAR --runThreadN 2 "
        "--genomeDir '{params.star_index}' --limitOutSJcollapsed 2000000 "  # max number of collapsed junctions
        "--readFilesCommand zcat --outSAMtype BAM Unsorted "
        "--readFilesIn '{input}' "
        "--outFileNamePrefix '{params.prefix}'."

rule D_dropEst_5:
    input:
        bam = 'tmp/alignment/{sample}.Aligned.out.bam'
    output:
        cells_tsv = 'tmp/dropest_out/{sample}.cells.tsv',
        genes_tsv = 'tmp/dropest_out/{sample}.genes.tsv',
        mtx = 'tmp/dropest_out/{sample}.mtx',
        rds = 'tmp/dropest_out/{sample}.rds'
    params:
        config_fileDest = config['dropest_config'],
        gtf_file = config['gtf_dir'],
        output_prefix = "tmp/dropest_out/{sample}",
        log_prefix = 'tmp/logs/D_dropEst_4.{sample}'
    log:
        'tmp/logs/D_dropEst_4.{sample}_est_main.log'
    run:
        shell("if [ ! -d tmp/dropest_out ]; then mkdir tmp/dropest_out; fi")
        shell("dropest -w -M -u -G 1 -L iIeEBA -u -m -F -l {params.log_prefix} "
        "-g '{params.gtf_file}' -c '{params.config_fileDest}' "
        "-o '{params.output_prefix}' "
        "'{input.bam}' > {log}")


