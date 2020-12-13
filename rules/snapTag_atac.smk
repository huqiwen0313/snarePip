rule snapTag_1:
    input:
        r1 = join(FASTQ_DIR, PATTERN_R1),
        r2 = join(FASTQ_DIR, PATTERN_R2),
        r3 = join(FASTQ_DIR, PATTERN_R3)

    output:
        fq1 = 'tmp/snapTag/{sample}.R1.fastq',
        fq2 = 'tmp/snapTag/{sample}.R3.fastq'

    params:
        prefix = '{sample}',
        dir = 'tmp/snapTag',
        script_dir = config['scripts_dir'],
        link_path = config['link_table']

    shell:
        "perl {params.script_dir}/get_barcodes_by_pos_SNARE2.pl {input.r1} {input.r2} {input.r3} {params.prefix} {params.dir} {params.link_path}"
