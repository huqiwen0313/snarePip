rule snap_pre_5:
    input:
        bam = "tmp/bam/{sample}.bam",
        index = rules.samtools_index_4.output
    output:
        "tmp/snap/{sample}.snap"
    log:
        "tmp/logs/snap_pre_5.{sample}.log"
    params:
        genome_name = config['genome_name'],
        genome_size = config['hs_csize']
    run:
        shell("snaptools snap-pre --input-file='{input.bam}' "
        "--output-snap='{output}' --genome-name='{params.genome_name}' "
        "--genome-size='{params.genome_size}' --min-mapq=30 --min-flen=0 --max-flen=1000 --keep-chrm=TRUE "
        "--keep-single=FALSE  --keep-secondary=FALSE --overwrite=True --min-cov=100 --verbose=True 2> {log}")
        shell("snaptools snap-add-bmat --snap-file {output} --bin-size-list 1000 5000 10000 --verbose TRUE")
