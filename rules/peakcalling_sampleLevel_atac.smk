rule sample_level_peaks_10:
    input:
        p2 = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/objects/{sample_name}.p2.obj.rds"),
    output:
        bed = join(ASSAY, '{tissue}/samples/{sample_name}/Sample_output/macs2/{sample_name}.combined.peaks.bed'),
    params:
        sampleID = '{sample_name}',
        path_to_macs = config['macs'],
        atacPath = CURRENT_DIR,
        rnaPath = config['rna_path'],
        gsize = "hs",
        scripts_path = config['scripts_dir'],
        assay_type = config['assay_type']
    shell:
        "Rscript {params.scripts_path}/sampleLevel.peakCalling.R {params.atacPath} {params.rnaPath} {params.sampleID} {params.path_to_macs} {params.gsize} {params.assay_type}"

rule merge_peaks_11:
    input:
        peaks = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/macs2/{sample_name}.combined.peaks.bed")
    output:
        merged_peaks = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/macs2/{sample_name}.combined.peaks.merged.bed")
    log:
        "tmp/logs/merge_peaks_11.{sample_name}.{tissue}.log"
    params:
        sort_peaks = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/macs2/{sample_name}.combined.peaks.sorted.bed"),
        experiment_peaks = join(ASSAY, "{tissue}/samples/{sample_name}/Experiment_output/macs/*.peaks.bed"),
        merged_peaks = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/macs2/{sample_name}.peaks.bed")
    run:
        shell("cat {params.experiment_peaks} {input.peaks} > {params.merged_peaks}")
        shell("sort -k1,1 -k2,2n {params.merged_peaks} > {params.sort_peaks}")
        shell("bedtools merge -i {params.sort_peaks} -d 300 > {output.merged_peaks}")

rule sample_snap_pre_12:
    input:
        peak = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/macs2/{sample_name}.combined.peaks.merged.bed")
    output:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/snap/{sample_name}.snap")
    params:
        genome_name = config['genome_name'],
        genome_size = config['hs_csize'],
        bam = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/bam/{sample_name}.bam")
    run:
        shell("samtools sort -n -@ 15 -m 50000000000 -o {params.bam}.sorted.bam {params.bam}")
        shell("snaptools snap-pre --input-file='{params.bam}.sorted.bam' "
        "--output-snap='{output}' --genome-name='{params.genome_name}' "
        "--genome-size='{params.genome_size}' --min-mapq=30 --min-flen=0 --max-flen=1000 --keep-chrm=TRUE "
        "--keep-single=FALSE  --keep-secondary=FALSE --overwrite=True --min-cov=100 --verbose=True ")
        shell("snaptools snap-add-bmat --snap-file {output} --bin-size-list 1000 5000 10000 --verbose TRUE")
        shell("snaptools snap-add-pmat --snap-file {output} --peak-file {input.peak}")
