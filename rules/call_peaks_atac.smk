rule call_peaks_6:
    input:
        bams = "tmp/bam/{sample}.sorted.bam",
        snap = "tmp/snap/{sample}.snap",
        index = rules.samtools_index_4.output
    output:
        bams = "tmp/bam/{sample}.rmsk.bam",
        peaks = "tmp/macs/{sample}.peaks.bed"
    params:
        blacklist = config['blacklist'],
        bedtools = config['bedtools'],
        macs = config['macs'],
        prefix = '{sample}',
        macs_dir = "tmp/macs",
        macs_peaks_output = "tmp/macs/{sample}_peaks.broadPeak"
    log:
        "logs/call_peaks_6.{sample}.log"
    run:
        shell("{params.bedtools} intersect -v -abam {input.bams} -b {params.blacklist} > {output.bams}")
        shell("{params.macs} callpeak -t {output.bams} -g hs --keep-dup all --nomodel --extsize 500 -q 0.1 -n {params.prefix} --outdir {params.macs_dir} --tempdir {params.macs_dir} --broad --broad-cutoff 0.1 2> {log}")
        shell("cut -f1-3 {params.macs_peaks_output} > {output.peaks}")
        shell("snaptools snap-add-pmat --snap-file {input.snap} --peak-file {output.peaks}")
