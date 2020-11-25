rule fastqc_7:
    input:
        bam = "tmp/bam/{sample}.bam"
    output:
        fastqcout = "tmp/fastqFiles/{sample}_fastqc.zip"
    params:
       fastqcout = "tmp/fastqFiles",
       fastqc_dir = config['fastqc_dir']
    run:
        shell("if [ ! -d tmp/fastqFiles ]; then mkdir tmp/fastqFiles; fi")
        shell("{params.fastqc_dir} -o {params.fastqcout} -f bam {input.bam}")

rule report_8:
    input:
        fastqcout = "tmp/fastqFiles/{sample}_fastqc.zip",
        peaks = "tmp/macs/{sample}.peaks.bed"
    output:
        rmdout = "tmp/reports/{sample}.qc.report.html"
    params:
        fastqc="tmp/fastqFiles",
        prefix = '{sample}',
        outdir = CURRENT_DIR + "/tmp",
        rmd = config['scripts_dir'] + "/ATAC.report.experiment.level.html.rmd",
        rmdout = CURRENT_DIR + "/tmp/reports/{sample}.qc.report.html",
        RNAdir = config['rna_path'],
        scripts_dir = config['scripts_dir']
    log:
        "tmp/logs/report_8.{sample}.log"
    run:
        shell("if [ ! -d tmp/QCs ]; then mkdir tmp/QCs; fi")
        shell("if [ ! -d tmp/pmats ]; then mkdir tmp/pmats; fi")
        shell("if [ ! -d tmp/objects ]; then mkdir tmp/objects; fi")
        shell("if [ ! -d tmp/reports ]; then mkdir tmp/reports; fi")
        shell("Rscript {params.scripts_dir}/generate.QC.report.3.R {params.outdir} {params.prefix} {params.RNAdir}")
        shell("Rscript -e \"rmarkdown::render('{params.rmd}', params = list(directory = '{params.outdir}', file='{params.prefix}'), output_file = '{params.rmdout}')\" 2> {log}")
