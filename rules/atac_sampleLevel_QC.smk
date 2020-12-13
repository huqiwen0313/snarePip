rule QC_feature_13:
    input:
        snap = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/snap/{sample_name}.snap")
    output:
        rmdout = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/QCs/{sample_name}.qc.txt")
    log:
        "tmp/logs/QC_feature_13.{sample_name}.{tissue}.log"
    params:
        atacDir = CURRENT_DIR,
        sampleID = '{sample_name}',
        RNAdir = config['rna_path'],
        scripts_path = config['scripts_dir'],
        reference_genome = config['genome_name'],
        assay_type = config['assay_type'],
        link_path = config['link_table'],
        ref_dir = config['ref_dir'],
        chromsize = config['chrom_size'],
        scrublet_path = config['scrublet_path']
    run:
        shell("Rscript {params.scripts_path}/sampleLevelQC.2.R {params.atacDir} {params.sampleID} {params.RNAdir} {params.reference_genome} {params.assay_type} {params.link_path} {params.ref_dir} {params.chromsize} {params.scrublet_path} > {log}")

rule sample_QC_report_14:
    input:
        qc = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/QCs/{sample_name}.qc.txt")
    output:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/{sample_name}.qc.report.html")
    log:
        "tmp/logs/sample_QC_report_14.{sample_name}.{tissue}.log"
    params:
        sampleID = '{sample_name}',
        sampleDir = CURRENT_DIR + "/" + join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output"),
        rmd = config['scripts_dir'] + "/ATAC.report.experiment.level.html.rmd",
        rmdout = CURRENT_DIR + join("/", ASSAY, "{tissue}/samples/{sample_name}/Sample_output/{sample_name}.qc.report.html")
    run:
        shell("Rscript -e \"rmarkdown::render('{params.rmd}', params = list(directory = '{params.sampleDir}', file='{params.sampleID}'), output_file = '{params.rmdout}')\" 2> {log}")
