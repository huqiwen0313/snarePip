rule dual_omic_report_15:
    input:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/{sample_name}.qc.report.html")
    output:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/dual_omics/{sample_name}.dual.report.html")
    log:
        "tmp/logs/dual_omic_report_15.{sample_name}.{tissue}.log"
    params:
        atacPath = CURRENT_DIR + "/" + join(ASSAY, "{tissue}/samples"),
        rnaPath =  config['rna_path'] + "/" + join(ASSAY, "{tissue}/samples"),
        sampleID = '{sample_name}',
        rmd = config['scripts_dir'] + "/dual.omics.report.rmd",
        rmdout = CURRENT_DIR + join("/", ASSAY, "{tissue}/samples/{sample_name}/Sample_output/dual_omics/{sample_name}.dual.report.html")
    run:
        shell("Rscript -e \"rmarkdown::render('{params.rmd}', params=list(ATACpath='{params.atacPath}', RNApath='{params.rnaPath}', sampleID='{params.sampleID}'), output_file = '{params.rmdout}')\" 2> {log}")
