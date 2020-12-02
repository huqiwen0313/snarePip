rule merge_to_sample_level_9:
    input:
        #'tmp/QCs/{sample}.qc.txt'
        get_sample_wildcards
    output:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/obj/{sample_name}.sample_matrix.rds")
    params:
        RNAdir = CURRENT_DIR,
        delete_tmp = "FALSE",
        type = "RNA",
        link = "TRUE",
        scripts_path = config['scripts_dir']
    run:
        shell("Rscript {params.scripts_path}/distr.sample.R {params.RNAdir} {params.delete_tmp} {params.type} {params.link}")

rule dropReport_QC_sample_10:
    input:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/obj/{sample_name}.sample_matrix.rds")
    output:
        QCs = join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/QCs/{sample_name}.qc.txt")
    params:
        sampleID = "{sample_name}",
        path = join(CURRENT_DIR, ASSAY, "{tissue}/samples/{sample_name}/Sample_output"),
        scripts_path = config['scripts_dir']
    log:
        "tmp/logs/dropReport_QC_sample_10.{sample_name}.{tissue}.log"
    run:
        shell("Rscript {params.scripts_path}/RNA.sample.QC.R {params.path} {params.sampleID} > {log}")
