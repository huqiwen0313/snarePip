rule merge_to_sample_level_9:
    input:
        get_sample_wildcards
    output:
        join(ASSAY, "{tissue}/samples/{sample_name}/Sample_output/objects/{sample_name}.p2.obj.rds")
    params:
        ATACdir = CURRENT_DIR,
        delete_tmp = "FALSE",
        type = "ATAC",
        link = "TRUE",
        script_dir = config['scripts_dir']
    run:
        shell("Rscript {params.scripts_path}/distr.sample.R {params.ATACdir} {params.delete_tmp} {params.type} {params.link}")
