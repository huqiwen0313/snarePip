rule H_R_stat_and_filtering_8:
    input:
        FQzip = 'tmp/FastQC/{sample}_R1_filtered_fastqc.zip',
        pagoda_RData = 'tmp/pagoda_RData/{sample}.filtered.RData'
    output:
        QCtxt2 = 'tmp/QCs/{sample}.qc.txt'
    params:
        sample= '{sample}',
        log_STAR2 = 'tmp/QCs/Align_stat/{sample}.Log.final.out',
        output_prefix_seurat = join(CURRENT_DIR, 'tmp/seurat_obj/{sample}'),
        scripts_path = config['scripts_dir']
    log:
        'tmp/log/H_R_stat_and_filtering_8_{sample}.Log.out'
    run:
        shell("if [ ! -d tmp/QCs ]; then mkdir tmp/QCs; fi")
        shell("if [ ! -d tmp/reports ]; then mkdir tmp/reports; fi")
        shell("Rscript {params.scripts_path}/QC_all_RNA_seq_by_sample_filtered_addLibID_LAPBUKMAP.r '{params.sample}' hg38 " "'{params.output_prefix_seurat}'")
