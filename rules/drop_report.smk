rule F_dropReport_7:
    input:
        pagoda_RData = 'tmp/pagoda_RData/{sample}.filtered.RData',
        seurat_obj = 'tmp/seurat_obj/{sample}.seurat.filtered.RData'
    output:
        html = 'tmp/reports/{sample}_dropEst_QC.html'
    params:
        output_prefix=join(CURRENT_DIR, 'tmp/reports/{sample}_dropEst_QC'),
        scripts_path = config['scripts_dir'] 
    log:
        'tmp/log/F_dropReport_7_{sample}.Log.out'
    run:
        shell("{params.scripts_path}/dropReport_lib_one.Rsc -o '{params.output_prefix}' {input.seurat_obj} {input.pagoda_RData} > {log}")
