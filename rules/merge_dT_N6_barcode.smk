rule E_merge_dT_N6_barcode_6:
    input:
        rds = 'tmp/dropest_out/{sample}.rds'
    output:
        seurat_obj = 'tmp/seurat_obj/{sample}.seurat.filtered.RData',
        pagoda_RData = 'tmp/pagoda_RData/{sample}.filtered.RData'
    params:
        output_prefix_seurat = join(CURRENT_DIR, 'tmp/seurat_obj/{sample}'),
        output_prefix_pagoda = join(CURRENT_DIR, 'tmp/pagoda_RData/{sample}'),
        scripts_path = config['scripts_dir'],
        link_table = config['link_table'],
        dt_n6 = config['dt_n6']
    log:
        'tmp/logs/E_merge_dT_N6_barcode_5_{sample}.Log.out'
    run:
        shell("if [ ! -d tmp/pagoda_RData ]; then mkdir tmp/pagoda_RData; fi")
        shell("if [ ! -d tmp/seurat_obj ]; then mkdir tmp/seurat_obj; fi")
        shell("Rscript {params.scripts_path}/addLibID_generate_seurat_dropest_d.r '{input.rds}' '{params.output_prefix_seurat}' '{params.output_prefix_pagoda}' '{params.link_table}' '{params.dt_n6}'  > {log}")
