import os
import gspread
import pandas as pd
import numpy as np
import re
from snarePip.tasks.metable import *
from snarePip.tasks.snakemakeRun import *


@inherits(snakemakeATAC)
class UpdateQC(Task):
    """Update QC statistics after a success pipeline run,
    the default structure of QC table contains 4 spreadsheet: ATAC-experimental_level QC,
    RNA-experimental_level QC, ATAC-sample_level QC, RNA-sample_level QC"""

    # record the first spreadsheet index in QC table
    sindex = IntParameter(default=0)
    qcTableName = Parameter(default="template.test.data")
    assayType = Parameter(default="snare_2")
    tmp_dir = "./temp"

    def requires(self):
        return {
            'atac_experiement': ContentMetaTable(sheetname=self.qcTableName, worksheet=self.sindex),
            'rna_experiement': ContentMetaTable(sheetname=self.qcTableName, worksheet=self.sindex+1),
            'atac_sample': ContentMetaTable(sheetname=self.qcTableName, worksheet=self.sindex + 2),
            'rna_sample': ContentMetaTable(sheetname=self.qcTableName, worksheet=self.sindex + 3),
            'sample_table': ContentMetaTable(sheetname=self.sheet_name, worksheet=self.worksheet),
            'pipeline': self.clone(snakemakeATAC)
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(self.tmp_dir, "run.logs.csv"))

    def run(self):
        runlist = pd.read_csv(os.path.join(self.tmp_dir, "runlist.csv"))
        # get finished sample ID
        samples = runlist['Experiment_ID_Short'].unique()

        # get experiment and sample leve QCs
        logs = []
        for sample in samples:
            tissue = runlist.loc[runlist['Experiment_ID_Short'] == sample]['Tissue'].unique()
            rna_path = os.path.join(self.RNAdir, self.assayType, tissue[0], "samples", sample)
            atac_path = os.path.join(self.ATACdir, self.assayType, tissue[0], "samples", sample)

            # get experiment/sample level RNA QCs
            rna_qc_experiment = get_qc_stat(path=rna_path, sample_id=sample, level="experiment")
            rna_qc_sample = get_qc_stat(path=rna_path, sample_id=sample, level="sample")

            # get experiment/sample level ATAC QCs
            atac_qc_experiment = get_qc_stat(path=atac_path, sample_id=sample, level="experiment")
            atac_qc_sample = get_qc_stat(path=atac_path, sample_id=sample, level="sample")

            # updating qc tables
            # for RNA
            rna_experiment_sheet = self.input()['rna_experiement']
            rna_sample_sheet = self.input()['rna_sample']
            # update table
            rna_sample_sheet.append_row(atac_qc_sample)
            logs.append(list(map(rna_experiment_sheet.append_row, rna_qc_experiment)))

            # for ATAC
            atac_experiment_sheet = self.input()['atac_experiement']
            atac_sample_sheet = self.input()['atac_sample']
            logs.append(atac_sample_sheet.append_row(atac_qc_sample))
            logs.append(list(map(atac_experiment_sheet.append_row, atac_qc_experiment)))

        # update sample link table
        sample_table = self.input()['sample_table']
        logs.append(sample_table.update([runlist.columns.values.tolist()] + runlist.values.tolist()))
