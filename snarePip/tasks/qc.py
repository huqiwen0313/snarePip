import os
import gspread
import pandas as pd
import numpy as np
import re
from snarePip.tasks.metable import *
from snarePip.tasks.snakemakeRun import *
from snarePip.io import prepare_submission


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

        return LocalTarget(os.path.join(self.tmp_dir, "qclogs.csv"))

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
        # need furthur modification to replace specific rows
        logs.append(sample_table.update([runlist.columns.values.tolist()] + runlist.values.tolist()))
        # update uploading status
        runlist = pd.read_csv(os.path.join(self.tmp_dir, "runlist.csv"), sep=",")
        runlist['qcid'] = 1
        runlist.to_csv(self.output().path, index=False, header=True)


@inherits(UpdateQC)
class GetUploadSamples(Task):
    """Get sample information and dirs that need to be uploaded """
    # names of submission table contain sample uploading dirs
    subTableName = Parameter(default="hubmap_submission")
    tmp_dir = "./temp"

    def requires(self):
        return{
            'submission': ContentMetaTable(sheetname=self.subTableName, worksheet=0),
            'qc_pipeline': self.clone(UpdateQC)
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(self.tmp_dir, "upload.logs.csv"))

    def run(self):
        temp_sample_file = os.path.join(self.tmp_dir, "runlist.csv")
        # get processed_ids
        processed_ids_df = pd.read_csv(temp_sample_file, sep=",")

        # get submission table
        submission_sheet = self.input()['submission']
        submission_df = pd.DataFrame(submission_sheet.get_all_records(head=True))
        submission_df['Globus location'] = [re.sub(r' ', r"\ ", i)
                                            for i in submission_df['Globus location']]
        submission_df = submission_df[['Experiment_ID_Short', 'Globus location']]

        # merge tables and get samples with submission dir registered
        submission_df_selected = pd.concat([processed_ids_df, submission_df], axis=1, join='inner')
        submission_df_selected = submission_df_selected.loc[submission_df_selected['Globus location'] != '']
        submission_df_selected.to_csv(self.output().path, index=False, header=True)


@inherits(GetUploadSamples)
class generateUploadFiles(Task):
    """prepare uploading folders"""
    # name of contributor table
    cTableName = Parameter(default="contributor")
    # path to submission_metatable
    #submission_meta = Parameter(default="./hubmap_meta.csv")
    tmp_dir = "./temp"

    def requires(self):
        return{
            'contributors': ContentMetaTable(sheetname=self.cTableName, worksheet=0),
            'upload_log': self.clone(GetUploadSamples)
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(self.tmp_dir, "upload.log.txt"))

    def run(self):
        # read contributor table
        csheet = self.input()['contributors']
        contributor_df = pd.DataFrame(csheet.get_all_records(head=True))

        # get upload list
        upload_df = pd.read_csv(self.input()['upload_log'].path, sep=",")
        upload_rna = upload_df.loc[upload_df['Type'] == "RNA"]
        upload_atac = upload_df.loc[upload_df['Type'] == "ATAC"]

        if upload_rna.shape[0] > 0:
            rna_path = self.RNAdir
            atac_path = self.ATACdir
            prepare_submission(rna_path, self.input()['upload_log'].path,
                               contributor_df, type="RNA")

        if upload_atac.shape[0] > 0:
            prepare_submission(atac_path, self.input()['upload_log'].path,
                               contributor_df, type="ATAC")

        upload_df["success"] = 1
        upload_df.to_csv(self.output().path, index=False, header=True)
