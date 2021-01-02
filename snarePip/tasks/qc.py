import os
import gspread
import pandas as pd
import numpy as np
import re
import time
from snarePip.tasks.metable import *
from snarePip.tasks.snakemakeRun import *
from snarePip.io import prepare_submission
import subprocess
from datetime import datetime

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
        runlist = pd.read_csv(os.path.join(self.tmp_dir, "runlist.csv"), na_filter=False)
        # get finished sample ID
        samples = runlist['Experiment_ID_Short'].unique()

        # get experiment and sample leve QCs
        logs = []
        count = 0;
        c = 0
        for sample in samples:
            tissue = runlist.loc[runlist['Experiment_ID_Short'] == sample]['Tissue'].unique()
            rna_path = os.path.join(self.RNAdir, self.assayType, tissue[0], "samples", sample)
            atac_path = os.path.join(self.ATACdir, self.assayType, tissue[0], "samples", sample)

            # get experiment/sample level RNA QCs
            rna_qc_experiment = get_qc_stat(path=rna_path, sample_id=sample, level="experiment")
            rna_qc_sample = get_qc_stat(path=rna_path, sample_id=sample, level="sample")

            # get experiment/sample level ATAC QCs
            atac_qc_experiment = get_qc_stat(path=atac_path, sample_id=sample, level="experiment", type="ATAC")
            atac_qc_sample = get_qc_stat(path=atac_path, sample_id=sample, level="sample", type="ATAC")

            # updating qc tables
            # for RNA
            rna_experiment_sheet = self.input()['rna_experiement']
            rna_sample_sheet = self.input()['rna_sample']
            # update table
            rna_sample_sheet.append_row(rna_qc_sample)
            logs.append(list(map(rna_experiment_sheet.append_row, rna_qc_experiment)))
            time.sleep(100)

            # for ATAC
            atac_experiment_sheet = self.input()['atac_experiement']
            atac_sample_sheet = self.input()['atac_sample']
            logs.append(atac_sample_sheet.append_row(atac_qc_sample))
            logs.append(list(map(atac_experiment_sheet.append_row, atac_qc_experiment)))
            time.sleep(100)


        # update sample link table
        sample_table = self.input()['sample_table']
        # update sample tables in corespondent columns
        sample_table_df = pd.read_csv(os.path.join(self.tmp_dir, "process_samplelist.csv"),
                                      sep="\t", na_filter=False)
        # set values for processed samples
        findex = sample_table_df[sample_table_df['Experiment_ID'].isin(runlist['Experiment_ID'])].index
        sample_table_df.loc[findex, 'flag'] = runlist['flag'].tolist()
        sample_table_df.loc[findex, 'runid'] = runlist['runid'].tolist()

        # need furthur modification to replace specific rows
        logs.append(sample_table.update([sample_table_df.columns.values.tolist()] + sample_table_df.values.tolist()))
        
        # update uploading status
        runlist = pd.read_csv(os.path.join(self.tmp_dir, "runlist.csv"), sep=",")
        runlist['qcid'] = 1
        runlist.to_csv(self.output().path, index=False, header=True)


@inherits(UpdateQC)
class generateUploadFiles(Task):
    """prepare uploading folders"""
    # name of contributor table
    cTableName = Parameter(default="contributor")
    tmp_dir = "./temp"

    def requires(self):
        return{
            'contributors': ContentMetaTable(sheetname=self.cTableName, worksheet=0),
            'qc_log': self.clone(UpdateQC)
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(self.tmp_dir, "processed.folder.log.csv"))

    def run(self):
        # read contributor table
        csheet = self.input()['contributors']
        contributor_df = pd.DataFrame(csheet.get_all_records(head=True))

        # get upload list
        temp_sample_file = os.path.join(self.tmp_dir, "runlist.csv")

        # get processed_ids
        processed_ids_df = pd.read_csv(temp_sample_file, sep=",")
        upload_rna = processed_ids_df.loc[processed_ids_df['Type'] == "RNA"]
        upload_atac = processed_ids_df.loc[processed_ids_df['Type'] == "ATAC"]

        if upload_rna.shape[0] > 0:
            rna_path = self.RNAdir
            atac_path = self.ATACdir
            prepare_submission(rna_path, temp_sample_file,
                               contributor_df, type="RNA")

        if upload_atac.shape[0] > 0:
            prepare_submission(atac_path, temp_sample_file,
                               contributor_df, type="ATAC")

        processed_ids_df["success"] = 1
        processed_ids_df.to_csv(self.output().path, index=False, header=True)


@inherits(generateUploadFiles)
class CleanFiles(Task):
    """Clean processed files """

    tmp_dir = "./temp"

    def requires(self):
        return{
            'processed_files': self.clone(generateUploadFiles)
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        date_time = datetime.now().strftime("%d_%m_l%Y_%H_%M")
        return LocalTarget(os.path.join(self.tmp_dir, date_time + "_processed.log.csv"))

    def run(self):
        # clean log files
        rna_path = self.RNAdir
        atac_path = self.ATACdir
        subprocess.call("mv %s %s" % (os.path.join(rna_path, "tmp", "logs") + "/*", os.path.join(rna_path, "logs")),
                        shell=True)
        subprocess.call("mv %s %s" % (os.path.join(atac_path, "tmp", "logs") + "/*", os.path.join(atac_path, "logs")),
                        shell=True)

        # clean raw fastq files
        subprocess.call("rm " + os.path.join(rna_path, "by_samples_fastq") + "/*",
                        shell=True)
        subprocess.call("rm " + os.path.join(atac_path, "by_samples_fastq") + "/*",
                        shell=True)

        # clean tmp folder
        subprocess.call("rm " + os.path.join(rna_path, "tmp") + "/*" + " -rf",
                        shell=True)
        subprocess.call("rm " + os.path.join(atac_path, "tmp") + "/*" + " -rf",
                        shell=True)

        # clean temp folder
        temp_sample_file = os.path.join(self.tmp_dir, "runlist.csv")
        process_df = pd.read_csv(temp_sample_file, sep=",")
        subprocess.call("rm " + self.tmp_dir + "/*",
                        shell=True)
        process_df.to_csv(self.output().path, index=False, header=True)

