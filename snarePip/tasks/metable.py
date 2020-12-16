import os
import gspread
import pandas as pd
import numpy as np
import re
from luigi.util import inherits
from luigi import ExternalTask, LocalTarget, Parameter, Task, IntParameter
from snarePip.io import jsonkey_key_path, googlesheet_connection, get_qc_stat
from snarePip.luigi.task import Requirement, Requires, TargetOutput, ListTarget


def getflag():
    return 1


class ContentMetaTable(ExternalTask):
    # get spreadsheet through googlesheet API
    sheetname = Parameter(default="template.test.data")
    worksheet = IntParameter(default=0)

    def output(self):
        # return the connected spreadsheet from google API
        # get credential json key
        jsonkey_file = jsonkey_key_path()
        sheet = googlesheet_connection(spreadsheetname=self.sheetname,
                                       jsonpath=jsonkey_file, worksheet=self.worksheet)
        # flag sheet object
        setattr(sheet, "exists", getflag)
        return sheet


class TargetFolder(ExternalTask):
    """returns lists of the raw fastqs in target folder"""
    folder_dir = Parameter()

    def output(self):
        # return list of raw fastq files
        fq_files = os.listdir(self.folder_dir)
        # get unique samples
        samples = np.unique(list(map(lambda st: re.sub("_S.*", "", st), fq_files)))
        # replace "P\d" for some samples
        return ListTarget((map(lambda st: re.sub("\.P[\d+]", "", st), samples)))


class GetNewSample(Task):
    # get new sample information to start data processing
    sheet_name = Parameter(default="template.test.data")
    worksheet = IntParameter(default=4)

    def requires(self):
        return ContentMetaTable(sheetname=self.sheet_name, worksheet=self.worksheet)

    def output(self):
        tmp_dir = "./temp"
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)

        return LocalTarget(os.path.join(tmp_dir, "process_samplelist.csv"))

    def run(self):
        sheet = self.input()
        samples_df = pd.DataFrame(sheet.get_all_records(head=True))

        # get new samples need to be process
        process_samples_df = samples_df.loc[samples_df['flag'] == '']
        process_samples_df['Experiment_ID'] = list(map(lambda st: re.sub("SNARE2.*_", "", st),
                                                       process_samples_df['Experiment_ID']))
        process_samples_df.to_csv(self.output().path, index=False, header=True)


class CheckSample(Task):
    """check if samples that need to process exist in the target folder"""
    fastq_folder = "by_samples_fastq"
    RNAdir = Parameter(default=os.path.join("./test_RNA_dir", fastq_folder))
    ATACdir = Parameter(default=os.path.join("./test_ATAC_dir", fastq_folder))
    sheet_name = Parameter(default="template.test.data")
    worksheet = IntParameter(default=0)

    def requires(self):
        return {
                'newSample': GetNewSample(sheet_name=self.sheet_name, worksheet=self.worksheet),
                'rnaSamples_folder': TargetFolder(folder_dir=os.path.join(self.RNAdir, self.fastq_folder)),
                'atacSamples_folder': TargetFolder(folder_dir=os.path.join(self.ATACdir, self.fastq_folder))
        }

    def output(self):
        tmp_dir = "./temp"
        if not os.path.exists(tmp_dir):
            raise Exception("temp file does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(tmp_dir, "runlist.csv"))

    def run(self):
        process_samples_df = pd.read_csv(self.input()['newSample'].path, sep=",")
        rna_sample_list = self.input()['rnaSamples_folder'].getobj()
        atac_sample_list = self.input()['atacSamples_folder'].getobj()

        # get RNA and ATAC samples that need to be processed
        rna_processed_list = process_samples_df.loc[process_samples_df['Type'] == "RNA"]["Experiment_ID"]
        rna_processed_list = list(set(rna_sample_list) & set(rna_processed_list))
        atac_processed_list = process_samples_df.loc[process_samples_df['Type'] == "ATAC"]["Experiment_ID"]
        atac_processed_list = list(set(atac_sample_list) & set(atac_processed_list))

        # output samples that need to be processed in the raw fastq folder
        samples_df = process_samples_df.loc[process_samples_df['Experiment_ID'].isin(rna_processed_list + atac_processed_list)]
        samples_df['runid'] = None
        samples_df.to_csv(self.output().path, index=False, header=True)

