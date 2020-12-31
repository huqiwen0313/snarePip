import os
import subprocess
from decouple import config
from luigi import ExternalTask, LocalTarget, Parameter, Task, IntParameter
from snarePip.io import *
from snarePip.tasks.metable import *


class GetUploadSamples(Task):
    """Get sample information and dirs that need to be uploaded """
    # names of submission table contain sample uploading dirs
    subTableName = Parameter(default="hubmap_submission")
    # names of submission meta-tables
    subTableMeta = Parameter(default="hubmap_submission_metatable")

    RNAdir = Parameter(default="./test_RNA_dir")
    ATACdir = Parameter(default="./test_atac_dir")

    tmp_dir = "./temp"
    submission_folder = "data_submission"

    def requires(self):
        return{
            'submission': ContentMetaTable(sheetname=self.subTableName, worksheet=0),
            'rna_submission_meta': ContentMetaTable(sheetname=self.subTableMeta, worksheet=0),
            'atac_submission_meta': ContentMetaTable(sheetname=self.subTableMeta, worksheet=1),
            'rna_local_folder': TargetFolder(folder_dir=os.path.join(self.RNAdir, self.submission_folder)),
            'atac_local_folder': TargetFolder(folder_dir=os.path.join(self.ATACdir, self.submission_folder))
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(self.tmp_dir, "upload.logs.csv"))

    def run(self):
        rna_upload_list = list(self.input()['rna_local_folder'].getobj())
        atac_upload_list = list(self.input()['atac_local_folder'].getobj())

        # get submission table
        submission_sheet = self.input()['submission']
        submission_df = pd.DataFrame(submission_sheet.get_all_records(head=True))
        # output submission sheet
        submission_df.to_csv(os.path.join(self.tmp_dir, "submission.csv"), index=False, header=True)

        submission_df['Globus location'] = [re.sub(r' ', r"\ ", i)
                                            for i in submission_df['Globus location']]
        submission_df = submission_df[['Experiment_ID_Short', 'Globus location', 'Type', 'flag']]

        # get samples with submission dir registered
        submission_df_selected = submission_df.loc[submission_df['Experiment_ID_Short'].isin(self.input()['rna_local_folder'].getobj())]
        submission_df_selected = submission_df_selected.loc[(submission_df_selected['Globus location'] != '') &
                                                            (submission_df_selected['flag'] != 1)]

        # distribute hubmap submission meta-table
        atac_submission_path = os.path.join(self.ATACdir, self.submission_folder)
        rna_submission_path = os.path.join(self.RNAdir, self.submission_folder)

        rna_submission_meta_sheet = self.input()['rna_submission_meta']
        rna_submission_meta_df = pd.DataFrame(rna_submission_meta_sheet.get_all_records(head=True))
        atac_submission_meta_sheet = self.input()['atac_submission_meta']
        atac_submission_meta_df = pd.DataFrame(atac_submission_meta_sheet.get_all_records(head=True))

        rna_submission_meta_samplelist = list(map(lambda st: re.sub("_N.*", "", st),
                                                  rna_submission_meta_df['data_path']))
        atac_submission_meta_samplelist = list(map(lambda st: re.sub("_N.*", "", st),
                                                  atac_submission_meta_df['data_path']))
        for s in list(rna_upload_list):
            rna_indexes = np.where(np.array(rna_submission_meta_samplelist) == s)
            rna_table = rna_submission_meta_df.iloc[rna_indexes[0].tolist(), ]
            rna_table.to_csv(os.path.join(rna_submission_path, s, s+".scrnaseq-metadata.tsv"),
                             index=False, header=True)
            atac_indexes = np.where(np.array(atac_submission_meta_samplelist) == s)
            atac_table = atac_submission_meta_df.iloc[atac_indexes[0].tolist(), ]
            atac_table.to_csv(os.path.join(atac_submission_path, s, s + ".scatacseq-metadata.tsv"),
                              index=False, header=True)

        submission_df_selected.to_csv(self.output().path, index=False, header=True)


@inherits(GetUploadSamples)
class UploadSamples(Task):
    """upload samples through globus """

    tmp_dir = "./temp"
    submission_folder = "data_submission"

    # get endpoint
    if os.path.exists(".env"):
        dpoint = config("dpoint")
        upoint = config("upoint")
    else:
        raise Exception(".env file does not exist, please specify dpoint and upoint in .env file")

    def requires(self):
        return {
            'upload_samples': self.clone(GetUploadSamples),
            'submission': ContentMetaTable(sheetname=self.subTableName, worksheet=0)
        }

    def output(self):
        if not os.path.exists(self.tmp_dir):
            raise Exception("temp folder does not exist, please check ancestor tasks")

        return LocalTarget(os.path.join(self.tmp_dir, "uploaded.samples.csv"))

    def run(self):
        upload_list_df = pd.read_csv(self.input()['upload_samples'].path, sep=",",
                                     na_filter=False)
        rna_upload_list_df = upload_list_df.loc[upload_list_df['Type'] == "RNA"]
        atac_upload_list_df = upload_list_df.loc[upload_list_df['Type'] == "ATAC"]

        for i in rna_upload_list_df['Experiment_ID_Short']:
            globdir = rna_upload_list_df.loc[rna_upload_list_df['Experiment_ID_Short'] == i]
            denst = self.dpoint + ":" + globdir['Globus location']

            upload_folder_path = os.path.join(self.RNAdir, self.submission_folder, i)
            upload_dir = self.upoint + ":" + upload_folder_path
            subprocess.call("globus transfer --recursive " + upload_dir + denst, shell=True)

        for i in atac_upload_list_df['Experiment_ID_Short']:
            globdir = atac_upload_list_df.loc[atac_upload_list_df['Experiment_ID_Short'] == i]
            denst = self.dpoint + ":" + globdir['Globus location']

            upload_folder_path = os.path.join(self.ATACdir, self.submission_folder, i)
            upload_dir = self.upoint + ":" + upload_folder_path
            print("globus transfer --recursive " + upload_dir + " " + denst)
            subprocess.call("globus transfer --recursive " + upload_dir + " " + denst, shell=True)

        # update submission table
        submission_df = pd.read_csv(os.path.join(self.tmp_dir, "submission.csv"),
                                    na_filter=False)
        findex = submission_df[submission_df['Experiment_ID_Short'].isin(rna_upload_list_df['Experiment_ID_Short'])].index
        submission_df.loc[findex, 'flag'] = 1

        # update status
        logs = []
        submission = self.input()['submission']
        logs.append(submission.update([submission_df.columns.values.tolist()] + submission_df.values.tolist()))

        upload_list_df.to_csv(self.output().path, index=False, header=True)
