import snakemake
import os
import io
import contextlib
import pandas as pd
from luigi import ExternalTask, LocalTarget, Parameter, Task, IntParameter
from luigi.util import inherits
from snarePip.luigi.task import Requirement, Requires, TargetOutput, ListTarget
from snarePip.tasks.metable import GetNewSample, CheckSample
from snarePip.io import generate_random_id


def runSnakemake(snakefile, cores=1, workdir=None, **kwargs):
    """function to run snakemake

    Args:
    snakefile (str): path to snakefile
    cores (int): the number of provided cores
    workdir (str): path to working directory
    """

    status = snakemake.snakemake(snakefile=snakefile, cores=cores, workdir=workdir,
                                 **kwargs)
    if status:
        print("Processing pipeline completed")
        return 1
    else:
        print("Pipeline stopped because an error occurs")
        return 0


@inherits(CheckSample)
class snakemakeRNA(Task):
    """run snakemake RNA"""

    snakefileRNA = Parameter(default="Snakefile.test")
    ncores = IntParameter(default=4)

    requires = Requires()
    check_sample = Requirement(CheckSample)

    # generate random ID associated with a run
    random_id = generate_random_id(length=6)

    def output(self):
        if not os.path.exists(os.path.join(self.RNAdir, "logs")):
            os.mkdir(os.path.join(self.RNAdir, "logs"))
        return LocalTarget(os.path.join(self.RNAdir, "logs",
                                        "snakemakeRNA_summary_" + self.random_id + ".csv"))

    def run(self):
        status = snakemake.snakemake(self.snakefileRNA, cores=self.ncores, workdir=self.RNAdir)
        if status:
            # get running status summary log file
            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                snakemake.snakemake(self.snakefileRNA, cores=self.ncores,
                                    workdir=self.RNAdir, detailed_summary=True)
            summary_report = f.getvalue().splitlines()
            summary_report = pd.DataFrame(list(map(lambda i: str.split(i, "\t"), summary_report)))
            summary_report.to_csv(self.output().path, index=False, header=False)

            # update sample list
            samples_run_df = pd.read_csv(self.input()['check_sample'].path)
            samples_run_df.loc[samples_run_df['Type'] == "RNA", 'flag'] = 1
            samples_run_df.loc[samples_run_df['Type'] == "RNA", 'runid'] = self.random_id
            samples_run_df.to_csv(self.input()['check_sample'].path, index=False, header=True)

        else:
            raise Exception("snakemake running errors, please check log files for details")


@inherits(snakemakeRNA)
class snakemakeATAC(Task):
    """run snakemake ATAC pipeline"""

    snakefileATAC = Parameter(default="Snakefile.test")

    requires = Requires()
    snake_RNA = Requirement(snakemakeRNA)
    check_sample = Requirement(CheckSample)

    # generate random ID associated with a run
    random_id = generate_random_id(length=6)

    def output(self):
        if not os.path.exists(os.path.join(self.ATACdir, "logs")):
            os.mkdir(os.path.join(self.ATACdir, "logs"))
        return LocalTarget(os.path.join(self.ATACdir, "logs",
                                        "snakemakeATAC_summary" + self.random_id + ".csv"))

    def run(self):
        status = snakemake.snakemake(self.snakefileATAC,
                                     cores=self.ncores, workdir=self.ATACdir)
        if status:
            # get running status summary log file
            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                snakemake.snakemake(self.snakefileATAC, cores=self.ncores,
                                    workdir=self.ATACdir, detailed_summary=True)
            summary_report = f.getvalue().splitlines()
            summary_report = pd.DataFrame(list(map(lambda i: str.split(i, "\t"), summary_report)))
            summary_report.to_csv(self.output().path, index=False, header=False)

            # update sample list
            samples_run_df = pd.read_csv(self.input()['check_sample'].path)
            samples_run_df.loc[samples_run_df['Type'] == "ATAC", 'flag'] = 1
            samples_run_df.loc[samples_run_df['Type'] == "ATAC", 'runid'] = self.random_id
            samples_run_df.to_csv(self.input()['check_sample'].path, index=False, header=True)

        else:
            raise Exception("snakemake running errors, please check log files for details")

