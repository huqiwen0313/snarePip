from snarePip.tasks.snakemakeRun import snakemakeRNA, snakemakeATAC
from snarePip.tasks.metable import TargetFolder
from snarePip.tasks.qc import *
from snarePip.tasks.upload import *
from luigi import build
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sampletable", default="template.test.data",
                    help="name of google sheet contains sample information")
parser.add_argument("-w", "--worksheet", default=0, type=int,
                    help="which worksheet to load from google spreadsheet")
parser.add_argument("-r", "--RNAdir", default="./by_samples_fastq",
                    help="path contains RNA raw fastq files")
parser.add_argument("-a", "--ATACdir", default="./by_samples_fastq",
                    help="path contains RNA raw fastq files")
parser.add_argument("-c", "--cores", default=1, type=int,
                    help="number of cores used to run snakemake pipeline")
parser.add_argument("-t", "--type", default="snare_2",
                    help="assay type, e.g. snare_2, tenX")
parser.add_argument("-sr", "--snakeRNA", default="Snakefile.RNA",
                    help="snakemake file for RNA processing")
parser.add_argument("-sa", "--snakeATAC", default="Snakefile.ATAC",
                    help="snakemake file for ATAC processing")
parser.add_argument("-sb", "--subtable", default="hubmap_submission",
                    help="name of data submission table")
parser.add_argument("-sc", "--ctable", default="contributor",
                    help="name of contributor table")
parser.add_argument("-sm", "--submeta", default="hubmap_submission_metatable",
                    help="name of hubmap meta table for submission")
parser.add_argument("-b", "--build", action='store_true')
parser.add_argument("-up", "--upload", action='store_true')
args = parser.parse_args()


def main(arg=None):
    if args.build:
        build([TargetFolder(
        folder_dir=args.RNAdir)], local_scheduler=True)
    elif args.upload:
        build([UploadSamples(
            RNAdir=args.RNAdir,
            ATACdir=args.ATACdir,
            subTableName=args.subtable,
            subTableMeta=args.submeta
            )], local_scheduler=True)
    else:
        build([generateUploadFiles(
            RNAdir=args.RNAdir,
            ATACdir=args.ATACdir,
            sheet_name=args.sampletable,
            worksheet=args.worksheet,
            ncores=args.cores,
            snakefileRNA=args.snakeRNA,
            snakefileATAC=args.snakeATAC,
            sindex=0,
            qcTableName=args.sampletable,
            assayType=args.type,
            cTableName=args.ctable
        )], local_scheduler=True)
        build([CleanFiles(
            RNAdir=args.RNAdir,
            ATACdir=args.ATACdir,
            sheet_name=args.sampletable,
            worksheet=args.worksheet,
            ncores=args.cores,
            snakefileRNA=args.snakeRNA,
            snakefileATAC=args.snakeATAC,
            sindex=0,
            qcTableName=args.sampletable,
            assayType=args.type,
            cTableName=args.ctable
        )], local_scheduler=True)

