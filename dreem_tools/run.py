import click
import pandas as pd
import glob
import subprocess
import shutil
import os
import yaml

from dreem_tools import logger, util, parse
from dreem_tools.util import load


def parse_demultiplex_log(lines):
    data = []
    for l in lines:
        if l.startswith("#"):
            continue
        spl = l.split()
        # check if last row has number of reads
        try:
            num = float(spl[-1])
        except:
            continue
        data.append(spl)
    df = pd.DataFrame(data, columns="id,tag,count".split(","))
    df["count"] = pd.to_numeric(df["count"])
    df.to_csv("demultiplex.csv", index=False)
    return df


def generate_barcode_file(df):
    """
    generates a .fa file for novobarcode
    :param df: a dataframe with that contains informmation of barcode sequences
    :return: None
    """
    log = logger.get_logger("RUN")
    expects = ["barcode", "barcode_seq", "construct"]
    for e in expects:
        if e not in df:
            raise ValueError(f"{e} column is required in the csv")
    s = "Distance\t4\nFormat\t5\n"
    seen = []
    warning = False
    for i, row in df.iterrows():
        if row["barcode"] in seen:
            log.warning(
                f"{row['barcode']} has been used more than once this may be an issue"
            )
            warning = True
            continue
        s += f"{row['barcode']}\t{row['barcode_seq']}\n"
        seen.append(row["barcode"])
    log.info(f"{len(seen)} unique barcodes found from csv file")
    if not warning:
        log.info("no barcode conflicts detected")
    f = open("rtb_barcodes.fa", "w")
    f.write(s)
    f.close()


@click.group()
def cli():
    pass


@cli.command(help="download a run using bs commandline tool")
@click.argument("run_name")
@click.option(
    "-d",
    "--dir",
    default=None,
    help="root directory of where data should be downloaded will default to "
    "$BASESPACE",
)
def download(run_name, dir):
    log = logger.setup_applevel_logger()
    if dir is None:
        log.info("-d/--dir was not suppled will use $BASESPACE")
        if os.getenv("BASESPACE") is None:
            log.error("$BASESPACE is not set! set it or use -d/--dir")
        log.info("$BASESPACE -> " + os.getenv("BASESPACE"))
        dir = os.getenv("BASESPACE")
    os.chdir(dir)
    os.makedirs(run_name, exist_ok=True)
    os.chdir(run_name)
    log.info(f"running: `bs download project --name {run_name}`")
    os.system(f"bs download project --name {run_name}")


@cli.command()
@click.argument("csv")
@click.option("--debug", is_flag=True)
def demultiplex(csv, debug):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    log = logger.setup_applevel_logger(file_name="demultplex.log")
    log.info("preparing rtb_barcodes.fa file for demultiplexing")
    df = pd.read_csv(csv)
    generate_barcode_file(df)
    r1_fastq = glob.glob("*R1*")[0]
    r2_fastq = glob.glob("*R2*")[0]
    log.info(f"renaming {r1_fastq} -> test_S1_L001_R1_001.fastq")
    log.info(f"renaming {r2_fastq} -> test_S1_L001_R2_001.fastq")
    shutil.move(r1_fastq, "test_S1_L001_R1_001.fastq")
    shutil.move(r2_fastq, "test_S1_L001_R2_001.fastq")
    output = subprocess.check_output(
        "novobarcode -b rtb_barcodes.fa -f test_S1_L001_R1_001.fastq test_S1_L001_R2_001.fastq",
        shell=True,
    )
    output = output.decode("UTF-8")
    log.info(f"output from novobarcode:\n{output}")
    lines = output.split("\n")
    df_demult = parse_demultiplex_log(lines)
    log.info(f"total number of reads: {df_demult['count'].sum()}")
    log.info(
        f"total number of data reads: "
        f"{df_demult['count'].sum() - df_demult.iloc[-1]['count']}"
    )
    if not os.path.isdir("NC"):
        log.error("nothing was generated STOPING NOW!")
        exit()
    if not debug:
        os.remove("test_S1_L001_R1_001.fastq")
        os.remove("test_S1_L001_R2_001.fastq")
        shutil.rmtree("NC")


@cli.command()
@click.argument("csv")
@click.argument("data_dir")
@click.argument("seq_path")
@click.option("--move-plots", is_flag=True)
def runmulti(csv, data_dir, seq_path, move_plots):
    log = logger.setup_applevel_logger(file_name="run_multi.log")
    log.info("creating processed/ all dreem runs will go here")
    log.info("creating analysis/ all finalized analysis will go here")
    os.makedirs("processed", exist_ok=True)
    os.makedirs("analysis", exist_ok=True)
    df = pd.read_csv(csv)
    # catch old column naming scheme
    if "name" in df:
        df = df.rename({"name": "construct"}, axis="columns")
    if "type" in df:
        df = df.rename({"type": "data_type"}, axis="columns")
    os.chdir("processed")
    if move_plots:
        os.makedirs("plots", exist_ok=True)
    dfs = []
    for i, row in df.iterrows():
        dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
        os.makedirs(dir_name, exist_ok=True)
        os.chdir(dir_name)
        util.run_dreem_from_row(data_dir, seq_path, row)
        mhs = None
        try:
            mhs = load("output/BitVector_Files/mutation_histos.p")
        except:
            os.chdir("..")
            continue
        df_sum = util.mutation_histos_to_dataframe(mhs)
        all_cols = list(df.columns)
        all_cols.remove("construct")
        df_sum["dir"] = dir_name
        df_sum["rna_name"] = row["construct"]
        for col in all_cols:
            df_sum[col] = row[col]
        dfs.append(df_sum)
        os.chdir("..")
    df_sum_final = pd.concat(dfs)
    df_sum_final.to_json(f"../analysis/summary.json", orient="records")


@cli.command()
@click.argument("json_file")
@click.argument("yml_file")
def parsedata(json_file, yml_file):
    log = logger.setup_applevel_logger()
    df = pd.read_json(json_file)
    f = open(yml_file)
    params = yaml.load(f, Loader=yaml.FullLoader)
    parse.parse(df, params)
    log.info("written results in processed.json")


if __name__ == "__main__":
    cli()
