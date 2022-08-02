import click
import pandas as pd
import glob
import subprocess
import shutil
import os
import re
import pickle

from dreem import bit_vector
from dreem_tools import logger, util


def process_seq_csv(path):
    df = pd.read_csv(path)
    seq = "".join(list(df["nuc"]))
    ss = "".join(list(df["struc"]))
    data = ";".join([str(x) for x in list(df["mismatches"])])
    return (seq, ss, data)


def load(pathname):
    """Safe loader for binary pickle files"""
    fh = open(pathname, "rb")
    data = pickle.load(fh)
    fh.close()
    return data


@click.group()
def cli():
    pass


@cli.command()
@click.argument("csv")
def demultiplex(csv):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    s = "Distance\t4\nFormat\t5\n"
    df = pd.read_csv(csv)
    seen = []
    for i, row in df.iterrows():
        if row["barcode"] in seen:
            continue
        s += f"{row['barcode']}\t{row['barcode_seq']}\n"
        seen.append(row["barcode"])
    f = open("rtb_barcodes.fa", "w")
    f.write(s)
    f.close()
    r1_fastq = glob.glob("*R1*")[0]
    r2_fastq = glob.glob("*R2*")[0]
    shutil.move(r1_fastq, "test_S1_L001_R1_001.fastq")
    shutil.move(r2_fastq, "test_S1_L001_R2_001.fastq")
    subprocess.call(
            "novobarcode -b rtb_barcodes.fa -f test_S1_L001_R1_001.fastq test_S1_L001_R2_001.fastq",
            shell=True,
    )
    if not os.path.isdir("NC"):
        print("nothing was generated STOPING NOW!")
        exit()
    os.remove("test_S1_L001_R1_001.fastq")
    os.remove("test_S1_L001_R2_001.fastq")
    shutil.rmtree("NC")


@cli.command()
@click.argument("csv")
@click.argument("data_dir")
@click.argument("seq_path")
@click.option("--move-plots", is_flag=True)
def runmulti(csv, data_dir, seq_path, move_plots):
    # TODO write infomation about the run here
    log = logger.setup_applevel_logger()
    os.makedirs("processed", exist_ok=True)
    os.makedirs("analysis", exist_ok=True)
    df = pd.read_csv(csv)
    # catch old column naming scheme
    if "name" in df:
        df = df.rename({"name" : "construct"},  axis='columns')
    if "type" in df:
        df = df.rename({"type" : "data_type"},  axis='columns')
    os.chdir("processed")
    if move_plots:
        os.makedirs("plots", exist_ok=True)
    dfs = []
    for i, row in df.iterrows():
        print(row['construct'], row['code'], row['data_type'])
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
@click.argument("pickle_files", nargs=-1)
@click.option("-o", "--output", default="output.p")
def merge(pickle_files, output):
    log = logger.setup_applevel_logger()
    log.info(f"{len(pickle_files)} pickle files supplied")
    merged_mh = None
    for pf in pickle_files:
        if merged_mh is None:
            merged_mh = load(pf)
            continue
        mh = load(pf)
        log.info(pf)
        bit_vector.merge_mutational_histogram_dicts(merged_mh, mh)
    log.info(f"writing merged file: {output}")
    pickle.dump(merged_mh, open(output, "wb"))


if __name__ == "__main__":
    cli()
