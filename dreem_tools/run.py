import click
import pandas as pd
import glob
import subprocess
import shutil
import os
import re
import json
import pickle
from tabulate import tabulate

from dreem import bit_vector
from dreem_tools import logger

log = logger.setup_applevel_logger()


def get_num(s):
    return float(re.findall(r"[-+]?\d*\.\d+|\d+", s)[0])


def unique_sec(all_dirs):
    org = all_dirs.pop(0)
    spl1 = org.split("_")
    diff_pos = []
    diff = []
    for d in all_dirs:
        spl2 = d.split("_")
        for i in range(0, len(spl1)):
            if spl1[i] == spl2[i]:
                continue
            diff_pos.append(i)
            diff.append(get_num(spl2[i]))
    if len(diff_pos) != len(all_dirs):
        print("cannot find unique section")
        exit()
    diff.insert(0, get_num(org.split("_")[diff_pos[0]]))
    all_dirs.insert(0, org)
    return diff


def process_seq_csv(path):
    df = pd.read_csv(path)
    seq = "".join(list(df["nuc"]))
    ss = "".join(list(df["struc"]))
    data = ";".join([str(x) for x in list(df["mismatches"])])
    return (seq, ss, data)


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
def runmulti(csv, data_dir, seq_path):
    if not os.path.isdir("processed"):
        os.mkdir("processed")
    if not os.path.isdir("analysis"):
        os.mkdir("analysis")
    df = pd.read_csv(csv)
    os.chdir("processed")
    dfs = []
    df_data = pd.DataFrame(columns="rna_name name code seq ss data".split())
    pos = 0
    for i, row in df.iterrows():
        name, code, barcode_seq, type = (
            row["name"],
            row["code"],
            row["barcode_seq"],
            row["type"],
        )
        dir_name = row["name"] + "_" + row["code"] + "_" + row["type"]
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        os.chdir(dir_name)
        fq1 = f"{data_dir}/{barcode_seq}/test_S1_L001_R2_001.fastq"
        fq2 = f"{data_dir}/{barcode_seq}/test_S1_L001_R1_001.fastq"
        fa = f"{seq_path}/fastas/{code}.fasta"
        db = f"{seq_path}/rna/{code}.csv"
        cmd = f"dreem -fa {fa} -fq1 {fq1} -fq2 {fq2} --dot_bracket {db} --plot_sequence "
        if "align_type" in row:
            if row["align_type"] == "NORM":
                cmd += "--map_score_cutoff 15 "
            elif row["align_type"] == "SHAPEMAPPER":
                cmd += (
                    '--map_score_cutoff 25 --bt2_alignment_args "--local;--no-unal;--no-mixed;--no-discordant;-i S,1,0.75;'
                    '--rdg 5,1;--rfg 5,1;--ignore-quals;--dpad 30;--maxins 800;-R 10;-N 0;-D 15;-p 16"'
                )
        else:
            cmd += "--map_score_cutoff 15"
        print(cmd)
        try:
            subprocess.call(cmd, shell=True)
        except:
            log.error("CMD did not run")
            continue
        df_sum = pd.DataFrame()
        try:
            df_sum = pd.read_csv("output/BitVector_Files/summary.csv")
        except:
            os.chdir("..")
        df_sum["dir"] = dir_name
        df_sum["code"] = code
        df_sum["rna_name"] = name
        df_sum["barcode"] = row["barcode"]
        df_sum["barcode_seq"] = barcode_seq
        df_sum["type"] = type
        datas = []
        seqs = []
        sss = []
        for j, row2 in df_sum.iterrows():
            path = f"output/BitVector_Files/{row2['name']}_*.csv"
            seq, ss, data = process_seq_csv(glob.glob(path)[0])
            datas.append(data)
            seqs.append(seq)
            sss.append(ss)
            pos += 1
        df_sum["sequence"] = seqs
        df_sum["structure"] = sss
        df_sum["data"] = datas
        dfs.append(df_sum)
        os.chdir("..")
    df_sum_final = pd.concat(dfs)
    df_sum_final.to_csv("../analysis/summary.csv", index=False)
    df_sum_final["data"] = [
        [round(float(x), 5) for x in d.split(";")] for d in df_sum_final["data"]
    ]
    df_sum_final.to_json(f"../analysis/summary.json", orient="records")


@cli.command()
@click.argument("data_path")
def group(data_path):
    # print("SUMMARY:\n" + tabulate(paired, ['condition', 'dir'], tablefmt="github"))
    df = pd.read_csv(data_path)
    for i, group in df.groupby(["code", "name"]):
        cols = ["condition"]
        secs = unique_sec(list(df["rna_name"]))
        group["condition"] = secs
        group = group.sort_values("condition")
        for j, row in group.iterrows():
            if len(cols) == 1:
                nucs = list(row["seq"])

            data = row["data"].split(";")


def load(pathname):
    """Safe loader for binary pickle files"""
    fh = open(pathname, "rb")
    data = pickle.load(fh)
    fh.close()
    return data


@cli.command()
@click.argument("pickle_files", nargs=-1)
@click.option("-o", "--output", default="output.p")
def merge(pickle_files, output):
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
