import pandas as pd
from collections import namedtuple


def get_data_from_pop_avg(df):
    seq = "".join(list(df["nuc"]))
    ss = "".join(list(df["struc"]))
    MutationData = namedtuple(
        "MutationData", ["sequence", "structure", "mismatches", "mismatch_del"]
    )
    return MutationData(seq, ss, list(df["mismatches"]), list(df["mismatch_del"]))


def get_avg_reactivity(df, p5_cut=20, p3_cut=20):
    avg = 0
    count = 0
    df_sub = df[p5_cut:-p3_cut]
    for i, row in df_sub.iterrows():
        if row["nuc"] == "A" or row["nuc"] == "C":
            avg += row["mismatches"]
            count += 1
    return avg / count


def get_avg_noise(df, p5_cut=20, p3_cut=20):
    avg = 0
    count = 0
    df_sub = df[p5_cut:-p3_cut]
    for i, row in df_sub.iterrows():
        if row["nuc"] == "T" or row["nuc"] == "G":
            avg += row["mismatches"]
            count += 1
    return avg / count
