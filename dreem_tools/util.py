import numpy as np
import pandas as pd
import subprocess
import pickle

import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from dreem_tools import logger
from dreem import run as run_dreem

log = logger.get_logger("UTIL")


def load(pathname):
    """Safe loader for binary pickle files"""
    fh = open(pathname, "rb")
    data = pickle.load(fh)
    fh.close()
    return data


def run_dreem_from_row(data_dir, seq_path, row):
    args = run_dreem.get_default_run_args()
    args["fasta"] = f"{seq_path}/fastas/{row['code']}.fasta"
    args[
        "fastq1"
    ] = f"{data_dir}/{row['barcode_seq']}/test_S1_L001_R2_001.fastq"
    args[
        "fastq2"
    ] = f"{data_dir}/{row['barcode_seq']}/test_S1_L001_R1_001.fastq"
    args["dot_bracket"] = f"{seq_path}/rna/{row['code']}.csv"
    args["plot_sequence"] = True
    # only if included forcing this was annoying!
    if "align_type" in row:
        if row["align_type"] == "SHAPEMAPPER":
            args["map_score_cutoff"] = 25
            args["bt2_alignment_args"] = (
                "--local;--no-unal;--no-mixed;--no-discordant;-i S,1,0.75;"
                "--rdg 5,1;--rfg 5,1;--ignore-quals;--dpad 30;--maxins 800;-R "
                "10;-N 0;-D 15;-p 16"
            )
    else:
        args["map_score_cutoff"] = 15
    run_dreem.run(args)

def pop_avg_from_mutation_histo(mh):
    mut_y = []
    for pos in range(mh.start, mh.end + 1):
        try:
            mut_frac = mh.mut_bases[pos] / mh.info_bases[pos]
        except:
            mut_frac = 0.0
        mut_y.append(round(mut_frac, 5))
    return mut_y


def mutation_histos_to_dataframe(mhs):
    cols = "name,sequence,structure,data_type,num_reads,num_aligns".split(",")
    cols += "data,no_mut,1_mut,2_mut,3_mut,3plus_mut,sn".split(",")
    all_data = []
    for name, mh in mhs.items():
        data = []
        struct = np.nan
        if mh.structure is not None:
            struct = mh.structure
        sequence = mh.sequence.replace("T", "U")
        data.extend([name, sequence, struct, mh.data_type, mh.num_reads])
        data.extend([mh.num_aligned, pop_avg_from_mutation_histo(mh)])
        data.extend(mh.get_percent_mutations())
        data.append(mh.get_signal_to_noise())
        all_data.append(data)
    df = pd.DataFrame(all_data, columns=cols)
    return df


def plot_population_avg(mh, file_base_name):
    xaxis_coordinates = [i for i in range(mh.start, mh.end)]
    popavg_filename = file_base_name + "popavg_reacts.csv"
    popavg_file = open(popavg_filename, "w")
    popavg_file.write("position,mismatches,mismatch_del,nuc")
    if mh.structure is not None:
        popavg_file.write(",struc")
    popavg_file.write("\n")
    delmut_y, mut_y = [], []
    for pos in range(mh.start, mh.end + 1):
        try:
            delmut_frac = (
                mh.del_bases[pos] + mh.mut_bases[pos]
            ) / mh.info_bases[pos]
            mut_frac = mh.mut_bases[pos] / mh.info_bases[pos]
        except:
            delmut_frac = 0.0
            mut_frac = 0.0
        delmut_y.append(delmut_frac)
        mut_y.append(mut_frac)
        mut_frac, delmut_frac = round(mut_frac, 5), round(delmut_frac, 5)
        s = "{},{},{},{}".format(
            pos, mut_frac, delmut_frac, mh.sequence[pos - 1]
        )
        popavg_file.write(s)
        if mh.structure is not None:
            popavg_file.write("," + mh.structure[pos - 1])
        popavg_file.write("\n")
    popavg_file.close()

    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    colors = []
    ref_bases = []
    for i in range(mh.start, mh.end):
        if i >= len(mh.sequence):
            continue
        colors.append(cmap[mh.sequence[i - 1]])
        ref_bases.append(mh.sequence[i - 1])
    mut_trace = go.Bar(
        x=xaxis_coordinates,
        y=mut_y,
        text=ref_bases,
        marker=dict(color=colors),
        showlegend=False,
    )
    mut_fig_layout = go.Layout(
        title=mh.name,
        xaxis=dict(title="Postion"),
        yaxis=dict(title="Fraction"),
        plot_bgcolor="white",
        height=400,
    )
    mut_fig = go.Figure(data=mut_trace, layout=mut_fig_layout)
    seqs = list(mh.sequence)
    if mh.structure is not None:
        db = list(mh.structure)
    else:
        db = " " * len(seqs)
    mut_fig.update_yaxes(
        gridcolor="lightgray", linewidth=1, linecolor="black", mirror=True
    )
    mut_fig.update_xaxes(linewidth=1, linecolor="black", mirror=True)
    plotly.offline.plot(
        mut_fig,
        filename=file_base_name + "pop_avg.html",
        auto_open=False,
    )
    mut_fig.write_image(file_base_name + "pop_avg.png", height=400, width=1000)
