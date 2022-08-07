import numpy as np
import pandas as pd
import subprocess
import pickle

from dreem_tools import logger

log = logger.get_logger("UTIL")


def load(pathname):
    """Safe loader for binary pickle files"""
    fh = open(pathname, "rb")
    data = pickle.load(fh)
    fh.close()
    return data


def run_dreem_from_row(data_dir, seq_path, row):
    fq1 = f"{data_dir}/{row['barcode_seq']}/test_S1_L001_R2_001.fastq"
    fq2 = f"{data_dir}/{row['barcode_seq']}/test_S1_L001_R1_001.fastq"
    fa = f"{seq_path}/fastas/{row['code']}.fasta"
    db = f"{seq_path}/rna/{row['code']}.csv"
    cmd = (
        f"dreem -fa {fa} -fq1 {fq1} -fq2 {fq2} --dot_bracket {db}"
        f" --plot_sequence "
    )
    # only if included forcing this was annoying!
    if "align_type" in row:
        if row["align_type"] == "NORM":
            cmd += "--map_score_cutoff 15 "
        elif row["align_type"] == "SHAPEMAPPER":
            cmd += (
                '--map_score_cutoff 25 --bt2_alignment_args "--local;--no-unal;'
                "--no-mixed;--no-discordant;-i S,1,0.75;"
                "--rdg 5,1;--rfg 5,1;--ignore-quals;--dpad 30;--maxins 800;-R "
                '10;-N 0;-D 15;-p 16"'
            )
    else:
        cmd += "--map_score_cutoff 15"
    log.info("running dreem as: ")
    log.info(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except:
        log.error("command returned errors!")
        return None


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
