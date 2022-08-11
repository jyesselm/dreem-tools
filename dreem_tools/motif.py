import pandas as pd
import rna_library as rl
import dacite
from dataclasses import dataclass

MotifType = rl.MotifType


@dataclass(frozen=True, order=True)
class MotifSearchParams(object):
    mtype: MotifType = None
    sequence: str = None
    structure: str = None
    min_pos: int = 0
    max_pos: int = 999
    min_id: int = 0
    max_id: int = 999
    extend: int = 0


def get_params_from_dict(d):
    if "type" in d:
        type_name = d.pop("type", None)
        if type_name == "ref_hp":
            d["sequence"] = "CGAGUAG"
        elif type_name == "gaaa_tetraloop":
            d["sequence"] = "GGAAAC"
        elif type_name == "tlr":
            d["sequence"] = "UAUG&CUAAG"
    if "mtype" in d:
        if type(d["mtype"]) == str:
            if d["mtype"].lower() == "twoway":
                d["mtype"] = MotifType.JUNCTION
            elif d["mtype"].lower() == "junction":
                d["mtype"] = MotifType.JUNCTION
            elif d["mtype"].lower() == "hairpin":
                d["mtype"] = MotifType.HAIRPIN
            elif d["mtype"].lower() == "singlestrand":
                d["mtype"] = MotifType.SINGLESTRAND
            else:
                raise ValueError(f"invalid mtype: {d['mtype']}")
    return dacite.from_dict(MotifSearchParams, d)


def _get_data_from_motif(m, row, extend):
    data = row["data"]
    strands = m.strands()
    if extend > 0:
        new_strands = []
        for s in strands:
            min_val, max_val = s[0], s[-1]
            if min_val - extend < 0:
                raise ValueError("cannot extend motif it will go beyond 0")
            if max_val + extend > len(data):
                raise ValueError("cannot extend motif it will go beyond length")
            r1 = list(range(min_val - extend, min_val))
            r2 = list(range(max_val + 1, max_val + extend + 1))
            new_strands.append(r1 + s + r2)
        strands = new_strands

    sequences = []
    structures = []
    data_strands = []
    for s in strands:
        data_strand = [round(data[x], 5) for x in s]
        data_strands.append(data_strand)
        sequences.append("".join([row["sequence"][x] for x in s]))
        structures.append("".join([row["structure"][x] for x in s]))
    return [
        "&".join(sequences),
        "&".join(structures),
        strands,
        data_strands,
    ]


def get_motifs(row, params: MotifSearchParams):
    s = rl.SecStruct(row["structure"], row["sequence"].replace("T", "U"))
    motifs = []
    for e in s:
        if e.id() < params.min_id:
            continue
        if e.id() > params.max_id:
            continue
        if params.mtype is not None and e.type() != params.mtype:
            continue
        if params.sequence is not None and e.sequence() != params.sequence:
            continue
        if params.structure is not None and e.structure() != params.structure:
            continue
        strands = e.strands()
        fail = False
        for s in strands:
            if s[-1] > params.max_pos:
                fail = True
            if s[0] < params.min_pos:
                fail = True
        if fail:
            continue
        motifs.append(e)
    return motifs


def get_motif_dataframe(
    df: pd.DataFrame, params: MotifSearchParams, mpos: int = -1, name="m"
):
    all_motifs = []
    all_data = []
    max_num = 0
    for i, row in df.iterrows():
        data = [row["name"]]
        motifs = get_motifs(row, params)
        if mpos != -1:
            motifs = [motifs[mpos]]
        all_motifs.append(motifs)
        for m in motifs:
            data.extend(_get_data_from_motif(m, row, params.extend))
        if len(motifs) > max_num:
            max_num = len(motifs)
        all_data.append(data)
    # setup columns
    all_cols = []
    if max_num > 1:
        for i in range(0, max_num):
            fname = name + "_" + str(i)
            cols = (
                f"{fname}_sequence,{fname}_structure,{fname}_ind,"
                f"{fname}_data".split(",")
            )
            all_cols.extend(cols)
    elif max_num == 1:
        fname = name
        all_cols = (
            f"{fname}_sequence,{fname}_structure,{fname}_ind,"
            f"{fname}_data".split(",")
        )
    else:
        raise ValueError(f"no motifs found with these params: {params}")
    all_cols.insert(0, "name")
    return pd.DataFrame(all_data, columns=all_cols)


# def get_ref_hp_avg(df, min_pos=0, max_pos=999):
#    df_ref = me.get_hairpin(df, "CGAGUAG", min_pos=min_pos, max_pos=max_pos)
#    return [(row[0][2] + row[0][5]) / 2 for row in df_ref["ref_hp_data"]]
