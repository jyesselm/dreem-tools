import pandas as pd
import numpy as np
import rna_library as rl


def trim(df: pd.DataFrame, p5_dist: int, p3_dist: int) -> pd.DataFrame:
    df["sequence"] = [seq[p5_dist:-p3_dist] for seq in df["sequence"]]
    df["structure"] = [ss[p5_dist:-p3_dist] for ss in df["structure"]]
    df["data"] = [d[p5_dist:-p3_dist] for d in df["data"]]
    return df


def get_avg_reactivity(df: pd.DataFrame):
    avgs = []
    for i, row in df.iterrows():
        avg = 0.0
        count = 0.0
        for d, e in zip(row["data"], row["sequence"]):
            if e == "A" or e == "C":
                avg += d
                count += 1
        avgs.append(avg / count)
    df["avg_reactivity"] = avgs


def get_avg_noise(df: pd.DataFrame):
    avgs = []
    for i, row in df.iterrows():
        avg = 0.0
        count = 0.0
        for d, e in zip(row["data"], row["sequence"]):
            if e == "G" or e == "U" or e == "T":
                avg += d
                count += 1
        avgs.append(avg / count)
    df["avg_noise"] = avgs


def get_data_from_motif(m, row, extend):
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


def get_motifs(row, mtype=None, sequence=None, min_pos=0, max_pos=999):
    s = rl.SecStruct(row["structure"], row["sequence"].replace("T", "U"))
    motifs = []
    for e in s:
        if mtype is not None and e.type() != mtype:
            continue
        if sequence is not None and e.sequence() != sequence:
            continue
        strands = e.strands()
        fail = False
        for s in strands:
            if s[-1] > max_pos:
                fail = True
            if s[0] < min_pos:
                fail = True
        if fail:
            continue
        motifs.append(e)
    return motifs


def get_motifs_by_pos(df, mtype, mpos):
    data = []
    for i, row in df.iterrows():
        motifs = get_motifs(row, mtype)
        if len(motifs) <= mpos:
            raise ValueError(
                f"cannot get motif by pos if len is too short {len(motifs)}"
            )
        m = motifs[mpos]
        print(get_data_from_motif(m, row, 0))


class MotifExtraction(object):
    def __init__(self):
        pass

    def __get_motifs(
        self, row, mtype=None, sequence=None, min_pos=0, max_pos=999
    ):
        s = rl.SecStruct(row["structure"], row["sequence"].replace("T", "U"))
        motifs = []
        for e in s:
            if mtype is not None and e.type() != mtype:
                continue
            if sequence is not None and e.sequence() != sequence:
                continue
            strands = e.strands()
            fail = False
            for s in strands:
                if s[-1] > max_pos:
                    fail = True
                if s[0] < min_pos:
                    fail = True
            if fail:
                continue
            motifs.append(e)
        return motifs

    def __get_motif_dataframe(
        self,
        df,
        mtype=None,
        sequence=None,
        name="motif",
        extend=0,
        min_pos=0,
        max_pos=999,
    ):
        df_cols = "name,reads,aligned,sn,rna_name,dir,code".split(",")
        df_cols_act = []
        for col in df_cols:
            if col in df:
                df_cols_act.append(col)
        df_m_cols = (
            df_cols_act
            + f"{name}_name,{name}_sequence,{name}_structure,{name}_ind,{name}_data".split(
                ","
            )
        )
        all_data = []
        for i, row in df.iterrows():
            motifs = self.__get_motifs(row, mtype, sequence, min_pos, max_pos)
            for m in motifs:
                data = row[df_cols_act].tolist()
                data.append(name)
                data.extend(get_data_from_motif(m, row, extend))
                all_data.append(data)
        return pd.DataFrame(all_data, columns=df_m_cols)

    def get_hairpin(
        self,
        df: pd.DataFrame,
        sequence: str,
        name: str = "hp",
        bp: int = 0,
        min_pos: int = 0,
        max_pos: int = 999,
    ):
        return self.__get_motif_dataframe(
            df, rl.MotifType.HAIRPIN, sequence, name, bp, min_pos, max_pos
        )

    def get_hairpins(
        self,
        df: pd.DataFrame,
        name: str = "hp",
        bp: int = 0,
        min_pos: int = 0,
        max_pos: int = 999,
    ):
        return self.__get_motif_dataframe(
            df, rl.MotifType.HAIRPIN, None, name, bp, min_pos, max_pos
        )

    def get_twoway(
        self,
        df: pd.DataFrame,
        sequence: str,
        name: str = "junc",
        bp: int = 0,
        min_pos: int = 0,
        max_pos: int = 999,
    ):
        return self.__get_motif_dataframe(
            df, rl.MotifType.JUNCTION, sequence, name, bp, min_pos, max_pos
        )

    def get_twoways(
        self,
        df: pd.DataFrame,
        name: str = "junc",
        bp: int = 0,
        min_pos: int = 0,
        max_pos: int = 999,
    ):
        return self.__get_motif_dataframe(
            df, rl.MotifType.JUNCTION, None, name, bp, min_pos, max_pos
        )

    def get_motif(
        self,
        df: pd.DataFrame,
        sequence: str,
        name: str = "m",
        bp: int = 0,
        min_pos: int = 0,
        max_pos: int = 999,
    ):
        return self.__get_motif_dataframe(
            df, None, sequence, name, bp, min_pos, max_pos
        )

    def get_motifs(
        self,
        df: pd.DataFrame,
        name: str = "m",
        bp: int = 0,
        min_pos: int = 0,
        max_pos: int = 999,
    ):
        return self.__get_motif_dataframe(
            df, None, None, name, bp, min_pos, max_pos
        )


def get_ref_hp_avg(df, min_pos=0, max_pos=999):
    me = MotifExtraction()
    df_ref = me.get_hairpin(
        df, "CGAGUAG", name="ref_hp", min_pos=min_pos, max_pos=max_pos
    )
    return [(row[0][2] + row[0][5]) / 2 for row in df_ref["ref_hp_data"]]


def get_gaaa_avg(df, min_pos=0, max_pos=999):
    me = MotifExtraction()
    df_gaaa = me.get_hairpin(
        df, "GGAAAC", name="ttr_loop", min_pos=min_pos, max_pos=max_pos
    )
    return [
        (row[0][2] + row[0][2] + row[0][3]) / 3
        for row in df_gaaa["ttr_loop_data"]
    ]


def get_tlr_first_a(df, min_pos=0, max_pos=999):
    me = MotifExtraction()
    df_tltlr = me.get_twoway(
        df, "UAUG&CUAAG", name="tltlr", min_pos=min_pos, max_pos=max_pos
    )
    return [row[0][1] for row in df_tltlr["tltlr_data"]]
