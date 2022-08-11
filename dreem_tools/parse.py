import pandas as pd
import copy

from dreem_tools import motif, dataframe


class ColumnParser(object):
    def __init__(self):
        pass

    def parse(self, df, params):
        params = copy.deepcopy(params)
        if "type" in params:
            return self.__parse_type(df, params)
        elif "pos" in params:
            pos = params["pos"]
            df_m = motif.get_motif_dataframe(
                df, motif.get_params_from_dict(params)
            )
            if "m_data" not in df_m:
                raise ValueError(
                    f"did not a single motif for column from {params}"
                )
            n_strand = df_m.iloc[0]["m_sequence"].count("&") + 1
            if n_strand != len(pos):
                raise ValueError(f"n_strand={n_strand} != len(pos)={len(pos)}")
            if len(pos) == 1:
                return [row[0][pos[0]] for row in df_m["m_data"]]
            elif len(pos) == 2:
                return [row[pos[0]][pos[1]] for row in df_m["m_data"]]
            elif len(pos) == 3:
                return [row[pos[0]][pos[1]][pos[2]] for row in df_m["m_data"]]
            else:
                raise ValueError(f"too many strands: {len(pos)}")
        else:
            raise ValueError("must contain member 'type' or 'pos'")

    def __parse_type(self, df, params):
        if params["type"] == "ref_hp_as":
            params["type"] = "ref_hp"
            df_m = motif.get_motif_dataframe(
                df, motif.get_params_from_dict(params)
            )
            return [(row[0][2] + row[0][5]) / 2 for row in df_m["m_data"]]
        elif params["type"] == "gaaa_avg":
            params["type"] = "gaaa_tetraloop"
            df_m = motif.get_motif_dataframe(
                df, motif.get_params_from_dict(params)
            )
            return [
                (row[0][2] + row[0][2] + row[0][3]) / 3
                for row in df_m["m_data"]
            ]
        elif params["type"] == "tlr_first_a":
            params["type"] = "tlr"
            df_m = motif.get_motif_dataframe(
                df, motif.get_params_from_dict(params)
            )
            return [row[0][1] for row in df_m["m_data"]]
        else:
            raise ValueError(f"unknown colum type: {params['type']}")


class MotifParser(object):
    def __init__(self):
        pass

    def parse(self, df, params):
        mpos = -1
        if "name" in params:
            name = params.pop("name")
        else:
            name = "m"
        if "mpos" in params:
            mpos = params.pop("mpos")
            df_m = motif.get_motif_dataframe(
                df, motif.get_params_from_dict(params), mpos=mpos, name=name
            )
        else:
            df_m = motif.get_motif_dataframe(
                df, motif.get_params_from_dict(params), name=name
            )
        # if requesting a specific motif there should only be one
        if mpos != -1 or "type" in params:
            col_name = f"{name}_sequence"
            if col_name not in df_m:
                raise ValueError(
                    f"did not get expected column name: {col_name} "
                )
        return df_m

    def __parse_type(self, df, params, name):
        df_m = motif.get_motif_dataframe(
            df, motif.get_params_from_dict(params), name=name
        )


def trim_exp_data(df, params):
    p5_trim = 0
    p3_trim = 0
    if "p5" in params["trim"]:
        p5_trim = params["trim"]["p5"]
    if "p3" in params["trim"]:
        p3_trim = params["trim"]["p3"]
    return dataframe.trim(df, p5_trim, p3_trim)


def parse_motif(df, params):
    mp = MotifParser()
    return mp.parse(df, params)


def parse_col(df, params):
    cp = ColumnParser()
    return cp.parse(df, params)


def parse_experiment(df, params):
    # remove any rows that do not have data
    mask = [True for _ in range(0, len(df))]
    for i, row in df.iterrows():
        if row["data"][0] is None:
            mask[i] = False
    df = df[mask]
    if "preset":
        pass
    if "trim" in params:
        df = trim_exp_data(df, params)
    if "cols" in params:
        col_params = params["cols"]
        for pname, cparams in col_params.items():
            df[pname] = parse_col(df, cparams)
    if "motifs" in params:
        m_params = params["motifs"]
        for pname, mparams in m_params.items():
            all_params = mparams.copy()
            all_params["name"] = pname
            df = df.merge(parse_motif(df, mparams))
    return df


def parse(df, params):
    dfs = []
    for exp_name, g in df.groupby("exp_name"):
        df_sub = pd.DataFrame(g[g["exp_name"] == exp_name])
        if exp_name in params:
            df_sub = parse_experiment(df_sub, params[exp_name])
        dfs.append(df_sub)
    df = pd.concat(dfs)
    df.to_json("processed.json", orient="records")
