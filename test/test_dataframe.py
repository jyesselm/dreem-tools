import os
import pandas as pd
from dreem_tools import dataframe
from seq_tools.data_frame import convert_to_rna


def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return base_dir


LIB_PATH = get_lib_path()
TEST_RESOURCES = get_lib_path() + "/test/resources/"


def get_test_dataframe() -> pd.DataFrame:
    df = pd.read_json(TEST_RESOURCES + "summary.json")
    dataframe.trim(df, 21, 20)
    convert_to_rna(df)
    df = df[0:5]
    return df


def test_trim():
    df = pd.read_json(TEST_RESOURCES + "summary.json")
    dataframe.trim(df, 21, 20)
    assert len(df.loc[0]["sequence"]) == 109


def test_get_avg_reactivity():
    df = pd.read_json(TEST_RESOURCES + "summary.json")
    dataframe.trim(df, 21, 20)
    dataframe.get_avg_reactivity(df)
    assert "avg_reactivity" in df


def test_get_hairpin():
    df = get_test_dataframe()
    me = dataframe.MotifExtraction()
    df_m = me.get_hairpin(df, "CGAGUAG")
    assert len(df_m) == len(df)


def test_get_hairpin_extend():
    df = get_test_dataframe()
    me = dataframe.MotifExtraction()
    df_m = me.get_hairpin(df, "CGAGUAG", bp=2)
    assert df_m.loc[0]["hp_sequence"] == "CCCGAGUAGGG"
    assert df_m.loc[0]["hp_structure"] == "(((.....)))"
    assert len(df_m) == len(df)


def test_get_hairpins():
    df = get_test_dataframe()
    me = dataframe.MotifExtraction()
    df_m = me.get_hairpins(df)
    assert len(df_m) == 10


def test_get_twoway():
    df = get_test_dataframe()
    me = dataframe.MotifExtraction()
    df_m = me.get_twoway(df, "UAUG&CUAAG")
    assert len(df_m) == len(df)


def test_get_twoways():
    df = get_test_dataframe()
    me = dataframe.MotifExtraction()
    df_m = me.get_twoways(df)
    df_m = df_m[df_m["name"] == "01_seq_4787"]
    # should have 3 junctions, tetraloop-loop receptor and two ires
    assert list(df_m["junc_sequence"]) == [
        "UAUG&CUAAG",
        "GAACUAC&GC",
        "GC&GAACUAC",
    ]


def test_get_motif():
    df = get_test_dataframe()
    me = dataframe.MotifExtraction()
    df_m = me.get_motif(df, "AC")
    assert len(df_m) == len(df)


def test_get_motifs():
    pass
