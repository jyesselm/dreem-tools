import pandas as pd
from dacite import from_dict
from dreem_tools.motif import *
from dreem_tools import paths


def test_params():
    d = {"mtype": MotifType.JUNCTION}
    params = get_params_from_dict(d)
    assert MotifType.JUNCTION == params.mtype
    # test preset tests
    d = {"type": "ref_hp"}
    params = get_params_from_dict(d)
    assert params.sequence == "CGAGUAG"

def test_params_types():
    d = {"mtype": "twoway"}
    params = get_params_from_dict(d)


def test_get_motifs():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path + "/summary.json")
    row = df.iloc[0]
    # check that can load all different types of motifs
    motifs = get_motifs(row, MotifSearchParams())
    assert len(motifs) == 16
    params = MotifSearchParams(mtype=MotifType.JUNCTION)
    juncs = get_motifs(row, params)
    assert len(juncs) == 3
    assert juncs[0].sequence() == "UAUG&CUAAG"
    params = MotifSearchParams(mtype=MotifType.SINGLESTRAND)
    ss_motifs = get_motifs(row, params)
    assert len(ss_motifs) == 4
    params = MotifSearchParams(mtype=MotifType.HAIRPIN)
    hps = get_motifs(row, params)
    assert len(hps) == 3
    # check if length constraints work
    params = MotifSearchParams(min_pos=20, max_pos=100)
    motifs = get_motifs(row, params)
    assert len(motifs) == 5


def test_get_motifs_by_seq_and_ss():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path + "/summary.json")
    row = df.iloc[0]
    juncs = get_motifs(row, MotifSearchParams(sequence="UAUG&CUAAG"))
    assert len(juncs) == 1
    juncs = get_motifs(row, MotifSearchParams(sequence=""))
    assert len(juncs) == 0
    hps = get_motifs(row, MotifSearchParams(structure="(....)"))
    assert len(hps) == 2


def test_get_motif_dataframe():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path + "/summary.json")
    df = df[0:2]
    df_m = get_motif_dataframe(df, MotifSearchParams(mtype=MotifType.JUNCTION))
    assert len(df_m.columns) == 13
    assert "m_0_sequence" in df_m
