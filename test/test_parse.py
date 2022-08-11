import pytest
import pandas as pd
from dreem_tools.parse import *
from dreem_tools import paths


def test_parse_col_types():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path + "/summary.json")
    df = pd.DataFrame(df[df["dir"] == "org_subset_36C_C006C_DMS"])
    params = {"type": "ref_hp_as"}
    data = parse_col(df, params)
    assert len(data) == len(df)
    params = {"type": "gaaa_avg"}
    data = parse_col(df, params)
    assert len(data) == len(df)
    params = {"type": "tlr_first_a"}
    data = parse_col(df, params)
    assert len(data) == len(df)


def test_parse_col_other():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path + "/summary.json")
    df = pd.DataFrame(df[df["dir"] == "org_subset_36C_C006C_DMS"])
    # get second A in tetraloop receptor
    params = {"sequence": "UAUG&CUAAG", "pos": [1, 1]}
    data = parse_col(df, params)
    assert len(data) == len(df)
    # catch wront pos size
    with pytest.raises(ValueError):
        params = {"sequence": "UAUG&CUAAG", "pos": [1]}
        data = parse_col(df, params)
    # get third A in GAAA
    params = {"sequence": "GGAAAC", "pos": [4]}
    data = parse_col(df, params)
    assert len(data) == len(df)
    # try getting more than one motif per row
    with pytest.raises(ValueError):
        params = {"mtype": "twoway", "pos": [0, 0]}
        data = parse_col(df, params)
    # didnt have 'type' or 'pos'
    with pytest.raises(ValueError):
        params = {"mtype": "twoway"}
        data = parse_col(df, params)


def test_parse_motif():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path + "/summary.json")
    df = pd.DataFrame(df[df["dir"] == "org_subset_36C_C006C_DMS"])
    params = {"type": "ref_hp"}
    df_m = parse_motif(df, params)
    # should have these columns
    assert "m_data" in df_m
    assert "m_sequence" in df_m
    assert "m_structure" in df_m
    assert "m_ind" in df_m

