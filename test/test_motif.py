import pandas as pd
from dreem_tools.motif import *
from dreem_tools import paths

def test_get_motifs():
    path = paths.TEST_RESOURCES
    df = pd.read_json(path+"/summary.json")
    row = df.iloc[0]
    # check that can load all different types of motifs
    motifs = get_motifs(row)
    assert len(motifs) == 16
    juncs = get_motifs(row, MotifType.JUNCTION)
    assert len(juncs) == 3
    assert juncs[0].sequence() == "UAUG&CUAAG"
    ss_motifs = get_motifs(row, MotifType.SINGLESTRAND)
    assert len(ss_motifs) == 4
    hps = get_motifs(row, MotifType.HAIRPIN)
    assert len(hps) == 3
    # check if length constraints work




    #motifs = get_motifs()


