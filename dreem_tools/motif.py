import pandas as pd
import rna_library as rl

MotifType = rl.MotifType

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
