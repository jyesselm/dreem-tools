import os
import pickle

from dreem_tools import util

def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return base_dir

def load(pathname):
    """Safe loader for binary pickle files"""
    fh = open(pathname, "rb")
    data = pickle.load(fh)
    fh.close()
    return data

LIB_PATH = get_lib_path()
TEST_RESOURCES = get_lib_path() + "/test/resources/"


def test_mutation_histos_to_dataframe():
    mhs = load(TEST_RESOURCES + "mutation_histos.p")
    df = util.mutation_histos_to_dataframe(mhs)
