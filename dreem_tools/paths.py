import os


def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return base_dir


LIB_PATH = get_lib_path()
TEST_RESOURCES = get_lib_path() + "/test/resources/"
