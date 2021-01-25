import os

def get_jaspardb_path(fn,sub_dir=None):
    """
    Return a sqlite file from the pyjaspar data directory.
    This code is adapted from https://github.com/daler/pybedtools

    """
    #print(data_dir())
    #sys.exit()
    if sub_dir:
      fn = os.path.join(data_dir(), sub_dir, fn)
    else:
      fn = os.path.join(data_dir(), fn)
    #print(fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


def data_dir():
    """
    Returns the data directory that contains sqlite files.
    """
    #data_path = os.path.dirname(intervene.__file__)
    #data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'example_data')
    #print(data_path)
    return os.path.join(os.path.dirname(__file__), 'data')
