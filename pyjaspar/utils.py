import os
import csv
from io import StringIO

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

def dict_list_to_tsv(dict_list, keys, separator='\t'):
    # Function to convert list values to semi-colon separated strings
    def stringify_value(value):
        if isinstance(value, list):
            return ';'.join(map(str, value))
        return value

    # Preprocess the dictionary list
    dict_list = [{k: stringify_value(v) for k, v in d.items()} for d in dict_list]

    # Create a StringIO object to write CSV data to a string
    output_buffer = StringIO()
    writer = csv.DictWriter(output_buffer, fieldnames=keys, delimiter=separator)

    # Write header row
    writer.writeheader()

    # Write data rows
    for row in dict_list:
        writer.writerow(row)

    # Get the content of the StringIO buffer as a string
    output_text = output_buffer.getvalue()

    # Close the StringIO buffer
    output_buffer.close()

    return output_text