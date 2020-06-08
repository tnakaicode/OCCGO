import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import scipy.constants as cnt


class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def write_json(filename="test.json", meta={}):
    json_str = json.dumps(meta, ensure_ascii=False, indent=4,
                          sort_keys=True, separators=(',', ': '), cls=NumpyEncoder)
    fp = open(filename, 'w')
    fp.write(json_str)
