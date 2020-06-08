import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import time
import os
import shutil
import datetime
import platform
from optparse import OptionParser

basepath = os.path.dirname(__file__)


class Model (object):

    def __init__(self, cfgfile="./cfg/model.json"):
        super().__init__()
        self.rood_dir = basepath + "/../"
        self.cfg = json.load(open(self.rood_dir + cfgfile, "r"))
