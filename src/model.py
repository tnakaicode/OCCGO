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

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Lin
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_Trsf
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

sys.path.append(os.path.join('../'))
from src.base import plotocc

basepath = os.path.dirname(__file__) + "/"


class model_base (object):

    def __init__(self, meta={"name": "name"}):
        super().__init__()
        self.meta = meta
        self.name = meta["name"]

        if "axs" in meta.keys():
            pnt = gp_Pnt(*meta["axs"]["xyz"])
            dir_x = gp_Dir(*meta["axs"]["dir_x"])
            dir_y = gp_Dir(*meta["axs"]["dir_y"])
            dir_z = gp_Dir(*meta["axs"]["dir_z"])
            self.axs = gp_Ax3(pnt, dir_z, dir_x)
        else:
            self.axs = gp_Ax3()

        self.dat = np.loadtxt(basepath + "model_dat.txt")
        self.pts = []
        for xyz in self.dat:
            self.pts.append(gp_Pnt(*xyz))
        self.rim = make_polygon(self.pts, closed=True)


class model (plotocc):

    def __init__(self, cfgfile="./cfg/model.json"):
        super().__init__()
        self.rood_dir = basepath + "../"
        self.cfg = json.load(open(self.rood_dir + cfgfile, "r"))

    def set_model(self, name="surf"):
        if name not in self.cfg.keys():
            meta = {}
            meta["name"] = name
        else:
            meta = self.cfg[name]
        return model_base(meta)
