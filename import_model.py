import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import time
import os

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Lin

from src.base import create_tempnum
from src.model import model, model_base
from src.jsonfile import write_json

obj = model(cfgfile="./cfg/model.json")
print(obj.cfg["name"])

md0 = obj.set_model(name="surf0")
md1 = obj.set_model(name="surf1")
md2 = obj.set_model(name="surf")
md2.pts, md2.rim = obj.make_PolyWire(skin=None)
md3 = model_base(meta={"name": "surf3"})
md3.axs.Translate(gp_Pnt(), gp_Pnt(0, 0, 100))

meta = {}
meta["surf1"] = md1.export_dict()
meta["surf2"] = md2.export_dict()
meta["surf3"] = md3.export_dict()
write_json(create_tempnum("model", obj.tmpdir, ext=".json"), meta)

obj.display.DisplayShape(md1.rim)
obj.show()
