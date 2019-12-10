import numpy as np
import matplotlib.pyplot as plt
import sys, time, os
sys.path.append(os.path.join("../"))

from src.pyocc.OCCQt import Viewer
from src.pyocc.OCCDisplay import OCCDisplay

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pln
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS  import TopoDS_Compound
from OCC.Core.BOPAlgo import BOPAlgo_MakerVolume, BOPAlgo_Builder
from OCC.Core.BRep    import BRep_Builder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.Graphic3d   import Graphic3d_EF_PDF
from OCCUtils.Topology  import Topo
from OCCUtils.Construct import make_box, make_face
from OCCUtils.Construct import vec_to_dir

class OCCObject (OCCDisplay):

    def __init__(self):
        super(OCCObject, self).__init__()
        print (self.app, self.vi, self.wi)
    
    def Test(self):
        self.display.DisplayShape(make_box(100, 100, 100))
        self.display.DisplayShape(gp_Pnt())
        self.display.FitAll()
        self.start_display()
    
    def Cap (self):
        menu_name = "ScreenCapture"
        self.add_menu (menu_name)
        self.add_function_to_menu (menu_name, self.export_cap)

if __name__ == "__main__":
    obj = OCCObject()
    obj.Cap ()
    obj.Test()