import sys
import time
import os
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopAbs import TopAbs_VERTEX
from OCC.Core.TopoDS import TopoDS_Iterator, topods_Vertex
from OCCUtils.Topology import shapeTypeString, dumpTopology
from PyQt5.QtWidgets import QApplication, qApp
from PyQt5.QtWidgets import QDialog, QCheckBox


class Viewer (object):

    def __init__(self):
        from OCC.Display.qtDisplay import qtViewer3d
        self.app = self.get_app()
        self.wi = self.app.topLevelWidgets()[0]
        self.vi = self.wi.findChild(qtViewer3d, "qt_viewer_3d")

    def get_app(self):
        app = QApplication.instance()
        #app = qApp
        # checks if QApplication already exists
        if not app:
            app = QApplication(sys.argv)
        return app

    def on_select(self):
        self.vi.sig_topods_selected.connect(self._on_select)

    def _on_select(self, shapes):
        """
        Parameters
        ----------
        shape : TopoDS_Shape
        """
        for shape in shapes:
            self.DumpTop(shape)

    def DumpTop(self, shape, level=0):
        """
        Print the details of an object from the top down
        """
        brt = BRep_Tool()
        s = shape.ShapeType()
        if s == TopAbs_VERTEX:
            pnt = brt.Pnt(topods_Vertex(shape))
            dmp = " " * level
            dmp += "%s - " % shapeTypeString(shape)
            dmp += "%.5e %.5e %.5e" % (pnt.X(), pnt.Y(), pnt.Z())
            print(dmp)
        else:
            dmp = " " * level
            dmp += shapeTypeString(shape)
            print(dmp)
        it = TopoDS_Iterator(shape)
        while it.More():
            shp = it.Value()
            it.Next()
            self.DumpTop(shp, level + 1)
