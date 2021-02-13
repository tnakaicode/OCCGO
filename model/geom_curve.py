import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import urllib
import scipy.constants as cnt
from optparse import OptionParser

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Pln, gp_Trsf, gp_Lin
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.gp import gp_Circ, gp_Elips
from OCC.Core.Geom import Geom_BezierCurve, Geom_BSplineCurve
from OCC.Core.Geom import Geom_Circle, Geom_Ellipse, Geom_Curve
from OCC.Core.Geom import Geom_BSplineSurface, Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d
from OCC.Core.TColGeom import TColGeom_Array1OfCurve
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire


def curvature(px, r=1000, s=0.0, t=10.0):
    """( x + sx )**2 / 2*rx + tx"""
    if (r == 0):
        py = np.zeros_like(px + s) + t
    else:
        py = (px + s)**2 / (2*r) + t
    return py


def pnt_eclips(pnt, wxy=[10, 20]):
    if wxy[0] >= wxy[1]:
        ax2 = gp_Ax2(pnt, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
        w_x = wxy[0]
        w_y = wxy[1]
    elif wxy[1] >= wxy[0]:
        ax2 = gp_Ax2(pnt, gp_Dir(0, 0, 1), gp_Dir(0, 1, 0))
        w_x = wxy[1]
        w_y = wxy[0]
    else:
        ax2 = gp_Ax2(pnt, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
        w_x = wxy[0]
        w_y = wxy[1]
    return Geom_Ellipse(ax2, w_x, w_y)


def get_airfol_data(filename="./dae51.dat"):
    fp = urllib.request.urlopen(filename)
    dat = fp.readlines()
    print(filename)
    print(dat[0].split()[0])
    print(dat[1].split())
    n0 = int(float(dat[1].split()[0]))
    n1 = int(float(dat[1].split()[1]))

    ns0 = 3
    ne0 = ns0+n0
    ns1 = 3+n0+1
    ne1 = ns1+n1
    upp = []
    bot = []
    for line in dat[ns0:ne0]:
        x, y = [float(v) for v in line.split()]
        upp.append(np.array([x, y, 0]))

    for line in dat[ns1:ne1]:
        x, y = [float(v) for v in line.split()]
        bot.append(np.array([x, y, 0]))
    return np.array(upp), np.array(bot)


def gen_data_spline(dat):
    num = dat.shape[0]
    pts = TColgp_Array1OfPnt(1, num)
    for i, xyz in enumerate(dat):
        pnt = gp_Pnt(*xyz)
        pts.SetValue(i+1, pnt)
    geo_spl = GeomAPI_PointsToBSpline(pts, 4)
    return geo_spl.Curve()


class Airfoil (object):

    def __init__(self, url):
        self.display, self.start_display, self.add_menu, self.add_function_to_menu = init_display()
        self.url = url

    def gen_circle(self, pz, pr0, pr1):
        num = pz.shape[0]
        self.crv = TColGeom_Array1OfCurve(1, num)
        for i, z in enumerate(pz):
            pnt = gp_Pnt(0, 0, z)
            wxy = [np.abs(pr0[i]), np.abs(pr1[i])]
            ecl = pnt_eclips(pnt, wxy)
            self.crv.SetValue(i+1, ecl)
            self.display.DisplayShape(ecl)

    def gen_through(self):
        obj = BRepOffsetAPI_ThruSections()

        ax2_1 = gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
        crl_1 = gp_Circ(ax2_1, 100)
        obj.AddWire(crl_1)

        ax2_2 = gp_Ax2(gp_Pnt(0, 0, 100), gp_Dir(0, 0, 1))
        crl_2 = gp_Circ(ax2_2, 200)
        obj.AddWire(crl_2)

        obj.Build()
        self.display.DisplayShape(obj.Shape())

    def get_airfoil(self, name="dae51"):
        filename = self.url + name + ".dat"
        self.upp, self.bot = get_airfol_data(filename)
        self.upp_spl = self.gen_data_spline(self.upp)
        self.bot_spl = self.gen_data_spline(self.bot)

    def gen_data_spline(self, dat):
        num = dat.shape[0]
        pts = TColgp_Array1OfPnt(1, num)
        for i, xyz in enumerate(dat):
            pnt = gp_Pnt(*xyz)
            pts.SetValue(i+1, pnt)
        geo_spl = GeomAPI_PointsToBSpline(pts)
        return geo_spl

    def OCC_Display(self):
        for i, xyz in enumerate(self.upp):
            pnt = gp_Pnt(*xyz)
            self.display.DisplayShape(pnt)

        for i, xyz in enumerate(self.bot):
            pnt = gp_Pnt(*xyz)
            self.display.DisplayShape(pnt)

        self.display.DisplayShape(self.upp_spl.Curve())
        self.display.DisplayShape(self.bot_spl.Curve())

        self.display.FitAll()
        # self.start_display()


if __name__ == "__main__":
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="dae51")
    opt, argc = parser.parse_args(argvs)
    print(argc, opt)

    cfg = json.load(open("airfoil.json", "r"))
    airfilo_url = cfg["url"]

    obj = Airfoil(airfilo_url)
    # obj.get_airfoil("dae51")

    """
    y = x^2 / 2r
    x = sqrt (2r * y)
    r = x^2 / 2y
    """

    pz = np.linspace(0, 1000, 100)
    p0 = curvature(pz, r=-(750**2/(2*500)), s=-500, t=500)
    p1 = curvature(pz, r=-(750**2/(2*750)), s=-500, t=750)
    #obj.gen_circle(pz, p0, p1)
    obj.gen_through()

    obj.start_display()
