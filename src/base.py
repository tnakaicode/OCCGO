import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import pickle
import json
import time
import os
import glob
import shutil
import datetime
import platform
import subprocess
from optparse import OptionParser
from scipy.spatial import ConvexHull, Delaunay
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from PyQt5.QtWidgets import QApplication, qApp
from PyQt5.QtWidgets import QDialog, QCheckBox
from PyQt5.QtWidgets import QFileDialog
# pip install PyQt5

from OCC import VERSION
from OCC.Display.backend import load_backend, get_qt_modules
from OCC.Display.OCCViewer import OffscreenRenderer
log = logging.getLogger(__name__)

from OCC.Core.gp import gp_Cylinder, gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Pln, gp_Lin, gp_Elips
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_Trsf
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SHAPE, TopAbs_VERTEX, TopAbs_SHELL, TopAbs_SOLID
from OCC.Core.TopoDS import TopoDS_Iterator, topods_Vertex
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.Core.Bnd import Bnd_Box2d, Bnd_Box, Bnd_OBB
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepFill import BRepFill_CurveConstraint
from OCC.Core.BRepTools import breptools_Read
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepTools import breptools_Write
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepMeshData import BRepMeshData_Face
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.Geom import Geom_Plane, Geom_Surface, Geom_BSplineSurface
from OCC.Core.Geom import Geom_Curve, Geom_Line, Geom_Ellipse, Geom_Circle
from OCC.Core.Geom import Geom_RectangularTrimmedSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCC.Core.GeomFill import GeomFill_BoundWithSurf
from OCC.Core.GeomFill import GeomFill_BSplineCurves
from OCC.Core.GeomFill import GeomFill_StretchStyle, GeomFill_CoonsStyle, GeomFill_CurvedStyle
from OCC.Core.AIS import AIS_Manipulator
from OCC.Core.V3d import V3d_SpotLight, V3d_XnegYnegZpos
from OCC.Core.Quantity import Quantity_Color, Quantity_NOC_WHITE, Quantity_NOC_CORAL2, Quantity_NOC_BROWN
from OCC.Core.Graphic3d import Graphic3d_RM_RAYTRACING
from OCC.Core.BRepTools import breptools_Write, breptools_Read
from OCC.Core.StlAPI import StlAPI_Reader, StlAPI_Writer
from OCC.Extend.DataExchange import write_step_file, read_step_file
from OCC.Extend.DataExchange import write_iges_file, read_iges_file
from OCC.Extend.DataExchange import write_stl_file, read_stl_file
from OCCUtils.Topology import Topo
from OCCUtils.Topology import shapeTypeString, dumpTopology
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def pnt_to_xyz(p):
    return p.X(), p.Y(), p.Z()


def float_to_string(number):
    if number == 0 or abs(np.log10(abs(number))) < 100:
        return ' {: 0.10E}'.format(number)
    else:
        return ' {: 0.10E}'.format(number).replace('E', '')


def occ_to_grasp_cor(axs, name="name", filename="pln.cor"):
    pnt = axs.Location()
    v_x = axs.XDirection()
    v_y = axs.YDirection()
    fp = open(filename, "w")
    fp.write(' {:s}\n'.format(name))
    fp.write(' {:s}\n'.format("mm"))
    fp.write(''.join([float_to_string(v) for v in pnt_to_xyz(pnt)]) + '\n')
    fp.write(''.join([float_to_string(v) for v in pnt_to_xyz(v_x)]) + '\n')
    fp.write(''.join([float_to_string(v) for v in pnt_to_xyz(v_y)]) + '\n')
    fp.close()


def split_filename(filename="./temp_20200408000/not_ignore.txt"):
    name = os.path.basename(filename)
    rootname, ext_name = os.path.splitext(name)
    return name, rootname


def create_tempdir(name="temp", flag=1):
    print(datetime.date.today())
    datenm = "{0:%Y%m%d}".format(datetime.date.today())
    dirnum = len(glob.glob("./{}_{}*/".format(name, datenm)))
    if flag == -1 or dirnum == 0:
        tmpdir = "./{}_{}{:03}/".format(name, datenm, dirnum)
        os.makedirs(tmpdir)
        fp = open(tmpdir + "not_ignore.txt", "w")
        fp.close()
    else:
        tmpdir = "./{}_{}{:03}/".format(name, datenm, dirnum - 1)
    return tmpdir


def create_tempnum(name, tmpdir="./", ext=".tar.gz"):
    num = len(glob.glob(tmpdir + name + "*" + ext)) + 1
    filename = '{}{}_{:03}{}'.format(tmpdir, name, num, ext)
    #print(num, filename)
    return filename


class SetDir (object):

    def __init__(self, temp=True):
        self.root_dir = os.getcwd()
        self.tempname = ""
        self.rootname = ""

        pyfile = sys.argv[0]
        self.filename = os.path.basename(pyfile)
        self.rootname, ext_name = os.path.splitext(self.filename)

        if temp == True:
            self.create_tempdir()
            self.tempname = self.tmpdir + self.rootname
            print(self.rootname)

    def init(self):
        self.tempname = self.tmpdir + self.rootname

    def create_tempdir(self, name="temp", flag=1):
        self.tmpdir = create_tempdir(name, flag)
        self.tempname = self.tmpdir + self.rootname
        print(self.tmpdir)

    def create_dir(self, name="temp"):
        os.makedirs(name, exist_ok=True)
        if os.path.isdir(name):
            os.makedirs(name, exist_ok=True)
            fp = open(name + "not_ignore.txt", "w")
            fp.close()
            print("make {}".format(name))
        else:
            print("already exist {}".format(name))
        return name

    def create_dirnum(self, name="./temp", flag=+1):
        dirnum = len(glob.glob("{}_*/".format(name))) + flag
        if dirnum < 0:
            dirnum = 0
        dirname = name + "_{:03}/".format(dirnum)
        os.makedirs(dirname, exist_ok=True)
        fp = open(dirname + "not_ignore.txt", "w")
        fp.close()
        print("make {}".format(dirname))
        return dirname

    def add_tempdir(self, dirname="./", name="temp", flag=1):
        self.tmpdir = dirname
        self.tmpdir = create_tempdir(self.tmpdir + name, flag)
        self.tempname = self.tmpdir + self.rootname
        print(self.tmpdir)

    def add_dir(self, name="temp"):
        dirnum = len(glob.glob("{}/{}/".format(self.tmpdir, name)))
        if dirnum == 0:
            tmpdir = "{}/{}/".format(self.tmpdir, name)
            os.makedirs(tmpdir)
            fp = open(tmpdir + "not_ignore.txt", "w")
            fp.close()
            print("make {}".format(tmpdir))
        else:
            tmpdir = "{}/{}/".format(self.tmpdir, name)
            print("already exist {}".format(tmpdir))
        return tmpdir

    def add_dir_num(self, name="temp", flag=-1):
        if flag == -1:
            num = len(glob.glob("{}/{}_*".format(self.tmpdir, name))) + 1
        else:
            num = len(glob.glob("{}/{}_*".format(self.tmpdir, name)))
        tmpdir = "{}/{}_{:03}/".format(self.tmpdir, name, num)
        os.makedirs(tmpdir, exist_ok=True)
        fp = open(tmpdir + "not_ignore.txt", "w")
        fp.close()
        print("make {}".format(tmpdir))
        return tmpdir

    def open_filemanager(self, path="."):
        if sys.platform == "win32":
            subprocess.run('explorer.exe {}'.format(path))
        elif sys.platform == "linux":
            subprocess.check_call(['xdg-open', path])
        else:
            subprocess.run('explorer.exe {}'.format(path))

    def open_tempir(self):
        path = os.path.abspath(self.tmpdir)
        self.open_filemanager(path)

    def open_newtempir(self):
        self.create_tempdir("temp", -1)
        path = os.path.abspath(self.tmpdir)
        self.open_filemanager(path)

    def exit_app(self):
        sys.exit()


class PlotBase(SetDir):

    def __init__(self, aspect="equal", *args, **kwargs):
        SetDir.__init__(self)
        self.dim = 2
        self.fig, self.axs = plt.subplots()

    def new_fig(self, aspect="equal", dim=None):
        if dim == None:
            self.new_fig(aspect=aspect, dim=self.dim)
        elif self.dim == 2:
            self.new_2Dfig(aspect=aspect)
        elif self.dim == 3:
            self.new_3Dfig(aspect=aspect)
        else:
            self.new_2Dfig(aspect=aspect)

    def new_2Dfig(self, aspect="equal", *args, **kwargs):
        self.fig, self.axs = plt.subplots(*args, **kwargs)
        self.axs.set_aspect(aspect)
        self.axs.xaxis.grid()
        self.axs.yaxis.grid()
        return self.fig, self.axs

    def new_3Dfig(self, aspect="equal", *args, **kwargs):
        self.fig = plt.figure()
        self.axs = self.fig.add_subplot(111, projection='3d', *args, **kwargs)
        #self.axs = self.fig.gca(projection='3d')
        # self.axs.set_aspect('equal')

        self.axs.set_xlabel('x')
        self.axs.set_ylabel('y')
        self.axs.set_zlabel('z')

        self.axs.xaxis.grid()
        self.axs.yaxis.grid()
        self.axs.zaxis.grid()
        return self.fig, self.axs

    def SavePng(self, pngname=None, *args, **kwargs):
        if pngname == None:
            pngname = self.tempname + ".png"
        self.fig.savefig(pngname, *args, **kwargs)
        return pngname

    def SavePng_Serial(self, pngname=None, *args, **kwargs):
        if pngname == None:
            pngname = self.rootname
            dirname = self.tmpdir
        else:
            dirname = os.path.dirname(pngname) + "/"
            basename = os.path.basename(pngname)
            pngname, extname = os.path.splitext(basename)
        pngname = create_tempnum(pngname, dirname, ".png")
        self.fig.savefig(pngname, *args, **kwargs)
        return pngname

    def Show(self):
        try:
            plt.show()
        except AttributeError:
            pass


class plot2d (PlotBase):

    def __init__(self, aspect="equal", *args, **kwargs):
        PlotBase.__init__(self, *args, **kwargs)
        self.dim = 2
        # self.new_2Dfig(aspect=aspect)
        self.new_fig(aspect=aspect)

    def add_axs(self, row=1, col=1, num=1, aspect="auto"):
        self.axs.set_axis_off()
        axs = self.fig.add_subplot(row, col, num)
        axs.set_aspect(aspect)
        axs.xaxis.grid()
        axs.yaxis.grid()
        return axs

    def div_axs(self):
        self.div = make_axes_locatable(self.axs)
        # self.axs.set_aspect('equal')

        self.ax_x = self.div.append_axes(
            "bottom", 1.0, pad=0.5, sharex=self.axs)
        self.ax_x.xaxis.grid(True, zorder=0)
        self.ax_x.yaxis.grid(True, zorder=0)

        self.ax_y = self.div.append_axes(
            "right", 1.0, pad=0.5, sharey=self.axs)
        self.ax_y.xaxis.grid(True, zorder=0)
        self.ax_y.yaxis.grid(True, zorder=0)

    def sxy_to_nxy(self, mesh, sxy=[0, 0]):
        sx, sy = sxy
        nx, ny = mesh[0].shape
        xs, ys = mesh[0][0, 0], mesh[1][0, 0]
        xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
        dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
        mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
        return [mx, my]

    def contourf_sub(self, mesh, func, sxy=[0, 0], pngname=None):
        self.new_fig()
        self.div_axs()
        nx, ny = mesh[0].shape
        sx, sy = sxy
        xs, xe = mesh[0][0, 0], mesh[0][0, -1]
        ys, ye = mesh[1][0, 0], mesh[1][-1, 0]
        mx = np.searchsorted(mesh[0][0, :], sx) - 1
        my = np.searchsorted(mesh[1][:, 0], sy) - 1

        self.ax_x.plot(mesh[0][mx, :], func[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))
        self.ax_y.plot(func[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))
        im = self.axs.contourf(*mesh, func, cmap="jet")
        self.fig.colorbar(im, ax=self.axs, shrink=0.9)
        self.fig.tight_layout()
        if pngname == None:
            self.SavePng_Serial(pngname)
        else:
            self.SavePng(pngname)

    def contourf_sub2(self, mesh, func, sxy=[0, 0], pngname=None):
        self.new_fig()
        self.div_axs()
        nx, ny = mesh[0].shape
        sx, sy = sxy
        xs, xe = mesh[0][0, 0], mesh[0][0, -1]
        ys, ye = mesh[1][0, 0], mesh[1][-1, 0]
        mx = np.searchsorted(mesh[0][:, 0], sx) - 1
        my = np.searchsorted(mesh[1][0, :], sy) - 1

        self.ax_x.plot(mesh[0][mx, :], func[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))
        self.ax_y.plot(func[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))
        im = self.axs.contourf(*mesh, func, cmap="jet")
        self.fig.colorbar(im, ax=self.axs, shrink=0.9)
        self.fig.tight_layout()
        if pngname == None:
            self.SavePng_Serial(pngname)
        else:
            self.SavePng(pngname)

    def contourf_tri(self, x, y, z):
        self.new_fig()
        self.axs.tricontourf(x, y, z, cmap="jet")

    def contourf_sub_xy(self, mesh, func, sxy=[0, 0], pngname=None):
        self.new_fig()
        self.div_axs()
        nx, ny = mesh[0].shape
        sx, sy = sxy
        xs, xe = mesh[0][0, 0], mesh[0][0, -1]
        ys, ye = mesh[1][0, 0], mesh[1][-1, 0]
        mx = np.searchsorted(mesh[0][:, 0], sx) - 1
        my = np.searchsorted(mesh[1][0, :], sy) - 1

        self.ax_x.plot(mesh[0][mx, :], func[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))
        self.ax_y.plot(func[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))
        im = self.axs.contourf(*mesh, func, cmap="jet")
        self.fig.colorbar(im, ax=self.axs, shrink=0.9)
        self.fig.tight_layout()
        if pngname == None:
            self.SavePng_Serial(pngname)
        else:
            self.SavePng(pngname)

    def contourf_tri(self, x, y, z, lim=[-1, 1, -1, 1], pngname=None):
        self.new_2Dfig()
        self.axs.tricontourf(x, y, z, cmap="jet")
        self.axs.set_xlim(lim[0], lim[1])
        self.axs.set_ylim(lim[2], lim[3])

        if pngname == None:
            pngname = self.SavePng_Serial(pngname)
        else:
            pngname = self.SavePng(pngname)
        png_root, _ = os.path.splitext(pngname)
        self.axs.scatter(x, y, 5.0)
        self.SavePng(png_root + "_dot.png")

        pnt = np.array([x, y]).T
        cov = ConvexHull(pnt)
        tri = Delaunay(pnt)
        for idx in tri.simplices:
            xi = pnt[idx, 0]
            yi = pnt[idx, 1]
            self.axs.plot(xi, yi, "k", lw=0.5)
        self.axs.plot(pnt[cov.vertices, 0],
                      pnt[cov.vertices, 1], "k", lw=1.0)
        self.SavePng(png_root + "_grd.png")

        self.new_2Dfig()
        self.axs.scatter(x, y, 5.0)
        self.axs.set_xlim(lim[0], lim[1])
        self.axs.set_ylim(lim[2], lim[3])
        for idx in tri.simplices:
            xi = pnt[idx, 0]
            yi = pnt[idx, 1]
            self.axs.plot(xi, yi, "k", lw=0.75)
        self.axs.plot(pnt[cov.vertices, 0],
                      pnt[cov.vertices, 1], "k", lw=1.5)
        self.SavePng(png_root + "_grid.png")

    def contourf_div(self, mesh, func, loc=[0, 0], txt="", title="name", pngname=None, level=None):
        sx, sy = loc
        nx, ny = func.shape
        xs, ys = mesh[0][0, 0], mesh[1][0, 0]
        xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
        dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
        mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
        tx, ty = 1.1, 0.0

        self.new_2Dfig()
        self.div_axs()
        self.ax_x.plot(mesh[0][mx, :], func[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))

        self.ax_y.plot(func[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))

        self.fig.text(tx, ty, txt, transform=self.ax_x.transAxes)
        im = self.axs.contourf(*mesh, func, cmap="jet", levels=level)
        self.axs.set_title(title)
        self.fig.colorbar(im, ax=self.axs, shrink=0.9)

        plt.tight_layout()
        if pngname == None:
            self.SavePng_Serial(pngname)
        else:
            self.SavePng(pngname)

    def contourf_div_auto(self, mesh, func, loc=[0, 0], txt="", title="name", pngname=None, level=None):
        sx, sy = loc
        nx, ny = func.shape
        xs, ys = mesh[0][0, 0], mesh[1][0, 0]
        xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
        dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
        mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
        tx, ty = 1.1, 0.0

        self.div_axs()
        self.axs.set_aspect('auto')
        self.ax_x.plot(mesh[0][mx, :], func[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))

        self.ax_y.plot(func[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))

        self.fig.text(tx, ty, txt, transform=self.ax_x.transAxes)
        im = self.axs.contourf(*mesh, func, cmap="jet", levels=level)
        self.axs.set_title(title)
        self.fig.colorbar(im, ax=self.axs, shrink=0.9)

        plt.tight_layout()
        plt.savefig(pngname + ".png")


class plotpolar (plot2d):

    def __init__(self, aspect="equal"):
        plot2d.__init__(self)
        self.dim = 2
        self.new_polar(aspect)

    def new_polar(self, aspect="equal"):
        self.new_fig(aspect=aspect)
        self.axs.set_axis_off()
        self.axs = self.fig.add_subplot(111, projection='polar')

    def plot_polar(self, px, py, arrow=True, **kwargs):
        plt.polar(px, py, **kwargs)

        if arrow == True:
            num = np.linspace(1, len(px) - 1, 6)
            for idx, n in enumerate(num):
                n = int(n)
                plt.arrow(
                    px[-n], py[-n],
                    (px[-n + 1] - px[-n]) / 100.,
                    (py[-n + 1] - py[-n]) / 100.,
                    head_width=0.1,
                    head_length=0.2,
                )
                plt.text(px[-n], py[-n], "n={:d}".format(idx))


class plot3d (PlotBase):

    def __init__(self):
        PlotBase.__init__(self)
        self.dim = 3
        self.new_fig()

    def set_axes_equal(self):
        '''
        Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

        Input
          ax: a matplotlib axis, e.g., as output from plt.gca().
        '''

        x_limits = self.axs.get_xlim3d()
        y_limits = self.axs.get_ylim3d()
        z_limits = self.axs.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        y_range = abs(y_limits[1] - y_limits[0])
        z_range = abs(z_limits[1] - z_limits[0])

        x_middle = np.mean(x_limits)
        y_middle = np.mean(y_limits)
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        self.axs.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        self.axs.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        self.axs.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    def plot_ball(self, rxyz=[1, 1, 1]):
        u = np.linspace(0, 1, 10) * 2 * np.pi
        v = np.linspace(0, 1, 10) * np.pi
        uu, vv = np.meshgrid(u, v)
        x = rxyz[0] * np.cos(uu) * np.sin(vv)
        y = rxyz[1] * np.sin(uu) * np.sin(vv)
        z = rxyz[2] * np.cos(vv)

        self.axs.plot_wireframe(x, y, z)
        self.set_axes_equal()
        #self.axs.set_xlim3d(-10, 10)
        #self.axs.set_ylim3d(-10, 10)
        #self.axs.set_zlim3d(-10, 10)


def which(program):
    """Run the Unix which command in Python."""
    import os

    def is_exe(fpath):
        """Check if file is executable."""
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def _gmsh_path():
    """Find Gmsh."""

    if os.name == "nt":
        gmp = which("gmsh.exe")
    else:
        gmp = which("gmsh")
    if gmp is None:
        print(
            "Could not find Gmsh."
            + "Interactive plotting and shapes module not available."
        )
    return gmp


def pnt_from_axs(axs=gp_Ax3(), length=100):
    vec = point_to_vector(axs.Location()) + \
        dir_to_vec(axs.Direction()) * length
    return vector_to_point(vec)


def line_from_axs(axs=gp_Ax3(), length=100):
    return make_edge(axs.Location(), pnt_from_axs(axs, length))


def pnt_trf_vec(pnt=gp_Pnt(), vec=gp_Vec()):
    v = point_to_vector(pnt)
    v.Add(vec)
    return vector_to_point(v)


def set_trf(ax1=gp_Ax3(), ax2=gp_Ax3()):
    trf = gp_Trsf()
    trf.SetTransformation(ax2, ax1)
    return trf


def set_loc(ax1=gp_Ax3(), ax2=gp_Ax3()):
    trf = set_trf(ax1, ax2)
    loc = TopLoc_Location(trf)
    return loc


def trsf_scale(axs=gp_Ax3(), scale=1):
    trf = gp_Trsf()
    trf.SetDisplacement(gp_Ax3(), axs)
    return trf


def gen_ellipsoid(axs=gp_Ax3(), rxyz=[10, 20, 30]):
    sphere = BRepPrimAPI_MakeSphere(gp_Ax2(), 1).Solid()
    loc = set_loc(gp_Ax3(), axs)
    mat = gp_Mat(
        rxyz[0], 0, 0,
        0, rxyz[1], 0,
        0, 0, rxyz[2]
    )
    gtrf = gp_GTrsf(mat, gp_XYZ(0, 0, 0))
    ellips = BRepBuilderAPI_GTransform(sphere, gtrf).Shape()
    ellips.Location(loc)
    return ellips


def spl_face(px, py, pz, axs=gp_Ax3()):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

    api = GeomAPI_PointsToBSplineSurface(pnt_2d, 3, 8, GeomAbs_G2, 0.001)
    api.Interpolate(pnt_2d)
    #surface = BRepBuilderAPI_MakeFace(curve, 1e-6)
    # return surface.Face()
    face = BRepBuilderAPI_MakeFace(api.Surface(), 1e-6).Face()
    face.Location(set_loc(gp_Ax3(), axs))
    return face


def spl_curv(px, py, pz):
    num = px.size
    pts = []
    p_array = TColgp_Array1OfPnt(1, num)
    for idx, t in enumerate(px):
        x = px[idx]
        y = py[idx]
        z = pz[idx]
        pnt = gp_Pnt(x, y, z)
        pts.append(pnt)
        p_array.SetValue(idx + 1, pnt)
    api = GeomAPI_PointsToBSpline(p_array)
    return p_array, api.Curve()


def spl_curv_pts(pts=[gp_Pnt()]):
    num = len(pts)
    p_array = TColgp_Array1OfPnt(1, num)
    for idx, pnt in enumerate(pts):
        p_array.SetValue(idx + 1, pnt)
    api = GeomAPI_PointsToBSpline(p_array)
    return p_array, api.Curve()


def rotate_xyz(axs=gp_Ax3(), deg=0.0, xyz="x"):
    if xyz == "x":
        ax1 = gp_Ax1(axs.Location(), axs.XDirection())
    elif xyz == "y":
        ax1 = gp_Ax1(axs.Location(), axs.YDirection())
    elif xyz == "z":
        ax1 = gp_Ax1(axs.Location(), axs.Direction())
    else:
        ax1 = gp_Ax1(axs.Location(), axs.Direction())
    axs.Rotate(ax1, np.deg2rad(deg))


def rotate_axs(axs=gp_Ax3(), deg=0.0, idx="x"):
    ax = axs.Axis()
    if idx == "x":
        ax.SetDirection(axs.XDirection())
    elif idx == "y":
        ax.SetDirection(axs.YDirection())
    elif idx == "z":
        ax.SetDirection(axs.Direction())
    else:
        ax.SetDirection(axs.Direction())
    axs.Rotate(ax, np.deg2rad(deg))


def trf_axs(axs=gp_Ax3(), pxyz=[0, 0, 0], rxyz=[0, 0, 0]):
    rotate_axs(axs, rxyz[0], "x")
    rotate_axs(axs, rxyz[1], "y")
    rotate_axs(axs, rxyz[2], "z")
    axs.SetLocation(gp_Pnt(*pxyz))


def write_stl_file_mesh1(a_shape, filename, mode="ascii", linear_deflection=0.9, angular_deflection=0.5):
    """ 
    export the shape to a STL file
    Be careful, the shape first need to be explicitely meshed using BRepMesh_IncrementalMesh

    a_shape: the topods_shape to export
    filename: the filename
    mode: optional, "ascii" by default. Can either be "binary"
    linear_deflection: optional, default to 0.001. Lower, more occurate mesh
    angular_deflection: optional, default to 0.5. Lower, more accurate_mesh
    """
    if a_shape.IsNull():
        raise AssertionError("Shape is null.")
    if mode not in ["ascii", "binary"]:
        raise AssertionError("mode should be either ascii or binary")
    if os.path.isfile(filename):
        print("Warning: %s file already exists and will be replaced" % filename)
    # first mesh the shape
    mesh = BRepMesh_IncrementalMesh(
        a_shape, linear_deflection, True, angular_deflection, True)
    mesh.SetParallelDefault(True)
    mesh.Perform()
    if not mesh.IsDone():
        raise AssertionError("Mesh is not done.")

    stl_exporter = StlAPI_Writer()
    if mode == "ascii":
        stl_exporter.SetASCIIMode(True)
    else:  # binary, just set the ASCII flag to False
        stl_exporter.SetASCIIMode(False)
    stl_exporter.Write(a_shape, filename)

    if not os.path.isfile(filename):
        raise IOError("File not written to disk.")


def write_stl_file_mesh2(a_shape, filename, mode="ascii", linear_deflection=0.9, angular_deflection=0.5):
    """ 
    export the shape to a STL file
    Be careful, the shape first need to be explicitely meshed using BRepMesh_IncrementalMesh

    a_shape: the topods_shape to export
    filename: the filename
    mode: optional, "ascii" by default. Can either be "binary"
    linear_deflection: optional, default to 0.001. Lower, more occurate mesh
    angular_deflection: optional, default to 0.5. Lower, more accurate_mesh
    """
    if a_shape.IsNull():
        raise AssertionError("Shape is null.")
    if mode not in ["ascii", "binary"]:
        raise AssertionError("mode should be either ascii or binary")
    if os.path.isfile(filename):
        print("Warning: %s file already exists and will be replaced" % filename)
    # first mesh the shape
    mesh = BRepMesh_IncrementalMesh(
        a_shape, linear_deflection, False, angular_deflection, True)
    mesh.SetParallelDefault(True)
    mesh.Perform()
    if not mesh.IsDone():
        raise AssertionError("Mesh is not done.")

    stl_exporter = StlAPI_Writer()
    if mode == "ascii":
        stl_exporter.SetASCIIMode(True)
    else:  # binary, just set the ASCII flag to False
        stl_exporter.SetASCIIMode(False)
    stl_exporter.Write(a_shape, filename)

    if not os.path.isfile(filename):
        raise IOError("File not written to disk.")


class OCCViewer (object):

    def __init__(self):
        from OCC.Display.qtDisplay import qtViewer3d
        #self.app = self.get_app()
        #self.wi = self.app.topLevelWidgets()[0]
        self.vi = self.wi.findChild(qtViewer3d, "qt_viewer_3d")
        self.selected_shape = []

    def get_app(self):
        app = QApplication.instance()
        #app = qApp
        # checks if QApplication already exists
        if not app:
            app = QApplication(sys.argv)
        return app

    def on_select(self):
        self.vi.sig_topods_selected.connect(self._on_select)

    def clear_selected(self):
        self.selected_shape = []

    def _on_select(self, shapes):
        """
        Parameters
        ----------
        shape : TopoDS_Shape
        """
        for shape in shapes:
            print()
            print(shape.Location().Transformation())
            self.selected_shape.append(shape)
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


class GenCompound (object):

    def __init__(self):
        self.builder = BRep_Builder()
        self.compound = TopoDS_Compound()
        self.builder.MakeCompound(self.compound)

    def add_shapes(self, shps=[]):
        for shp in shps:
            self.builder.Add(self.compound, shp)


def check_callable(_callable):
    if not callable(_callable):
        raise AssertionError("The function supplied is not callable")


def init_qtdisplay(backend_str=None,
                   size=(1024, 768),
                   display_triedron=True,
                   background_gradient_color1=[206, 215, 222],
                   background_gradient_color2=[128, 128, 128]):
    """ This function loads and initialize a GUI using either wx, pyq4, pyqt5 or pyside.
    If ever the environment variable PYTHONOCC_OFFSCREEN_RENDERER, then the GUI is simply
    ignored and an offscreen renderer is returned.
    init_display returns 4 objects :
    * display : an instance of Viewer3d ;
    * start_display : a function (the GUI mainloop) ;
    * add_menu : a function that creates a menu in the GUI
    * add_function_to_menu : adds a menu option

    In case an offscreen renderer is returned, start_display and add_menu are ignored, i.e.
    an empty function is returned (named do_nothing). add_function_to_menu just execute the
    function taken as a paramter.

    Note : the offscreen renderer is used on the travis side.
    """
    if os.getenv("PYTHONOCC_OFFSCREEN_RENDERER") == "1":
        # create the offscreen renderer
        offscreen_renderer = OffscreenRenderer()

        def do_nothing(*kargs, **kwargs):
            """ takes as many parameters as you want,
            ans does nothing
            """
            pass

        def call_function(s, func):
            """ A function that calls another function.
            Helpfull to bypass add_function_to_menu. s should be a string
            """
            check_callable(func)
            log.info("Execute %s :: %s menu fonction" % (s, func.__name__))
            func()
            log.info("done")

        # returns empty classes and functions
        return offscreen_renderer, do_nothing, do_nothing, call_function
    used_backend = load_backend(backend_str)
    log.info("GUI backend set to: %s", used_backend)
    # Qt based simple GUI
    from OCC.Display.qtDisplay import qtViewer3d
    QtCore, QtGui, QtWidgets, QtOpenGL = get_qt_modules()

    class MainWindow(QtWidgets.QMainWindow):

        def __init__(self, *args):
            QtWidgets.QMainWindow.__init__(self, *args)
            self.canva = qtViewer3d(self)
            self.setWindowTitle(
                "pythonOCC-%s 3d viewer ('%s' backend)" % (VERSION, used_backend))
            self.setCentralWidget(self.canva)
            if sys.platform != 'darwin':
                self.menu_bar = self.menuBar()
            else:
                # create a parentless menubar
                # see: http://stackoverflow.com/questions/11375176/qmenubar-and-qmenu-doesnt-show-in-mac-os-x?lq=1
                # noticeable is that the menu ( alas ) is created in the
                # topleft of the screen, just
                # next to the apple icon
                # still does ugly things like showing the "Python" menu in
                # bold
                self.menu_bar = QtWidgets.QMenuBar()
            self._menus = {}
            self._menu_methods = {}
            # place the window in the center of the screen, at half the
            # screen size
            self.centerOnScreen()

        def centerOnScreen(self):
            '''Centers the window on the screen.'''
            resolution = QtWidgets.QApplication.desktop().screenGeometry()
            x = (resolution.width() - self.frameSize().width()) / 2
            y = (resolution.height() - self.frameSize().height()) / 2
            self.move(x, y)

        def add_menu(self, menu_name):
            _menu = self.menu_bar.addMenu("&" + menu_name)
            self._menus[menu_name] = _menu

        def add_function_to_menu(self, menu_name, _callable):
            check_callable(_callable)
            try:
                _action = QtWidgets.QAction(
                    _callable.__name__.replace('_', ' ').lower(), self)
                # if not, the "exit" action is now shown...
                _action.setMenuRole(QtWidgets.QAction.NoRole)
                _action.triggered.connect(_callable)
                self._menus[menu_name].addAction(_action)
            except KeyError:
                raise ValueError(
                    'the menu item %s does not exist' % menu_name)

    # following couple of lines is a tweak to enable ipython --gui='qt'
    # checks if QApplication already exists
    app = QtWidgets.QApplication.instance()
    if not app:  # create QApplication if it doesnt exist
        app = QtWidgets.QApplication(sys.argv)

    win = MainWindow()
    win.resize(size[0] - 1, size[1] - 1)
    win.show()
    win.centerOnScreen()
    win.canva.InitDriver()
    win.resize(size[0], size[1])
    win.canva.qApp = app
    display = win.canva._display

    def add_menu(*args, **kwargs):
        win.add_menu(*args, **kwargs)

    def add_function_to_menu(*args, **kwargs):
        win.add_function_to_menu(*args, **kwargs)

    def start_display():
        win.raise_()  # make the application float to the top
        app.exec_()

    if display_triedron:
        display.display_triedron()

    if background_gradient_color1 and background_gradient_color2:
        # background gradient
        display.set_bg_gradient_color(
            background_gradient_color1, background_gradient_color2)

    return display, start_display, add_menu, add_function_to_menu, win


class Viewer (object):

    def __init__(self):
        from OCC.Display.qtDisplay import qtViewer3d
        #self.app = self.get_app()
        #self.wis = self.app.topLevelWidgets()
        #self.wi = self.app.topLevelWidgets()[-1]
        self.vi = self.wi.findChild(qtViewer3d, "qt_viewer_3d")
        self.selected_shape = []

    def get_app(self):
        app = QApplication.instance()
        #app = qApp
        # checks if QApplication already exists
        if not app:
            app = QApplication(sys.argv)
        return app

    def on_select(self):
        self.vi.sig_topods_selected.connect(self._on_select)

    def clear_selected(self):
        self.selected_shape = []

    def _on_select(self, shapes):
        """
        Parameters
        ----------
        shape : TopoDS_Shape
        """
        for shape in shapes:
            print()
            self.selected_shape.append(shape)
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


class plotocc (SetDir, Viewer):

    def __init__(self, temp=True, disp=True, touch=False):
        SetDir.__init__(self, temp)
        self.base_axs = gp_Ax3()
        self.disp = disp
        self.touch = touch
        self.colors = ["BLUE", "RED", "GREEN",
                       "YELLOW", "BLACK", "WHITE", "BROWN"]

        # OCC Viewer
        if disp == True:
            self.display, self.start_display, self.add_menu, self.add_function, self.wi = init_qtdisplay()
            if touch == True:
                Viewer.__init__(self)
                self.on_select()

            self.SaveMenu()
            self.ViewMenu()
            self.SelectMenu()

        # GMSH
        self.gmsh = _gmsh_path()

    def show_box(self, axs=gp_Ax3(), lxyz=[100, 100, 100]):
        box = make_box(*lxyz)
        ax1 = gp_Ax3(
            gp_Pnt(-lxyz[0] / 2, -lxyz[1] / 2, -lxyz[2] / 2),
            gp_Dir(0, 0, 1)
        )
        trf = gp_Trsf()
        trf.SetTransformation(axs, gp_Ax3())
        trf.SetTransformation(ax1, gp_Ax3())
        box.Location(TopLoc_Location(trf))
        self.display.DisplayShape(axs.Location())
        self.show_axs_pln(axs, scale=lxyz[0])
        self.display.DisplayShape(box, transparency=0.7)

    def show_pnt(self, xyz=[0, 0, 0]):
        self.display.DisplayShape(gp_Pnt(*xyz))

    def show_pts(self, pts=[gp_Pnt()], num=1):
        for p in pts[::num]:
            self.display.DisplayShape(p)
        self.display.DisplayShape(make_polygon(pts))

    def show_ball(self, scale=100, trans=0.5):
        shape = BRepPrimAPI_MakeSphere(scale).Shape()
        self.display.DisplayShape(shape, transparency=trans)

    def show_vec(self, beam=gp_Ax3(), scale=1.0):
        pnt = beam.Location()
        vec = dir_to_vec(beam.Direction()).Scaled(scale)
        print(vec.Magnitude())
        self.display.DisplayVector(vec, pnt)

    def show_ellipsoid(self, axs=gp_Ax3(), rxyz=[10., 10., 10.], trans=0.5):
        shape = gen_ellipsoid(axs, rxyz)
        self.display.DisplayShape(shape, transparency=trans, color="BLUE")
        return shape

    def show_axs_vec(self, beam=gp_Ax3(), scale=1.0):
        pnt = beam.Location()
        vec = dir_to_vec(beam.Direction()).Scaled(scale)
        self.display.DisplayVector(vec, pnt)

    def show_axs_pln(self, axs=gp_Ax3(), scale=100, name=None):
        pnt = axs.Location()
        dx = axs.XDirection()
        dy = axs.YDirection()
        dz = axs.Direction()
        vx = dir_to_vec(dx).Scaled(1 * scale)
        vy = dir_to_vec(dy).Scaled(1 * scale)
        vz = dir_to_vec(dz).Scaled(1 * scale)

        pnt_x = pnt_trf_vec(pnt, vx)
        pnt_y = pnt_trf_vec(pnt, vy)
        pnt_z = pnt_trf_vec(pnt, vz)
        lx, ly, lz = make_line(pnt, pnt_x), make_line(
            pnt, pnt_y), make_line(pnt, pnt_z)
        self.display.DisplayShape(pnt)
        self.display.DisplayShape(lx, color="RED")
        self.display.DisplayShape(ly, color="GREEN")
        self.display.DisplayShape(lz, color="BLUE")
        if name != None:
            self.display.DisplayMessage(axs.Location(), name)
        return [lx, ly, lz]

    def show_plane(self, axs=gp_Ax3(), scale=100):
        pnt = axs.Location()
        vec = dir_to_vec(axs.Direction())
        pln = make_plane(pnt, vec, -scale, scale, -scale, scale)
        self.display.DisplayShape(pln)

    def show_wire(self, pts=[], axs=gp_Ax3()):
        poly = make_polygon(pts)
        poly.Location(set_loc(gp_Ax3(), axs))

        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        # n_sided.Build()
        #face = n_sided.Face()
        self.display.DisplayShape(poly)
        return poly

    def prop_axs(self, axs=gp_Ax3(), scale=100, xyz="z"):
        if xyz == "x":
            vec = dir_to_vec(axs.XDirection()).Scaled(scale)
        elif xyz == "y":
            vec = dir_to_vec(axs.YDirection()).Scaled(scale)
        elif xyz == "z":
            vec = dir_to_vec(axs.Direction()).Scaled(scale)
        else:
            vec = dir_to_vec(axs.Direction()).Scaled(scale)
        return axs.Translated(vec)

    def make_plane_axs(self, axs=gp_Ax3(), rx=[0, 500], ry=[0, 500]):
        pln = BRepBuilderAPI_MakeFace(
            gp_Pln(axs),
            rx[0], rx[1], ry[0], ry[1]
        ).Face()
        return pln

    def make_circle(self, axs=gp_Ax3(), radi=100):
        return make_wire(make_edge(Geom_Circle(axs.Ax2(), radi)))

    def make_EllipWire(self, rxy=[1.0, 1.0], shft=0.0, axs=gp_Ax3()):
        rx, ry = rxy
        if rx > ry:
            major_radi = rx
            minor_radi = ry
            axis = gp_Ax2()
            axis.SetXDirection(axis.XDirection())
        else:
            major_radi = ry
            minor_radi = rx
            axis = gp_Ax2()
            axis.SetXDirection(axis.YDirection())
        axis.Rotate(axis.Axis(), np.deg2rad(shft))
        elip = make_edge(gp_Elips(axis, major_radi, minor_radi))
        poly = make_wire(elip)
        poly.Location(set_loc(gp_Ax3(), axs))
        return poly

    def make_PolyWire(self, num=6, radi=1.0, shft=0.0, axs=gp_Ax3(), skin=None):
        lxy = radi
        pnts = []
        angl = 360 / num
        for i in range(num):
            thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            pnts.append(gp_Pnt(x, y, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)
        poly.Location(set_loc(gp_Ax3(), axs))

        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        n_sided.Build()
        face = n_sided.Face()
        if skin == None:
            return poly
        elif skin == 0:
            return face
        else:
            solid = BRepOffset_MakeOffset(
                face, skin, 1.0E-5, BRepOffset_Skin, False, True, GeomAbs_Arc, True, True)
            return solid.Shape()

    def make_StarWire(self, num=5, radi=[2.0, 1.0], shft=0.0, axs=gp_Ax3(), skin=None):
        lxy = radi
        pnts = []
        angl = 360 / num
        for i in range(num):
            a_thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            ax, ay = radi[0] * np.sin(a_thet), radi[0] * np.cos(a_thet)
            pnts.append(gp_Pnt(ax, ay, 0))
            b_thet = a_thet + np.deg2rad(angl) / 2
            bx, by = radi[1] * np.sin(b_thet), radi[1] * np.cos(b_thet)
            pnts.append(gp_Pnt(bx, by, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)
        poly.Location(set_loc(gp_Ax3(), axs))

        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        n_sided.Build()
        face = n_sided.Face()
        if skin == None:
            return poly
        elif skin == 0:
            return face
        else:
            solid = BRepOffset_MakeOffset(
                face, skin, 1.0E-5, BRepOffset_Skin, False, True, GeomAbs_Arc, True, True)
            return solid.Shape()

    def make_skin(self, pts=[], axs=gp_Ax3(), skin=1.0):
        poly = make_polygon(pts, closed=True)
        poly.Location(set_loc(gp_Ax3(), axs))

        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        n_sided.Build()
        face = n_sided.Face()
        solid = BRepOffset_MakeOffset(
            face, skin, 1.0E-5, BRepOffset_Skin, False, True, GeomAbs_Arc, True, True)
        return solid.Shape()

    def make_skin_wire(self, poly, axs=gp_Ax3(), skin=1.0):
        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        n_sided.Build()
        face = n_sided.Face()
        if skin == None:
            return poly
        elif skin == 0:
            return face
        else:
            solid = BRepOffset_MakeOffset(
                face, skin, 1.0E-5, BRepOffset_Skin, False, True, GeomAbs_Arc, True, True)
            return solid.Shape()

    def make_bnd_box(self, shp, axs=gp_Ax3()):
        bbox = Bnd_Box()
        # bbox.SetGap(1.0E-6)
        # bbox.Add(axs.Direction())
        brepbndlib_Add(shp, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        p0 = gp_Pnt(xmin, ymin, zmin)
        p1 = gp_Pnt(xmax, ymax, zmax)
        print(p0, p1)
        box = make_box(p0, p1)
        return box

    def make_FaceByOrder(self, pts=[]):
        dat = []
        for p in pts:
            dat.append([p.X(), p.Y(), p.Z()])

        dat = np.array(dat)
        cov = ConvexHull(dat, qhull_options='QJ')

        #pts_ord = []
        # print(cov)
        # print(cov.simplices)
        # print(cov.vertices)
        # for idx in cov.vertices:
        #    print(idx, pnt[idx])
        #    pts_ord.append(gp_Pnt(*pnt[idx]))

        #poly = make_polygon(pts_ord)
        poly = make_polygon(pts)
        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        n_sided.Build()
        face = n_sided.Face()
        return dat, face

    def make_trimmedcylinder(self, axs=gp_Ax1(), radii=700, hight=500, rng=[0, 0.1], xyz="y"):
        loc = self.prop_axs(axs, radii, "z")
        if xyz == "y":
            rotate_xyz(loc, deg=90, xyz="y")
        elif xyz == "x":
            rotate_xyz(loc, deg=90, xyz="x")
            rotate_xyz(loc, deg=-90, xyz="z")
        else:
            loc = self.prop_axs(loc, -radii, "z")
            loc = self.prop_axs(loc, -radii, "x")

        face = BRepBuilderAPI_MakeFace(
            gp_Cylinder(loc, radii),
            rng[0], rng[1],
            -hight / 2, hight / 2
        ).Face()
        return face

    def proj_rim_pln(self, wire, surf, axs=gp_Ax3()):
        proj = BRepProj_Projection(wire, surf, axs.Direction())
        return proj.Current()

    def proj_pnt_pln(self, pnt, surf, axs=gp_Ax3()):
        lin = gp_Lin(pnt, axs.Direction())
        api = GeomAPI_IntCS(Geom_Line(lin), BRep_Tool.Surface(surf))
        dst = np.inf
        sxy = pnt
        num = api.NbPoints()
        for i in range(num):
            p0 = api.Point(i + 1)
            dst0 = pnt.Distance(p0)
            if dst0 < dst:
                dst = dst0
                sxy = p0
        return sxy

    def proj_pln_show(self, face, nxy=[10, 10], ux=[0, 1], uy=[0, 1], axs=gp_Ax3()):
        trf = set_trf(gp_Ax3(), axs)
        pln = self.make_plane_axs(axs, [-1000, 1000], [-1000, 1000])
        surf = BRep_Tool.Surface(face)
        for px in np.linspace(ux[0], ux[1], nxy[0]):
            for py in np.linspace(uy[0], uy[1], nxy[1]):
                p0 = surf.Value(px, py)
                p1 = self.proj_pnt_pln(p0, pln, axs)
                self.display.DisplayShape(p0)
                self.display.DisplayShape(p1)

    def proj_pln_showup(self, face, nxy=[10, 10], lx=[-10, 10], ly=[-10, 10], axs=gp_Ax3()):
        trf = set_trf(gp_Ax3(), axs)
        nx, ny = nxy
        xs, xe = lx
        ys, ye = ly
        plnx = np.linspace(xs, xe, nx)
        plny = np.linspace(ys, ye, ny)
        mesh = np.meshgrid(plnx, plny)
        data = np.zeros_like(mesh[0])
        for (ix, iy), x in np.ndenumerate(data):
            px, py = mesh[0][ix, iy], mesh[1][ix, iy]
            p0 = gp_Pnt(px, py, 0).Transformed(trf)
            p1 = self.proj_pnt_pln(p0, face, axs)
            z = p0.Distance(p1)
            data[ix, iy] = z
            self.display.DisplayShape(p0)
            self.display.DisplayShape(p1)
        return mesh, data

    def AddManipulator(self):
        self.manip = AIS_Manipulator(self.base_axs.Ax2())
        ais_shp = self.display.DisplayShape(
            self.base_axs.Location(),
            update=True
        )
        self.manip.Attach(ais_shp)

    def SaveMenu(self):
        self.add_menu("File")
        self.add_function("File", self.import_cadfile)
        self.add_function("File", self.export_cap)
        self.add_function("File", self.display.FitAll)
        if self.touch == True:
            self.add_function("File", self.export_stp_selected)
            self.add_function("File", self.export_stl_selected)
            self.add_function("File", self.export_igs_selected)
            self.add_function("File", self.export_brep_selected)
            self.add_function("File", self.clear_selected)
        self.add_function("File", self.open_newtempir)
        self.add_function("File", self.open_tempir)
        self.add_function("File", self.exit_win)

    def ViewMenu(self):
        self.add_menu("View")
        self.add_function("View", self.display.FitAll)
        self.add_function("View", self.display.View_Top)  # XY-Plane(+)
        self.add_function("View", self.display.View_Bottom)  # XY-Plane(-)
        self.add_function("View", self.display.View_Rear)  # XZ-Plane(+)
        self.add_function("View", self.display.View_Front)  # XZ-Plane(-)
        self.add_function("View", self.display.View_Right)  # YZ-Plane(+)
        self.add_function("View", self.display.View_Left)  # YZ-Plane(-)
        self.add_function("View", self.ray_tracing_mode)
        self.add_function("View", self.display.SetRasterizationMode)

    def SelectMenu(self):
        self.add_menu("Select")
        self.add_function("Select", self.display.SetSelectionModeNeutral)
        self.add_function("Select", self.display.SetSelectionModeVertex)
        self.add_function("Select", self.display.SetSelectionModeEdge)
        self.add_function("Select", self.display.SetSelectionModeFace)
        self.add_function("Select", self.SetSelectionModeShape)

    def SetSelectionModeShape(self):
        self.display.SetSelectionMode(TopAbs_SHAPE)

    def ray_tracing_mode(self):
        # Raytracing Parameter
        # https://dev.opencascade.org/doc/occt-7.5.0/refman/html/class_graphic3d___c_view.html
        # Rendering Parameter List
        # https://dev.opencascade.org/doc/occt-7.5.0/refman/html/class_graphic3d___rendering_params.html
        # create one spotlight
        spot_light = V3d_SpotLight(
            gp_Pnt(-100, -100, 100), V3d_XnegYnegZpos, Quantity_Color(Quantity_NOC_WHITE))
        # display the spotlight in rasterized mode
        self.display.Viewer.AddLight(spot_light)
        self.display.View.SetLightOn()
        self.display.SetRaytracingMode(depth=8)

    def import_cadfile(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self.wi, 'QFileDialog.getOpenFileName()', '',
                                                  'CAD files (*.stp *.step *.stl *.igs *.iges, *.brep, *.geo)', options=options)
        self.read_cadfile(fileName)

    def read_geofile(self, geofile):
        # msh1, msh2, msh22, msh3, msh4, msh40, msh41, msh,
        # unv, vtk, wrl, mail, stl, p3d, mesh, bdf, cgns, med,
        # diff, ir3, inp, ply2, celum, su2, x3d, dat, neu, m, key
        filename = os.path.basename(geofile)
        rootname, ext_name = os.path.splitext(filename)
        geo_name = create_tempnum(self.rootname, self.tmpdir, ".geo")
        shutil.copyfile(geofile, geo_name)
        shutil.copyfile(geofile, self.tmpdir + filename)
        stl_name, _ = os.path.splitext(geo_name)
        stl_name += ".stl"
        txt_name, _ = os.path.splitext(geo_name)
        txt_name += "_gmsh.txt"

        gmsh_run = self.gmsh
        gmsh_run += " -tol 0.1"
        gmsh_run += " " + geo_name
        gmsh_run += " -3 -o"
        gmsh_run += " {} -format stl".format(stl_name)
        gmsh_run += " -log {}".format(txt_name)
        gmsh_success = os.system(gmsh_run)
        shpe = read_stl_file(stl_name)
        return shpe

    def read_cadfile(self, fileName, axs=gp_Ax3(), dsp=True):
        print(fileName)
        base_dir = os.path.dirname(fileName)
        basename = os.path.basename(fileName)
        rootname, extname = os.path.splitext(fileName)
        if extname in [".stp", ".step"]:
            shpe = read_step_file(fileName)
        elif extname in [".igs", ".iges"]:
            shpe = read_iges_file(fileName)
        elif extname in [".stl"]:
            shpe = read_stl_file(fileName)
        elif extname in [".brep"]:
            shpe = TopoDS_Shape()
            builder = BRep_Builder()
            breptools_Read(shpe, fileName, builder)
        elif extname in [".geo"]:
            shpe = self.read_geofile(fileName)
        else:
            print("Incorrect file index")
            # sys.exit(0)
        shpe.Location(set_loc(ax2=axs))
        if dsp == True:
            self.display.DisplayShape(shpe, update=True)
            self.export_cap()
        return shpe

    def export_cap(self):
        pngname = create_tempnum(self.rootname, self.tmpdir, ".png")
        self.display.View.Dump(pngname)

    def export_cap_name(self, pngname=None):
        if pngname == None:
            pngname = create_tempnum(self.rootname, self.tmpdir, ".png")
        self.display.View.Dump(pngname)

    def export_stp(self, shp):
        stpname = create_tempnum(self.rootname, self.tmpdir, ".stp")
        write_step_file(shp, stpname)

    def export_stp_selected(self):
        if self.touch == True:
            builder = BRep_Builder()
            compound = TopoDS_Compound()
            builder.MakeCompound(compound)
            for shp in self.selected_shape:
                print(shp)
                builder.Add(compound, shp)
            self.export_stp(compound)

    def export_stl(self, shp):
        stpname = create_tempnum(self.rootname, self.tmpdir, ".stl")
        write_stl_file(shp, stpname)

    def export_stl_selected(self):
        if self.touch == True:
            builder = BRep_Builder()
            compound = TopoDS_Compound()
            builder.MakeCompound(compound)
            for shp in self.selected_shape:
                print(shp)
                builder.Add(compound, shp)
            stlname = create_tempnum(self.rootname, self.tmpdir, ".stl")
            write_stl_file(compound, stlname)

    def export_igs_selected(self):
        if self.touch == True:
            builder = BRep_Builder()
            compound = TopoDS_Compound()
            builder.MakeCompound(compound)
            for shp in self.selected_shape:
                print(shp)
                builder.Add(compound, shp)
            igsname = create_tempnum(self.rootname, self.tmpdir, ".igs")
            write_iges_file(compound, igsname)

    def export_brep_selected(self):
        if self.touch == True:
            builder = BRep_Builder()
            compound = TopoDS_Compound()
            builder.MakeCompound(compound)
            for shp in self.selected_shape:
                print(shp)
                builder.Add(compound, shp)
            brepname = create_tempnum(self.rootname, self.tmpdir, ".brep")
            breptools_Write(compound, brepname)

    def export_stp_name(self, shp, stpname=None):
        if stpname == None:
            stpname = create_tempnum(self.rootname, self.tmpdir, ".stp")
        write_step_file(shp, stpname)

    def DumpJson(self, shp):
        jsonname = create_tempnum(self.rootname, self.tmpdir, ".json")
        fp = open(jsonname, "w")
        shp.DumpJson(fp)
        fp.close()

    def exit_win(self):
        self.wi.close()

    def show(self):
        self.display.FitAll()
        self.display.View.Dump(self.tempname + ".png")
        self.start_display()


class LineDrawer(object):

    def __init__(self, dirname="./tmp/", txtname="plot_data"):
        self.trajectory = None
        self.xx = []
        self.yy = []
        self.id = 0
        self.fg = 0

        self.dirname = dirname
        self.txtname = dirname + txtname
        self.fp = open(self.txtname + ".txt", "w")

        self.init_fig()

    def run_base(self):
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        animation.FuncAnimation(
            self.fig, self.anim_animate, init_func=self.anim_init, frames=30, interval=100, blit=True)

    def init_fig(self):
        self.fig, self.axs = plt.subplots()
        self.axs.set_aspect('equal')
        self.axs.xaxis.grid()
        self.axs.yaxis.grid()
        self.divider = make_axes_locatable(self.axs)

        self.traj_line, = self.axs.plot([], [], 'o', markersize=4, mew=4)
        self.record_line, = self.axs.plot(
            [], [], 'o', markersize=4, mew=4, color='m')
        self.empty, = self.axs.plot([], [])

    def onclick(self, event):
        txt = ""

        # get mouse position and scale appropriately to convert to (x,y)
        if event.xdata is not None:
            self.trajectory = np.array([event.xdata, event.ydata])
            txt += "event.x_f {:.3f} ".format(event.xdata)
            txt += "event.y_f {:.3f} ".format(event.ydata)

        if event.button == 1:
            self.id += 1
            txt += "event {:d} ".format(event.button)
            txt += "event.x_d {:d} ".format(event.x)
            txt += "event.y_d {:d} ".format(event.y)
            txt += "flag {:d} {:d}".format(self.id % 2, self.id)
            self.xx.append(event.xdata)
            self.yy.append(event.ydata)
        elif event.button == 3:
            dat_txt = "data-{:d} num {:d}\n".format(self.fg, self.id)
            for i in range(self.id):
                dat_txt += "{:d} ".format(i)
                dat_txt += "{:.3f} ".format(self.xx[i])
                dat_txt += "{:.3f} ".format(self.yy[i])
                dat_txt += "\n"

            self.fp.write(dat_txt)
            self.fp.write("\n\n")

            self.fq = open(self.txtname + "-{:d}.txt".format(self.fg), "w")
            self.fq.write(dat_txt)
            self.fq.close()

            self.xx = []
            self.yy = []
            self.id = 0
            self.fg += 1

        print(txt)

    def onkey(self, event):
        print("onkey", event, event.xdata)
        # Record
        if event.key == 'r':
            traj = np.array([self.xx, self.yy])
            with open('traj.pickle', 'w') as f:
                pickle.dump(traj, f)
                # f.close()

    def anim_init(self):
        self.traj_line.set_data([], [])
        self.record_line.set_data([], [])
        self.empty.set_data([], [])

        return self.traj_line, self.record_line, self.empty

    def anim_animate(self, i):
        if self.trajectory is not None:
            self.traj_line.set_data(self.trajectory)

        if self.xx is not None:
            self.record_line.set_data(self.xx, self.yy)

        self.empty.set_data([], [])

        return self.traj_line, self.record_line, self.empty

    def show(self):
        try:
            plt.show()
        except AttributeError:
            pass


if __name__ == '__main__':
    obj = SetDir()
    obj.create_tempdir(flag=-1)
    print(obj.basepath)
    # for i in range(5):
    #    name = "temp{:d}".format(i)
    #    obj.add_dir(name)
    # obj.add_dir(name)
    # obj.tmpdir = obj.add_dir("temp")
    # obj.add_dir("./temp/{}/".format(name))
