import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy as sp
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
from scipy import ndimage
from scipy.spatial import ConvexHull, Delaunay
import argparse
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('parso').setLevel(logging.ERROR)


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

    def __init__(self):
        self.root_dir = os.getcwd()
        self.tempname = ""
        self.rootname = ""
        self.create_tempdir()

        pyfile = sys.argv[0]
        self.filename = os.path.basename(pyfile)
        self.rootname, ext_name = os.path.splitext(self.filename)
        self.tempname = self.tmpdir + self.rootname
        print(self.rootname)

    def init(self):
        self.tempname = self.tmpdir + self.rootname

    def create_tempdir(self, name="temp", flag=1):
        self.tmpdir = create_tempdir(name, flag)
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

    def open_file(self, filename=""):
        subprocess.Popen('explorer.exe {}'.format(filename))

    def open_filemanager(self, path="."):
        if sys.platform == "win32":
            subprocess.Popen('explorer.exe {}'.format(path))
        elif sys.platform == "linux":
            subprocess.check_call(['xdg-open', path])
        else:
            subprocess.run('explorer.exe {}'.format(path))

    def open_tempdir(self):
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


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


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

    def add_twin(self, aspect="auto", side="right", out=0):
        axt = self.axs.twinx()
        axt.set_aspect(aspect)
        axt.xaxis.grid()
        axt.yaxis.grid()
        axt.spines[side].set_position(('axes', out))
        make_patch_spines_invisible(axt)
        axt.spines[side].set_visible(True)
        return axt

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

    def contourf_sub_xy1(self, mesh, func, sxy=[0, 0], pngname=None):
        self.new_fig()
        self.div_axs()
        nx, ny = mesh[0].shape
        sx, sy = sxy
        xs, xe = mesh[0][0, 0], mesh[0][0, -1]
        ys, ye = mesh[1][0, 0], mesh[1][-1, 0]
        x_min, x_max = np.min(mesh[0][0, :]),np.max(mesh[0][0, :])
        y_min, y_max = np.min(mesh[1][:, 0]),np.max(mesh[1][:, 0])
        mx = np.searchsorted(mesh[0][:, 0], sx) - 1
        my = np.searchsorted(mesh[1][0, :], sy) - 1

        self.ax_x.plot(mesh[0][mx, :], func[mx, :])
        self.ax_x.set_title("y = {:.2f}".format(sy))
        self.ax_y.plot(func[:, my], mesh[1][:, my])
        self.ax_y.set_title("x = {:.2f}".format(sx))
        im = self.axs.contourf(*mesh, func, cmap="jet")
        
        #x_lin = np.linspace(x_min, x_max, 100)
        #y_lin = np.linspace(y_min, y_max, 100)
        #zi = ndimage.map_coordinates(func, np.vstack((x_lin, y_lin)))
        #
        #num_levels = len(im.allsegs)
        #num_element = len(im.allsegs[0])  # in level 0
        #num_vertices = len(im.allsegs[0][0])  # of element 0, in level 0
        #num_coord = len(im.allsegs[0][0][0])  # of vertex 0, in element 0, in level 0
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
        if pngname == None:
            self.SavePng_Serial(pngname)
        else:
            self.SavePng(pngname)


class plotpolar (plot2d):

    def __init__(self, aspect="equal", *args, **kwargs):
        plot2d.__init__(self, *args, **kwargs)
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

    def __init__(self, aspect="equal", *args, **kwargs):
        PlotBase.__init__(self, *args, **kwargs)
        self.dim = 3
        self.new_fig()

    def set_axes_equal(self, axis="xyz"):
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

        for i in axis:
            if i == "x":
                self.axs.set_xlim3d(
                    [x_middle - plot_radius, x_middle + plot_radius])
            elif i == "y":
                self.axs.set_ylim3d(
                    [y_middle - plot_radius, y_middle + plot_radius])
            elif i == "z":
                self.axs.set_zlim3d(
                    [z_middle - plot_radius, z_middle + plot_radius])
            else:
                self.axs.set_zlim3d(
                    [z_middle - plot_radius, z_middle + plot_radius])

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
    # for i in range(5):
    #    name = "temp{:d}".format(i)
    #    obj.add_dir(name)
    # obj.add_dir(name)
    # obj.tmpdir = obj.add_dir("temp")
    # obj.add_dir("./temp/{}/".format(name))
