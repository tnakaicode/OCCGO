import numpy as np
import matplotlib.pyplot as plt
import json
import sys
import time
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from linecache import getline, clearcache
from scipy.integrate import simps
from scipy.special import erf
from scipy.special import eval_hermite
from math import factorial as fact
from unwrap.unwrap import unwrap

from .geomtory import grad_pln, curvature
from .outer import v_outer, v_normal


def integrate_simps(mesh, func):
    nx, ny = func.shape
    px, py = mesh[0][int(nx / 2), :], mesh[1][:, int(ny / 2)]
    val = simps(simps(func, px), py)
    return val


def normalize_integrate(mesh, func):
    return func / integrate_simps(mesh, func)


def normalize_max(mesh, func):
    return func / func.max()


def moment(mesh, func, index):
    ix, iy = index[0], index[1]
    g_func = normalize_integrate(mesh, func)
    fxy = g_func * mesh[0]**ix * mesh[1]**iy
    val = integrate_simps(mesh, fxy)
    return val


def moment_seq(mesh, func, num):
    seq = np.empty([num, num])
    for ix in range(num):
        for iy in range(num):
            seq[ix, iy] = moment(mesh, func, [ix, iy])
    return seq


def get_centroid(mesh, func):
    dx = moment(mesh, func, (1, 0))
    dy = moment(mesh, func, (0, 1))
    return dx, dy


def get_weight(mesh, func, dxy):
    g_mesh = np.array((mesh[0] - dxy[0], mesh[1] - dxy[1]))
    lx = moment(g_mesh, func, (2, 0))
    ly = moment(g_mesh, func, (0, 2))
    return np.sqrt(lx) * np.sqrt(2), np.sqrt(ly) * np.sqrt(2)


def get_cov(mesh, func, dxy):
    g_mesh = [mesh[0] - dxy[0], mesh[1] - dxy[1]]
    Mxx = moment(g_mesh, func, (2, 0))
    Myy = moment(g_mesh, func, (0, 2))
    Mxy = moment(g_mesh, func, (1, 1))
    mat = np.array([
        [Mxx, Mxy],
        [Mxy, Myy],
    ])
    return mat


def get_wxy(mesh, func):
    sxy = get_centroid(mesh, func)
    cov = get_cov(mesh, func, sxy)
    wxy, mat = np.linalg.eig(cov)
    wxy = np.sqrt(wxy)
    rot = np.arcsin(mat[0, 1])
    #print (wxy, np.rad2deg(-np.arcsin(mat[0, 1])), np.rad2deg(np.arccos(mat[0, 0])))
    """if wxy[0] > wxy[1]:
        rot = np.pi/2 - (-np.arcsin(mat[0, 1]))
        wxy = wxy[1], wxy[0]
    else:
        rot = -np.arcsin(mat[0, 1]) 
        wxy = wxy[0], wxy[1]"""
    g_func = gaussian_func(mesh, sxy, wxy, rot)
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    g_mesh = [px, py]
    for i in range(20):
        wx = np.sqrt(2) * np.sqrt(integrate_simps(g_mesh,
                                                  g_mesh[0]**2 * g_func * func) / integrate_simps(g_mesh, g_func * func))
        wy = np.sqrt(2) * np.sqrt(integrate_simps(g_mesh,
                                                  g_mesh[1]**2 * g_func * func) / integrate_simps(g_mesh, g_func * func))
        wxy = [wx, wy]
        #print (wxy, gcf_calc(mesh, func, g_func))
        g_func = gaussian_func(mesh, sxy, wxy, rot)
    #print (wxy, np.rad2deg(rot))
    return wxy, rot


def profile_out(meta, ampl, filename, dir_name="./"):
    mesh = meta["MESH"]
    seq_num = 5

    fp = open(dir_name + filename, "w")
    pw = integrate_simps(mesh, ampl)
    sxy = get_centroid(mesh, ampl)
    wxy, rot = get_wxy(mesh, ampl)
    seq = moment_seq(mesh, ampl, seq_num)

    fp.write('{:5d} {:5d}\n'.format(meta["NX"], meta["NY"]))
    fp.write('{:0.5E} {:0.5E}\n'.format(meta["XS"], meta["XE"]))
    fp.write('{:0.5E} {:0.5E}\n'.format(meta["YS"], meta["YE"]))
    fp.write("\n")
    fp.write('{} \n'.format(meta["NAME"]))
    fp.write('pw {:0.5E} MeW\n'.format(pw))
    fp.write('sx {:0.5E} mm \n'.format(sxy[0]))
    fp.write('sy {:0.5E} mm \n'.format(sxy[1]))
    fp.write('wx {:0.5E} mm \n'.format(wxy[0]))
    fp.write('wy {:0.5E} mm \n'.format(wxy[1]))
    fp.write("\n")
    for ix in range(seq_num):
        for iy in range(seq_num):
            fp.write('{:d} {:d} {:0.5E} \n'.format(ix, iy, seq[ix, iy]))


def rot_mesh(mesh, sxy=[0, 0], rot=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    return [px, py]


def gauss_1d(px, sx=0, wx=10):
    py = np.exp(-0.5 * ((px - sx) / wx)**2)
    return py


def cdf_1d(px, sx=0, wx=10, kx=2):
    return (1 + erf(kx * (px - sx) / wx / np.sqrt(wx))) / 2


def gauss_1d_skew(px, sx=0, wx=10, kx=2):
    py = gauss_1d(px, sx, wx)
    py *= cdf_1d(px, sx, wx, kx)
    return py


def gauss_2d(mesh, sxy=[0, 0], wxy=[10, 10], deg=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    rot = np.deg2rad(deg)
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = np.exp(-0.5 * (px / wxy[0])**2)
    fy = np.exp(-0.5 * (py / wxy[1])**2)
    return fx * fy


def gauss_2d_skew(mesh, sxy=[0, 0], wxy=[10, 10], kxy=[2, 2], deg=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    rot = np.deg2rad(deg)
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = gauss_1d_skew(px, 0, wxy[0], kxy[0])
    fy = gauss_1d_skew(py, 0, wxy[1], kxy[1])
    return fx * fy


def gaussian_func(mesh, sxy=[0, 0], wxy=[10, 10], rot=0.0):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = np.exp(-0.5 * (px / wxy[0])**2)
    fy = np.exp(-0.5 * (py / wxy[1])**2)
    return fx * fy


def hermite_gauss_1d(px, p, w):
    co = 1 / np.sqrt(2**p * fact(p) * w * np.sqrt(np.pi / 2))
    hx = co * eval_hermite(p, np.sqrt(2) / w * px) * np.exp(-0.5 * (px / w)**2)
    return hx


def hermite_gauss_func(mesh, sxy=[0, 0], wxy=[10, 10], rot=0.0, pxy=[0, 0]):
    x, y = mesh[0] - sxy[0], mesh[1] - sxy[1]
    px = x * np.cos(rot) - y * np.sin(rot)
    py = y * np.cos(rot) + x * np.sin(rot)
    fx = hermite_gauss_1d(px, pxy[0], wxy[0])
    fy = hermite_gauss_1d(py, pxy[1], wxy[1])
    return fx * fy


def mask_eclip(mesh, sxy=[0, 0], rxy=[10, 10]):
    func = np.zeros_like(mesh[0])
    nx, ny = mesh[0].shape
    sx, sy = sxy
    rx, ry = rxy
    for ix in range(nx):
        for iy in range(ny):
            x, y = mesh[0][ix, iy] - sx, mesh[1][ix, iy] - sy
            if ((x / rx)**2 + (y / ry)**2) < 1:
                func[ix, iy] = 1.0
            else:
                func[ix, iy] = 0.0
    return func


def offset_phas(mesh, phas, sxy=[0, 0]):
    dx, dy = grad_pln(mesh, phas, sxy)
    g_phas = np.tan(dx) * mesh[0] + np.tan(dy) * mesh[1]
    g_func = np.ones_like(mesh[0]) * np.exp(-1j * (phas - g_phas))
    # return unwrap(np.angle(g_func))
    # return phas - g_phas
    return np.angle(g_func)


def offset_phas_deg(mesh, phas, dxy=[0, 0]):
    dx, dy = dxy
    g_phas = np.tan(dx) * mesh[0] + np.tan(dy) * mesh[1]
    g_func = np.ones_like(mesh[0]) * np.exp(-1j * (phas - g_phas))
    return np.angle(g_func)
    # return phas-g_phas


def gcf_calc(mesh, func, g_func):
    e_rg = integrate_simps(mesh, func * g_func)
    e_r2 = integrate_simps(mesh, np.abs(func)**2)
    e_g2 = integrate_simps(mesh, np.abs(g_func)**2)
    return np.abs(e_rg) / np.sqrt(e_r2 * e_g2)


def coupling_calc(mesh, func, g_func):
    e_rg = integrate_simps(mesh, func * g_func)
    e_r2 = integrate_simps(mesh, func * np.conjugate(func))
    e_g2 = integrate_simps(mesh, g_func * np.conjugate(g_func))
    return np.abs(e_rg)**2 / np.abs(e_r2 * e_g2)


def gcf_itr(mesh, func, g_ampl, dxy=[0, 0], mxy=[5, 5], scale=1.0):
    dx, dy = dxy
    mx, my = mxy
    g_gcf2 = np.empty([3, mx * my])
    g_dx = np.linspace(-1, 1, mx) * scale + np.rad2deg(dx)
    g_dy = np.linspace(-1, 1, my) * scale + np.rad2deg(dy)
    for ix, dgx in enumerate(g_dx):
        for iy, dgy in enumerate(g_dy):
            g_phas = np.tan(np.deg2rad(dgx)) * \
                mesh[0] + np.tan(np.deg2rad(dgy)) * mesh[1]
            g_func = g_ampl * np.exp(-1j * g_phas)
            g_gcf2[:, ix * my +
                   iy] = np.deg2rad(dgx), np.deg2rad(dgy), gcf_calc(mesh, func, g_func)
    iz = np.argmax(g_gcf2[2])
    dx, dy, v = g_gcf2[:, iz]
    print(np.rad2deg(dx), np.rad2deg(dy), v)
    return dx, dy


def get_phas_prof(mesh, surf, sxy=[0, 0], wxy=[10, 10]):
    sx, sy = sxy
    wx, wy = wxy
    nx, ny = surf.shape
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    nx_c, ny_c = int((sy - ys) / dy), int((sx - xs) / dx)
    nx_w, ny_w = int(wy / dy / 2), int(wx / dx / 2)
    nx_s, nx_e = nx_c - nx_w, nx_c + nx_w
    ny_s, ny_e = ny_c - ny_w, ny_c + ny_w
    X = mesh[0][nx_s:nx_e, ny_s:ny_e].flatten()
    Y = mesh[1][nx_s:nx_e, ny_s:ny_e].flatten()
    Z = surf[nx_s:nx_e, ny_s:ny_e].flatten()
    A = [np.ones_like(X)]
    for i in range(1, 3):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))

    A = np.array(A).T
    coeff, r, rank, s = np.linalg.lstsq(A, Z, rcond=-1)
    #print (coeff)
    #print (r)
    #print (s)
    dxy = np.arctan(coeff[2]) % (np.pi / 2), np.arctan(coeff[1]) % (np.pi / 2)
    #dxy = np.arctan(coeff[2]), np.arctan(coeff[1])
    if coeff[5] == 0:
        rx = 0
    else:
        rx = 1 / coeff[5] / 2
    if coeff[3] == 0:
        ry = 0
    else:
        ry = 1 / coeff[3] / 2
    rxy = rx, ry
    return dxy, rxy


def get_kvector(mesh, phas, sxy=[0, 0], wxy=[10, 10]):
    sx, sy = sxy
    wx, wy = wxy
    nx, ny = phas.shape
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    d, rdi = get_phas_prof(mesh, phas, sxy, wxy)
    surf = phas - curvature(mesh[0], rdi[0], sx) - \
        curvature(mesh[1], rdi[1], sy)
    nx_c, ny_c = int((sy - ys) / dy), int((sx - xs) / dx)
    grdx = np.gradient(mesh[0], axis=1)
    grdy = np.gradient(mesh[1], axis=0)
    surf_gx = np.gradient(surf, axis=1) / grdx
    surf_gy = np.gradient(surf, axis=0) / grdy
    ax, ay = surf_gx[nx_c, ny_c], surf_gy[nx_c, ny_c]
    dx, dy = np.arctan(ax), np.arctan(ay)
    #nx = np.array([np.cos(dx)*np.ones_like(phas), np.zeros_like(phas), np.sin(dx)*np.ones_like(phas)])
    nx = np.array([np.cos(dx), 0, np.sin(dx)])
    nx = v_normal(nx)
    #ny = np.array([np.zeros_like(phas), np.cos(dy)*np.ones_like(phas), np.sin(dy)*np.ones_like(phas)])
    ny = np.array([0, np.cos(dy), np.sin(dy)])
    ny = v_normal(ny)
    return v_outer(nx, ny)
