import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import scipy as sci
import scipy.optimize
import scipy.special
import scipy.sparse
import scipy.signal
import argparse
from linecache import getline, clearcache

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def pol2cart(rho, theta):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return rho, theta


def findbopt(V, l):
    dlta = 0.00001
    b = np.arange(0, 1, dlta)
    zeros = scipy.special.jn_zeros(l, 100)
    zeros = zeros[zeros < V]
    bzeros = 1 - (zeros / V)**2
    right = V * np.sqrt(b) * sci.special.kv(l + 1, V *
                                            np.sqrt(b)) / sci.special.kv(l, V * np.sqrt(b))
    left = V * np.sqrt(1 - b) * sci.special.jv(l + 1, V *
                                               np.sqrt(1 - b)) / sci.special.jv(l, V * np.sqrt(1 - b))
    d = left - right
    s = 0
    bopt = np.zeros(100)
    L = np.zeros(100) - 1
    for i in range(0, len(d) - 1):
        if(d[i] == 0 or d[i] * d[i + 1] < 0) and not np.isclose(bzeros, i * dlta, atol=0.001).any():
            bopt[s] = i * dlta
            L[s] = l
            s += 1
    bopt = bopt[bopt > 0]
    bopt = bopt[::-1]
    L = L[L > -0.5]
    return bopt, L


def singleLP(V, l, m, mode=1, a=1, clad=0.1, npix=512):
    x = np.linspace(-(a + clad), (a + clad), npix)
    X, Y = np.meshgrid(x, x)
    RHO_core, PHI = cart2pol(X, Y)
    RHO_core[RHO_core > a] = 0.0
    RHO_clad, PHI = cart2pol(X, Y)
    RHO_clad[RHO_clad <= a] = 0.0
    bopt, L = findbopt(V, l)
    if m > len(bopt):
        print("This mode does not exist")
        Psicos = np.zeros(PHI.shape)
        Psisin = np.zeros(PHI.shape)
    else:
        bopt = bopt[m - 1]
        U = V * np.sqrt(1 - bopt)
        W = V * np.sqrt(bopt)
        R1 = sci.special.jv(l, U * RHO_core / a) / sci.special.jv(l, U)
        R2 = sci.special.kv(l, W * RHO_clad / a) / sci.special.kv(l, W)
        R2[np.isnan(R2)] = 0.0
        R = R1 + R2
        Psicos = R * np.cos(l * PHI)
        Psisin = R * np.sin(l * PHI)

    if mode == 1:
        Psi = Psicos
    else:
        Psi = Psisin
    return Psi


def list1(V, lmax):
    bopt = 0
    L = 0
    for l in range(0, lmax):
        boptnew, Lnew = findbopt(V, l)
        bopt = np.append(bopt, boptnew)
        L = np.append(L, Lnew)
    bopt = bopt[1::]
    L = L[1::]
    return bopt, L


def even(bopt, L, mode, V, a=1, clad=0.1, npix=512):
    x = np.linspace(-(a + clad), (a + clad), npix)
    X, Y = np.meshgrid(x, x)
    RHO_core, PHI = cart2pol(X, Y)
    RHO_core[RHO_core >= a] = 0.0
    RHO_clad, PHI = cart2pol(X, Y)
    RHO_clad[RHO_clad < a] = 0.0
    U = V * np.sqrt(1 - bopt[mode])
    W = V * np.sqrt(bopt[mode])
    R1 = sci.special.jv(L[mode], U * RHO_core / a) / sci.special.jv(L[mode], U)
    R2 = sci.special.kv(L[mode], W * RHO_clad / a) / sci.special.kv(L[mode], W)
    R2[np.isnan(R2)] = 0.0
    R = R1 + R2
    Psicos = R * np.cos(L[mode] * PHI)
    Psicos /= np.linalg.norm(Psicos, "fro")
    return Psicos


def odd(bopt, L, mode, V, a=1, clad=0.1, npix=512):
    x = np.linspace(-(a + clad), (a + clad), npix)
    X, Y = np.meshgrid(x, x)
    RHO_core, PHI = cart2pol(X, Y)
    RHO_core[RHO_core >= a] = 0.0
    RHO_clad, PHI = cart2pol(X, Y)
    RHO_clad[RHO_clad < a] = 0.0
    U = V * np.sqrt(1 - bopt[mode])
    W = V * np.sqrt(bopt[mode])
    R1 = sci.special.jv(L[mode], U * RHO_core / a) / sci.special.jv(L[mode], U)
    R2 = sci.special.kv(L[mode], W * RHO_clad / a) / sci.special.kv(L[mode], W)
    R2[np.isnan(R2)] = 0.0
    R = R1 + R2
    Psisin = R * np.sin(L[mode] * PHI)
    if Psisin.all() == 0:
        Psisin = np.zeros(Psisin.shape)
    else:
        Psisin /= np.linalg.norm(Psisin, "fro")
    return Psisin


def field2(bopt, L, coeff, V, a=1, clad=0.1, npix=512):
    Psi = np.zeros((npix, npix), dtype=complex)
    for i in range(0, len(bopt)):
        Psicos = even(bopt, L, i, V, a, clad, npix)
        Psisin = odd(bopt, L, i, V, a, clad, npix)
        Psi += coeff[2 * i] * Psicos + coeff[2 * i + 1] * Psisin
    return Psi


def MatPur(Psi, bopt, L, V, N):
    coeff_hat = np.zeros(2 * len(bopt), dtype=complex)
    residual = Psi
    i = 0
    idx = []
    while np.sum(np.abs(residual)) > 0.1 and i <= N:
        print(i)
        co = 0.0
        for n in range(0, 2 * len(bopt)):
            if n not in idx:
                if n % 2 == 0:
                    mode = even(bopt, L, n / 2, V)
                else:
                    mode = odd(bopt, L, n / 2, V)
                inprod = np.sum(residual * mode.conj())
                if np.abs(inprod) > np.abs(co):
                    co = inprod
                    ci = n
                    atom = mode
        coeff_hat[ci] = co
        idx.append(ci)
        residual -= co * atom
        i += 1
    return coeff_hat, residual


def fastMatPur(Psi, bopt, L, V, N):
    coeff_hat = np.zeros(2 * len(bopt), dtype=complex)
    residual = Psi
    i = 0
    idx = []
    while np.sum(np.abs(residual)) > 1 and i < N:
        co = 0.0
        corr = np.zeros(2 * len(bopt), dtype=complex)
        for n in range(0, 2 * len(bopt)):
            if n not in idx:
                if n % 2 == 0:
                    mode = even(bopt, L, n / 2, V)
                else:
                    mode = odd(bopt, L, n / 2, V)
                inprod = np.sum(residual * mode.conj())
                corr[n] = inprod
                if np.abs(inprod) > np.abs(co):
                    co = inprod
                    ci = n
                    atom = mode
        coeff_hat[ci] = co
        idx.append(ci)
        residual -= co * atom
        i += 1
        print(i)
        # Continue with a subset of atoms with high correlation
        corr[np.abs(corr / max(corr)) < 0.9] = 0.0
        nn = np.nonzero(corr)
        for ii in np.nditer(nn):
            if ii not in idx:
                if ii % 2 == 0:
                    mode = even(bopt, L, ii / 2, V)
                else:
                    mode = odd(bopt, L, ii / 2, V)
                inprod = np.sum(residual * mode.conj())
                if np.abs(inprod) > np.abs(co):
                    co = inprod
                    ci = ii
                    atom = mode
        coeff_hat[ci] = co
        idx.append(ci)
        residual -= co * atom
        i += 1
        print(i)
    return coeff_hat, residual


def inprod(Psi, bopt, L, V, N):
    coeff_hat = np.zeros(2 * len(bopt), dtype=complex)
    residual = Psi
    for n in range(0, 2 * len(bopt)):
        print(n)
        if n % 2 == 0:
            mode = even(bopt, L, int(n / 2), V)
        else:
            mode = odd(bopt, L, int(n / 2), V)
        coeff_hat[n] = np.sum(residual * mode.conj())
        residual -= coeff_hat[n] * mode
    return coeff_hat, residual


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    # Identification and control of light propagation in optical waveguide

    V = 8  # V-number
    lmax = 5  # maximum number of l modes
    k = 0.8  # amount of zeros 0<k<1

    # compute b-values
    bopt, L = list1(V, lmax)

    # generate random modal coefficients
    coeff = np.random.normal(size=2 * len(bopt)) + \
        1j * np.random.normal(size=2 * len(bopt))
    coeff /= np.linalg.norm(coeff, 2)

    # compute corresponding field
    N = np.size(np.nonzero(coeff))
    F = field2(bopt, L, coeff, V)

    ## -- Decompose field into modes -- ##
    t1 = time.time()
    coeff_hat, residual = inprod(F, bopt, L, V, N)
    t = time.time() - t1
    print("time elapsed:", t)

    # compute corresponding field
    Psi_hat = field2(bopt, L, coeff_hat, V)
    F = field2(bopt, L, coeff, V)

    obj = plot2d(aspect="equal")
    ax1 = obj.add_axs(1, 3, 1)
    ax1.contourf(np.abs(F))
    ax2 = obj.add_axs(1, 3, 2)
    ax2.contourf(np.abs(Psi_hat))
    ax3 = obj.add_axs(1, 3, 3)
    ax3.contourf(np.abs(residual))
    obj.SavePng()

    obj.new_2Dfig()
    obj.axs.set_axis_off()
    obj.axs = obj.fig.add_subplot(111, projection='polar')
    for c in coeff:
        plt.polar(np.angle(c), np.abs(c), marker='o')
    for c in coeff_hat:
        plt.polar(np.angle(c), np.abs(c), marker='*')
    obj.SavePng(obj.tempname + "_polar.png")
