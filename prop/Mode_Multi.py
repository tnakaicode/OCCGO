import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import scipy.fftpack
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def ft(a):
    A = scipy.fftpack.fftshift(scipy.fftpack.fft2(scipy.fftpack.fftshift(a)))
    r = np.sum(np.abs(a)**2) / np.sum(np.abs(A)**2)
    return r**A


def ift(a):
    A = scipy.fftpack.ifftshift(
        scipy.fftpack.ifft2(scipy.fftpack.ifftshift(a)))
    r = np.sum(np.abs(a)**2) / np.sum(np.abs(A)**2)
    return r * A


def ls(H, F):
    M, N, N = H.shape
    A = np.zeros((N, N), dtype="complex128")
    B = np.zeros((N, N), dtype="complex128")
    G = np.zeros((N, N), dtype="complex128")
    for m in range(0, M):
        A += H[m] .conj() * F[m]
        B += np.abs(H[m])**2
    G[B != 0] = A[B != 0] / B[B != 0]
    return G


def fibre_image(u1, NZ):
    for i in range(0, NZ):
        u1 = u1 * np.exp(1j * 2 * np.pi * fibre_profile[i] * LZ / lmbda / NZ)
        u1 = prop(u1, LZ / NZ)
    u1 = prop(u1, LZ / NZ)
    return u1


def prop(u1, z):
    if dx >= (lmbda * z / L):
        fx = np.arange(-1 / (2 * dx), 1 / (2 * dx), 1.0)
        fy = np.arange(-1 / (2 * dx), 1 / (2 * dx), 1.0)
        FX, FY = np.meshgrid(fy, fx)
        H = np.exp(-1j * np.pi * lmbda * z * (FX**2 + FY**2))
        U1 = ft(u1)
        U2 = H * U1
        u2 = ift(U2)
    else:
        h = 1 / (1j * lmbda * z) * np.exp(1j * k / (2 * z) * (PX**2 + PY**2))
        H = ft(h)
        U1 = ft(u1)
        U2 = H * U1
        u2 = ift(U2)
    return u2


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    # Identification and control of light propagation in optical waveguide

    ## -- Code from here -- ##
    N = 512 * 2  # Grid size
    maxit = 20  # Number of iterations
    D = 20  # Number of images

    # fibre parameters
    NZ = 2  # Propagation length
    lmbda = 1.0  # Wavelength(normalised)
    L = 1.0  # Support length(normalised)
    LZ = 0.0002  # Propagation step length
    dx = L / N  # Sample size
    k = 2 * np.pi / lmbda  # Wave number

    ## -- fibre specifications -- ##
    # Fibreprofile specifications
    px = np.linspace(-1.0, 1.0, N)
    pz = np.linspace(-1.0, 1.0, NZ)
    PX, PY = np.meshgrid(px, px)
    RHO = np.sqrt(PX**2 + PY**2)
    THETA = np.arctan2(PY, PX)

    # Define the fibre
    inner_diameter = 0.50
    outer_diameter = 0.70
    phase_error = 10000.0
    pp = (np.ones((N, N)) * (RHO < inner_diameter)).astype(float)

    # Refractive index of the fibre
    fibre_profile = np.zeros((NZ, N, N), dtype="float64")
    fibre_profile[:, RHO > outer_diameter] = 1.00
    fibre_profile[:, RHO < outer_diameter] = 1.40
    fibre_profile[:, RHO < inner_diameter] = 1.45

    # Imperfections in the fibre
    for i in range(0, NZ):
        error = np.abs(ift(ft(np.random.randn(N, N)) *
                           np.exp(-(PX**2 + PY**2) / 1.0**2)))
        error = error / error.max()
        fibre_profile[i] = (fibre_profile[i] + phase_error * error)

    ## -- Multiframe Deconvolution -- ##

    # Initialisation
    y = np.zeros((D, N, N), dtype=complex)
    x = np.zeros((D, N, N), dtype=complex)
    X = np.zeros((D, N, N), dtype=complex)
    y_hat = np.zeros((D, N, N), dtype=complex)
    Y_HAT = np.zeros((D, N, N), dtype=complex)
    O_HAT = np.random.randn(N, N) + 1j * np.random.randn(N, N)

    # Acquire D input-output image pairs
    for d in range(0, D):
        print(d)
        x[d] = pp * np.exp(1j * np.random.randn(N, N))
        X[d] = ft(x[d])
        y[d] = fibre_image(x[d], NZ)

    # Find O via deconvolution
    for i in range(0, maxit):
        for d in range(0, D):
            # magnitude constraint
            y_hat[d] = np.abs(y[d]) * np.exp(1j * np.angle(ift(O_HAT * X[d])))
            Y_HAT[d] = ft(y_hat[d])
            O_HAT = ls(X, Y_HAT)
            O_HAT = np.exp(1j * np.angle(O_HAT))  # phase only constraint
            print(i, d)

    ## -- Imaging -- ##

    # Define target intensity (TUDelft flame)
    # img = pil.open(’TUflame.png’)
    # img = np.array(img)
    # img = img[:,:,0]
    # img = 1.0/(img+0.0001)
    # img[img == np.min(img)] = 0.0
    # img = downsample(downsample(downsample(img)))
    # pad = (N-np.shape(img)[0])/2
    # img = np.pad(img,((pad,pad),(pad,pad)),’edge’)
    # y_target = np.copy(img)
    # Y_TARGET = ft(y_target)

    # Define target intensity (spot)
    y_target = np.zeros((N, N))
    y_target[int(N / 2), int(N / 2)] = 1.0
    Y_TARGET = ft(y_target)

    # Compute input/output (PSF model)
    Y_OUT_HAT = np.copy(Y_TARGET)
    x_target_hat = pp * ift(Y_OUT_HAT / O_HAT)
    X_TARGET_HAT = ft(pp * np.exp(1j * np.angle(x_target_hat)))

    y_out_hat = ift(O_HAT * X_TARGET_HAT)
    y_out_hat /= np.sum(np.abs(y_out_hat)**2)

    # Copute input/output (Fibre simulation)
    X_TARGET = Y_TARGET / O_HAT
    x_target = pp * np.exp(1j * np.angle(ift(X_TARGET)))
    y_out = fibre_image(x_target, NZ)
    y_out /= np.sum(np.abs(y_out)**2)
    Y_OUT = ft(np.abs(y_out))

    ## -- Plot results -- ##
    phi = np.zeros((N, N))
    rad = pp * (np.angle(x_target) + np.pi) / (2.0 * np.pi)
    phi[255, 255:769] = 1.0
    phi[769, 255:769] = 1.0
    phi[255:769, 255] = 1.0
    phi[255:769, 769] = 1.0
    phi += rad
    phi = phi[254:N - 254, 254:N - 254]
    images = [np.abs(y_out_hat[254:N - 254, 254:N - 254])**2,
              np.abs(y_out[254:N - 254, 254:N - 254])**2, phi]
    title = ["Model intensity", "Fibre simulator intensity",
             "SLM phase pattern(rad)"]

    obj = plot2d()
    ax1 = obj.add_axs(1, 3, 1)
    im = ax1.imshow(images[0], cmap="jet")
    a = plt.axes([.2, .3, .15, .15])
    a.imshow(np.abs(y_out_hat[254 + 225:N - 254 - 225,
                              254 + 225:N - 254 - 225])**2, cmap="jet")
    plt.setp(a, xticks=[], yticks=[])
    a = plt.axes([.42, .3, .15, .15])
    a.imshow(np.abs(y_out[254 + 225:N - 254 - 225,
                          254 + 225:N - 254 - 225])**2, cmap="jet")
    plt.setp(a, xticks=[], yticks=[])
    ax1.set_title(title[0])

    ax2 = obj.add_axs(1, 3, 2)
    im = ax2.imshow(images[0], cmap="jet")
    a = plt.axes([.2, .3, .15, .15])
    a.imshow(np.abs(y_out_hat[254 + 225:N - 254 - 225,
                              254 + 225:N - 254 - 225])**2, cmap="jet")
    plt.setp(a, xticks=[], yticks=[])
    a = plt.axes([.42, .3, .15, .15])
    a.imshow(np.abs(y_out[254 + 225:N - 254 - 225,
                          254 + 225:N - 254 - 225])**2, cmap="jet")
    plt.setp(a, xticks=[], yticks=[])
    ax2.set_title(title[1])

    ax3 = obj.add_axs(1, 3, 3)
    im = ax3.imshow(images[0], cmap="jet")
    a = plt.axes([.2, .3, .15, .15])
    a.imshow(np.abs(y_out_hat[254 + 225:N - 254 - 225,
                              254 + 225:N - 254 - 225])**2, cmap="jet")
    plt.setp(a, xticks=[], yticks=[])
    a = plt.axes([.42, .3, .15, .15])
    a.imshow(np.abs(y_out[254 + 225:N - 254 - 225,
                          254 + 225:N - 254 - 225])**2, cmap="jet")
    plt.setp(a, xticks=[], yticks=[])
    ax3.set_title(title[2])

    obj.SavePng()
