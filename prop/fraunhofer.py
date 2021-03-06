"""


fraunhofer.py: calculates 2D Fraunhofer diffraction  (via Fourier Transform)


"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2016"


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import math
from scipy.special import jv
from timeit import default_timer as timer

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d
from src.fraunhofer import propagator2d


def wavefront_initialize(pixelsize_h=1e-6, pixelsize_v=1e-6, npixels_h=1024, npixels_v=1024, amplitude_value=0.0):
    #
    # wavefront definitions
    #
    # create array at object (aperture) plane
    #
    amplitude = np.zeros((npixels_h, npixels_v))  # amplitude map
    amplitude += amplitude_value
    p_i_h = np.arange(npixels_h) * pixelsize_h
    p_x = (p_i_h - 0.5 * (p_i_h[-1] - p_i_h[0]))
    p_i_v = np.arange(npixels_v) * pixelsize_v
    p_y = (p_i_v - 0.5 * (p_i_v[-1] - p_i_v[0]))
    return p_x, p_y, amplitude


def wavefront_aperture(p_x, p_y, amplitude, diameter=40e-6, type=0):
    # aperture_type: 0=circular, 1=Square, 2=Gaussian
    p_xx = p_x[:, np.newaxis]
    p_yy = p_y[np.newaxis, :]

    filter = np.zeros_like(amplitude)
    if type == 0:  # Circular aperture
        radius = (diameter / 2)
        print("radius=%f um" % (1e6 * radius))
        filter_illuminated_indices = np.where(p_xx**2 + p_yy**2 < radius**2)
        if filter_illuminated_indices[0].size == 0:
            print("Warning: wavefront_aperture(): Nothing goes trough the aperture")
        else:
            filter[filter_illuminated_indices] = 1.0
    elif type == 1:  # square
        radius = (diameter / 2)
        print("radius=%f um" % (1e6 * radius))
        filter_illuminated_indices = np.where(
            (np.abs(p_xx) < radius) & (np.abs(p_yy) < radius))
        if filter_illuminated_indices[0].size == 0:
            print("Warning: wavefront_aperture(): Nothing goes trough the aperture")
        else:
            filter[filter_illuminated_indices] = 1.0
    elif type == 2:  # Gaussian
        sigma = diameter / 2.35
        print("source sigma=%f um" % (1e6 * sigma))
        rho2 = p_xx**2 + p_yy**2
        # TODO: add Gaussian amplitude
        # Gaussian in intensity, so srrt for amplitude
        filter = np.sqrt(np.exp(-rho2 / 2 / sigma**2))
        filter = np.exp(-rho2 / 2 / sigma**2)  # Gaussian amplitude
    else:
        raise ValueError("Aperture type (shape) not valid")

    return p_x, p_y, amplitude * filter


def line_image(image, horizontal_or_vertical='H'):
    if horizontal_or_vertical == "H":
        npixels = image.shape[0]
        tmp = image[:, int(image.shape[1] / 2)]
    else:
        npixels = image.shape[1]
        tmp = image[int(image.shape[0] / 2), :]
    return tmp


def line_fwhm(line):
    #
    # CALCULATE fwhm in number of abscissas bins (supposed on a regular grid)
    #
    tt = np.where(line >= max(line) * 0.5)
    if line[tt].size > 1:
        # binSize = x[1]-x[0]
        FWHM = (tt[0][-1] - tt[0][0])
        return FWHM
    else:
        return -1


def plot_image(mymode, theta, psi, title="TITLE", xtitle=r"X [$\mu m$]", ytitle=r"Y [$\mu m$]", name="name.png"):
    fig = plt.figure()
    plt.imshow(mymode.T, origin='lower', extent=[
               theta[0], theta[-1], psi[0], psi[-1]], cmap="jet")
    plt.colorbar()
    ax = fig.gca()
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    plt.title(title)
    plt.savefig(name)


def plot(*positional_parameters, title="", xtitle="", ytitle="", name="name.png", legend=None, color=None):

    n_arguments = len(positional_parameters)
    if n_arguments == 0:
        return

    fig = plt.figure()

    if n_arguments == 1:
        y = positional_parameters[0]
        x = np.arange(y.size)
        plt.plot(x, y, label=legend)
    elif n_arguments == 2:
        x = positional_parameters[0]
        y = positional_parameters[1]
        plt.plot(x, y, label=legend, color=color)
    elif n_arguments == 4:
        x1 = positional_parameters[0]
        y1 = positional_parameters[1]
        x2 = positional_parameters[2]
        y2 = positional_parameters[3]
        if legend != None:
            legend1 = legend[0]
            legend2 = legend[1]
        else:
            legend1 = None
            legend2 = None
        if color != None:
            color1 = color[0]
            color2 = color[1]
        else:
            color1 = None
            color2 = None
        plt.plot(x1, y1, label=legend1, color=color1)
        plt.plot(x2, y2, label=legend2, color=color2)
    else:
        "Incorrect number of arguments, plotting only two first arguments"
        x = positional_parameters[0]
        y = positional_parameters[1]
        plt.plot(x, y, label=legend)

    if legend != None:
        ax = plt.subplot(111)
        ax.legend(bbox_to_anchor=(1.1, 1.05))

    plt.title(title)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.savefig(name)


if __name__ == "__main__":
    #
    # inputs (in SI)
    #

    wavelength = 1.24e-10

    aperture_diameter = 40e-6
    # if Gaussian, aperture_diameter = 2.35*sigma
    # 0=circular, 1=Square, 2=Gaussian (sigma = diameter/2.35)
    aperture_type = 2

    pixelsize_x = 1e-6
    pixelsize_y = 0.5e-6
    npixels_x = 1024
    npixels_y = 2048

    propagation_distance = 3.0

    method = "fourier_convolution"
    #method = "fraunhofer"
    # method = "integral"
    # method = "srw"

    #
    # calculations
    #

    # get a wavefront
    p_x, p_y, amplitude = wavefront_initialize(
        pixelsize_x, pixelsize_y, npixels_x, npixels_y, amplitude_value=1.0)
    # set aperture
    p_x, p_y, amplitude = wavefront_aperture(
        p_x, p_y, amplitude, diameter=aperture_diameter, type=aperture_type)
    # plot aperture
    plot_image(np.abs(amplitude)**2, p_x * 1e6, p_y * 1e6,
               title="aperture intensity, Diameter=%5.1f um" % (
                   1e6 * aperture_diameter),
               xtitle="X [um]", ytitle="Y [um]", name="fraunhofer_aper.png")

    #
    # propagation
    #
    angle_x, angle_y, amplitude_propagated = propagator2d(p_x, p_y, amplitude, method=method, wavelength=wavelength,
                                                          propagation_distance=propagation_distance, return_angles=1)

    # angle_x, angle_y, amplitude_propagated = propagator2d_fraunhoffer(p_x,p_y,amplitude,wavelength=wavelength)
    if method == "fraunhofer":
        print("Fraunhoffer diffraction valid for distances > > a^2/lambda = %f m" %
              ((aperture_diameter / 2)**2 / wavelength))

    plot_image(np.abs(amplitude_propagated)**2, angle_x * 1e6, angle_y * 1e6,
               title="Diffracted intensity (%s)" % method,
               xtitle="X [urad]", ytitle="Y [urad]", name="fraunhofer_prop.png")

    #
    # extract profiles and calculate theoretical ones
    #

    # retrieve H and V profiles
    horizontal_intensity_profile = line_image(
        np.abs(amplitude_propagated)**2, horizontal_or_vertical='H')
    horizontal_intensity_profile /= horizontal_intensity_profile.max()

    vertical_intensity_profile = line_image(
        np.abs(amplitude_propagated)**2, horizontal_or_vertical='V')
    vertical_intensity_profile /= vertical_intensity_profile.max()

    # theoretical profile
    if aperture_type == 0:
        # circular, also display analytical values
        x = (2 * np.pi / wavelength) * (aperture_diameter / 2) * angle_x
        y = (2 * np.pi / wavelength) * (aperture_diameter / 2) * angle_y
        U_vs_theta_x = 2 * jv(1, x) / x
        U_vs_theta_y = 2 * jv(1, y) / y
        I_vs_theta_x = U_vs_theta_x**2
        I_vs_theta_y = U_vs_theta_y**2
    elif aperture_type == 1:
        # square
        x = (2 * np.pi / wavelength) * (aperture_diameter / 2) * angle_x
        y = (2 * np.pi / wavelength) * (aperture_diameter / 2) * angle_y
        U_vs_theta_x = 2 * np.sin(x) / x
        U_vs_theta_y = 2 * np.sin(y) / y
        I_vs_theta_x = U_vs_theta_x**2
        I_vs_theta_y = U_vs_theta_y**2
        I_vs_theta_x /= I_vs_theta_x.max()
        I_vs_theta_y /= I_vs_theta_y.max()
    elif aperture_type == 2:
        # Gaussian
        sigma = aperture_diameter / 2.35
        sigma_ft = 1.0 / sigma * wavelength / (2.0 * np.pi)
        # Factor 2.0 is because we wwant intensity (amplitude**2)
        I_vs_theta_x = np.exp(-2.0 * (angle_x**2 / sigma_ft**2 / 2))
        I_vs_theta_y = np.exp(-2.0 * (angle_y**2 / sigma_ft**2 / 2))
    else:
        # Gaussian
        sigma = aperture_diameter / 2.35
        sigma_ft = 1.0 / sigma * wavelength / (2.0 * np.pi)
        # Factor 2.0 is because we wwant intensity (amplitude**2)
        I_vs_theta_x = np.exp(-2.0 * (angle_x**2 / sigma_ft**2 / 2))
        I_vs_theta_y = np.exp(-2.0 * (angle_y**2 / sigma_ft**2 / 2))

    fwhm_intensity_profile_horizontal = line_fwhm(
        horizontal_intensity_profile) * (angle_x[1] - angle_x[0])
    fwhm_intensity_profile_vertical = line_fwhm(
        vertical_intensity_profile) * (angle_y[1] - angle_y[0])
    fwhm_theoretical_profile_horizontal = line_fwhm(
        I_vs_theta_x) * (angle_x[1] - angle_x[0])
    fwhm_theoretical_profile_vertical = line_fwhm(
        I_vs_theta_y) * (angle_y[1] - angle_y[0])

    #
    # calculate widths
    #
    print("HORIZONTAL FWHM (%s) : %f urad, FWHM theoretical: %f urad, 1.22*wavelength/Diameter: %f urad" % (
        method,
        1e6 * fwhm_intensity_profile_horizontal, 1e6 * fwhm_theoretical_profile_horizontal, 1e6 * 1.22 * wavelength / aperture_diameter))
    print("VERTICAL FWHM (%s) : %f urad, FWHM theoretical: %f urad, 1.22*wavelength/Diameter: %f urad" % (
        method,
        1e6 * fwhm_intensity_profile_vertical, 1e6 * fwhm_theoretical_profile_vertical, 1e6 * 1.22 * wavelength / aperture_diameter))
    print("HORIZONTAL (4pi/lambda) sigma sigma' : (%s): %f, theoretical: %f " % (method,
                                                                                 4 * np.pi / wavelength * fwhm_intensity_profile_horizontal /
                                                                                 2.35 * aperture_diameter / 2.35,
                                                                                 4 * np.pi / wavelength * fwhm_theoretical_profile_horizontal / 2.35 * aperture_diameter / 2.35))

    # plot profiles
    plot(angle_x * 1e6, horizontal_intensity_profile, angle_x * 1e6, I_vs_theta_x,
         legend=["profile", "theory"], color=["red", "black"],
         title="Horizontal profile of diffracted intensity (%s)" % method, xtitle='theta [urad]', ytitle='Diffracted intensity [a.u.]', name="./fraunhofer_prof1.png")
    plot(angle_y * 1e6, vertical_intensity_profile, angle_y * 1e6, I_vs_theta_y,
         legend=["profile", "theory"], color=["red", "black"],
         title="Vertical profile of diffracted intensity (%s)" % method, xtitle='theta [urad]', ytitle='Diffracted intensity [a.u.]', name="./fraunhofer_prof2.png")
