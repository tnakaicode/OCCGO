"""

fresnel: 

        functions: 
             goFromTo: calculates the phase shift matrix
 

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2012"

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import math

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d


def goFromTo(source, image, distance=1.0, lensF=None, wavelength=1e-10):
    distance = np.array(distance)
    x1 = np.outer(source, np.ones(image.size))
    x2 = np.outer(np.ones(source.size), image)
    r = np.sqrt(np.power(x1 - x2, 2) + np.power(distance, 2))
    # add lens at the image plane
    if lensF != None:
        r = r - np.power(x1 - x2, 2) / lensF
    wavenumber = np.pi * 2 / wavelength
    return np.exp(1.j * wavenumber * r)


if __name__ == '__main__':

    # wavelength   =   1e-10
    # aperture_diameter   =   10e-6
    # detector_size = 0.8e-3
    # #wavelength   =   500e-9
    # #aperture_diameter   =   1e-3
    # #detector_size = 4e-3
    #
    # sourcepoints = 1000
    # detpoints =  1000
    # distance =   1.00
    # lensF        =   None

    # wavelength   =   5000e-10
    # sourcesize   =   500e-6
    # detector_size = 0.008
    #wavelength   =   500e-9
    #aperture_diameter   =   1e-3
    #detector_size = 4e-3

    wavelength = 1.24e-10  # 10keV
    aperture_diameter = 1e-6 # 40e-6  # 1e-3 # 1e-6
    detector_size = 800e-6
    distance = 4.0

    sourcepoints = 2000
    detpoints = 2000
    lensF = None

    sourcesize = aperture_diameter

    position1x = np.linspace(-sourcesize / 2, sourcesize / 2, sourcepoints)
    position2x = np.linspace(-detector_size / 2, detector_size / 2, detpoints)

    fields12 = goFromTo(position1x, position2x, distance,
                        lensF=lensF, wavelength=wavelength)
    print("Shape of fields12: ", fields12.shape)

    # prepare results
    fieldComplexAmplitude = np.dot(np.ones(sourcepoints), fields12)
    print("Shape of Complex U: ", fieldComplexAmplitude.shape)
    print("Shape of position1x: ", position1x.shape)
    fieldIntensity = np.power(np.abs(fieldComplexAmplitude), 2)
    fieldPhase = np.arctan2(np.real(fieldComplexAmplitude),
                            np.imag(fieldComplexAmplitude))

    #
    # write spec formatted file
    #
    out_file = "fresnel_kirchhoff_1D.spec"
    f = open(out_file, 'w')
    header = "#F %s \n\n#S  1 fresnel-kirchhoff diffraction integral\n#N 3 \n#L X[m]  intensity  phase\n" % out_file

    f.write(header)

    for i in range(detpoints):
        out = np.array((position2x[i], fieldIntensity[i], fieldPhase[i]))
        f.write(("%20.11e " * out.size + "\n") % tuple(out.tolist()))

    f.close()
    print("File written to disk: %s" % out_file)

    obj = plot2d("auto")
    obj.axs.plot(position2x * 1e6, fieldIntensity)
    obj.axs.set_title("Fresnel-Kirchhoff Diffraction")
    obj.axs.set_xlabel("X [um]")
    obj.axs.set_ylabel("Intensity [a.u.]")
    obj.SavePng()
