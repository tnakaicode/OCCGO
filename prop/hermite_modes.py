"""


hermite_modes.py: calculates Gauss-Schell modes

see e.g. https://commons.wikimedia.org/wiki/File:Hermite-gaussian.png

"""

__author__ = "M. Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2016"


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite


def plot_image(*positional_parameters, title="TITLE", xtitle=r"X", ytitle=r"Y",
               xrange=None, yrange=None,
               cmap=None, aspect=None, name="name.png",
               add_colorbar=True, figsize=None):

    n_arguments = len(positional_parameters)
    if n_arguments == 1:
        z = positional_parameters[0]
        x = np.arange(0, z.shape[0])
        y = np.arange(0, z.shape[1])
    elif n_arguments == 2:
        z = positional_parameters[0]
        x = positional_parameters[1]
        y = positional_parameters[1]
    elif n_arguments == 3:
        z = positional_parameters[0]
        x = positional_parameters[1]
        y = positional_parameters[2]
    else:
        raise Exception("Bad number of inputs")

    fig = plt.figure(figsize=figsize)

    plt.imshow(z.T, origin='lower', extent=[
               x[0], x[-1], y[0], y[-1]], cmap="jet", aspect="equal")
    if add_colorbar:
        plt.colorbar()
    ax = fig.gca()
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    plt.title(title)

    plt.xlim(xrange)
    plt.ylim(yrange)
    plt.savefig(name)


if __name__ == '__main__':
    size_x = 1.0
    size_y = 1.0
    n_x = 620
    n_y = 320
    w = size_x / 8.
    m = 1
    n = 1

    X = np.linspace(-size_x / 2, size_x / 2, n_x)
    Y = np.linspace(-size_y / 2, size_y / 2, n_y)

    XX = np.outer(X, np.ones_like(Y))
    YY = np.outer(np.ones_like(X), Y)

    out = (hermite(m)(np.sqrt(2) * XX / w) * np.exp(-XX**2 / w**2))**2 \
        * (hermite(n)(np.sqrt(2) * YY / w) * np.exp(-YY**2 / w**2))**2

    plot_image(out, X, Y, name="hermite_modes.png")
