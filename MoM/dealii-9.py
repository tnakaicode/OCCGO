import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import scipy.constants as cnt
sys.path.append(os.path.join('../'))

from src.plot import plot_contour_sub
from src.profile import mask_eclip

from pymor.algorithms import genericsolvers

"""

https://www.dealii.org/current/doxygen/deal.II/step_9.html

GMRES solver
Generalized Minimum RESidual method
https://www.dealii.org/current/doxygen/deal.II/classSolverGMRES.html

solve for u: b * Delta(u) = f
u = g on Boundary Omeg

d    = 2
omeg = [−1, 1] ^ d
b(x) = (2, 1 + 4 / 5 sin(8πx))
s    = 0.1
f(x) = 1 / (10 * s ^ d) for | x - x_0 | < s
f(x) = 0 else
x_o  = (-3 / 4, -3 / 4)
g    = exp(5 (1−| x | ^2) * sin(16pi | x | ^2)

* A simple refinement criterion
* error estimator first developed by Kelly
ita_K = {h_K/24 integral (del K) (del u_h)^2 ds}^(1/2)
ita_K ~ C * h * del^2(u)_K

* Equation data declaration
* AdvectionProblem class declaration
* GradientEstimation class declaration
* AdvectionProblem class implementation
GMRES solver
* GradientEstimation class implementation
* Main function

"""


class Problem (object):

    def __init__(self):
        self.dim = 2
        self.x_0 = [-3 / 4, -3 / 4]
        self.s_0 = 0.1

    def b_x(self, x, y):
        return 2 * np.ones_like(x), 1 + 4 / 5 * np.sin(8 * np.pi * x)

    def f_x(self, x, y):
        px = self.x_0[0] - x
        py = self.x_0[1] - y
        mask = mask_eclip([x, y], self.x_0, [self.s_0, self.s_0])
        func = 1 / (10 * self.s_0**self.dim)
        func *= mask
        return func

    def g_x(self, x, y):
        abs_x = np.sqrt(x**2 + y**2)
        val = np.exp(5 * (1 - abs_x**2)) * np.sin(16 * np.pi * abs_x**2)
        return val


if __name__ == '__main__':
    prob = Problem()
    print(prob.x_0)
    print(prob.b_x(0, 0))

    px = np.linspace(-1, 1, 100)
    py = np.linspace(-1, 1, 100)
    mesh = np.meshgrid(px, py)

    func = prob.f_x(*mesh)

    plot_contour_sub(mesh, prob.b_x(*mesh)[0], dirname="bx_0")
    plot_contour_sub(mesh, prob.b_x(*mesh)[1], dirname="bx_1")
    plot_contour_sub(mesh, prob.f_x(*mesh), dirname="func")
    plot_contour_sub(mesh, prob.g_x(*mesh), dirname="g_func")
