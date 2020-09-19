"""

Field Computation by Moment Methods - Harrington
Sec 1.4 Example

- d^2 f / d x^2 = 1 + x^2
f(0) = f(1) = 0

Solution:

f(x) = 5x/6 - x^2/2 - x^4/3

Power Series Solution

L = - d^2 / dx^2

f_n(x) = x - x^(n+1)

- d^2 f / d x^2 = 1 + x^2 = sum (n=1;N) (a_n * L(f_n(x)))


l_mn = n * (n + 1) * (m / (N + 1))^(n-1)
g_m  = 1 + 4 * (m / (N + 1))^2
w_m  = delta (x - x_m)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

sys.path.append(os.path.join('../'))
from src.base import plot2d


def gen_coef(dim=1):
    vec = np.empty((dim))
    for i in range(dim):
        m = i + 1
        val = 1 + 4 * (m / (dim + 1))**2
        vec[i] = val
    return vec


def gen_mat(dim=1):
    mat = np.empty((dim, dim))
    for (i0, i1), val in np.ndenumerate(mat):
        m, n = i0 + 1, i1 + 1
        val = n * (n + 1) * (m / (dim + 1))**(n - 1)
        mat[i0, i1] = val
    return mat


def gen_fmn(px, dim=1):
    mat = gen_mat(dim)
    vec = gen_coef(dim)
    coef = np.dot(vec, np.linalg.inv(mat))

    print("dim", dim)
    print(mat)
    print(vec)
    print(coef)

    py = np.empty_like(px)
    for i in range(dim):
        n = i + 1
        py += coef[i] * (px - px**(n))
    return py


if __name__ == '__main__':
    px = np.linspace(0, 1, 100)

    f1 = gen_fmn(px, 1)
    f2 = gen_fmn(px, 2)
    f3 = gen_fmn(px, 3)
    f4 = gen_fmn(px, 4)
    f5 = gen_fmn(px, 5)
    f6 = gen_fmn(px, 6)

    f_ex = - px**4 / 3 - px**2 / 2 + 5 * px / 6

    obj = plot2d("auto")
    obj.axs.plot(px, f_ex)
    #obj.axs.plot(px, f1)
    #obj.axs.plot(px, f2)
    obj.axs.plot(px, f3)
    obj.axs.plot(px, f4)
    obj.axs.plot(px, f5)
    obj.axs.plot(px, f6)

    obj.axs.plot(px, gen_fmn(px, 4))
    obj.axs.plot(px, gen_fmn(px, 5))
    obj.SavePng()
