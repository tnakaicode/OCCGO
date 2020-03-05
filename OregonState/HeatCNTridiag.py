""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

# HeatCNTridiag.py:  solution of heat eqtn via CN method

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


Max = 51
n = 50
m = 50
Ta = np.zeros((Max), float)
Tb = np.zeros((Max), float)
Tc = np.zeros((Max), float)
Td = np.zeros((Max), float)
a = np.zeros((Max), float)
b = np.zeros((Max), float)
c = np.zeros((Max), float)
d = np.zeros((Max), float)
x = np.zeros((Max), float)
t = np.zeros((Max, Max), float)


def Tridiag(a, d, c, b, Ta, Td, Tc, Tb, x, n):
    Max = 51
    h = np.zeros((Max), float)
    p = np.zeros((Max), float)
    for i in range(1, n + 1):
        a[i] = Ta[i]
        b[i] = Tb[i]
        c[i] = Tc[i]
        d[i] = Td[i]
    h[1] = c[1] / d[1]
    p[1] = b[1] / d[1]
    for i in range(2, n + 1):
        h[i] = c[i] / (d[i] - a[i] * h[i - 1])
        p[i] = (b[i] - a[i] * p[i - 1]) / (d[i] - a[i] * h[i - 1])
    x[n] = p[n]
    for i in range(n - 1, 1, -1):
        x[i] = p[i] - h[i] * x[i + 1]


width = 1.0
height = 0.1
ct = 1.0
for i in range(0, n):
    t[i, 0] = 0.0
for i in range(1, m):
    t[0][i] = 0.0
h = width / (n - 1)
k = height / (m - 1)
r = ct * ct * k / (h * h)

for j in range(1, m + 1):
    t[1, j] = 0.0
    t[n, j] = 0.0
    # BCs
for i in range(2, n):
    t[i][1] = np.sin(np.pi * h * i)
    # ICs
for i in range(1, n + 1):
    Td[i] = 2. + 2. / r
Td[1] = 1.
Td[n] = 1.
for i in range(1, n):
    Ta[i] = -1.0
    Tc[i] = -1.0
    # Off diagonal
Ta[n - 1] = 0.0
Tc[1] = 0.0
Tb[1] = 0.0
Tb[n] = 0.0
print("I'm working hard, wait for fig while I count to 50")

for j in range(2, m + 1):
    print(j)
    for i in range(2, n):
        Tb[i] = t[i - 1][j - 1] + t[i + 1][j - 1] \
            + (2 / r - 2) * t[i][j - 1]
    Tridiag(a, d, c, b, Ta, Td, Tc, Tb, x, n)
    for i in range(1, n + 1):
        t[i][j] = x[i]
print("Finished")
x = list(range(1, m + 1))
y = list(range(1, n + 1))
X, Y = np.meshgrid(x, y)


def functz(t):
    # Potential
    z = t[X, Y]
    return z


Z = functz(t)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X, Y, Z, color='r')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('T')
plt.show()
