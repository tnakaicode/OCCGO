""" From "COMPUTATIONAL PHYSICS" & "COMPUTER PROBLEMS in PHYSICS"
    by RH Landau, MJ Paez, and CC Bordeianu (deceased)
    Copyright R Landau, Oregon State Unv, MJ Paez, Univ Antioquia, 
    C Bordeianu, Univ Bucharest, 2018. 
    Please respect copyright & acknowledge our work."""

# LaplaceLine.py:  Solve Laplace's eqtn within square

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Nmax = 100
Niter = 50
V = np.zeros((Nmax, Nmax), float)
print("Working hard, wait for the figure while I count to 60")

for k in range(25, 75):
    V[0, k] = 100.0      # Line at 100V
for iter in range(Niter):
    if iter % 10 == 0:
        print(iter)
    for i in range(1, Nmax-2):
        for j in range(1, Nmax-2):
            V[i, j] = 0.25*(V[i+1, j]+V[i-1, j]+V[i, j+1]+V[i, j-1])
    print("iter, V[Nmax/5,Nmax/5]", iter, V[Nmax//5, Nmax//5])
x = range(0, 50, 2)
y = range(0, 50, 2)
X, Y = np.meshgrid(x, y)
Z = V[X, Y]

fig = plt.figure()                            # Create figure
ax = Axes3D(fig)                                # Plot axes
ax.plot_wireframe(X, Y, Z, color='r')      # Red wireframe
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('V(x,y)')
ax.set_title('Potential within Square V(x=0)=100V (Rotatable)')
plt.show()                                          # Show fig
