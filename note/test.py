#%%
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


%ls 
%pwd 
%timeit
x = np.linspace(0, 50, 100)
plt.plot(x, np.sin(x))
plt.show()

print ("hello")
%time

#%%

%ls 
x = np.linspace(0, 100, 100)
plt.plot(x, np.sin(x))
plt.show()

print ("hello")
%time
