from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

a13 = 2.5
a21 = 1.5
a31 = 0.2
r2 = 0.6
r3 = 4.5
k3 = 1
d3 = 0.5

def cancer_model(x, y, z, a12):
    x_dot = x - x*x - a12*x*y - a13*x*z
    y_dot = r2*y - r2*y*y - a21*x*y
    z_dot = (r3*x*z)/(x+k3) - a31*x*z - d3*z
    return x_dot, y_dot, z_dot


dr = 0.001  # parameter step size
a12 = np.arange(0, 3, dr)  # parameter range
dt = 0.001  # time step
t = np.arange(0, 50, dt)  # time range

# Solution arrays
xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
zs = np.empty(len(t) + 1)

# Initial values
xs[0], ys[0], zs[0] = (0.1, 0.1, 0.1)

z_maxes = []
r_mins = []
z_mins = []

for R in a12:
    for i in range(len(t)):
        # Approximate solutions
        x_dot, y_dot, z_dot = cancer_model(xs[i], ys[i], zs[i], R)
        xs[i + 1] = xs[i] + (x_dot * dt)
        ys[i + 1] = ys[i] + (y_dot * dt)
        zs[i + 1] = zs[i] + (z_dot * dt)
    # Compute z
    for i in range(1, len(zs) - 1):
        # Max
        if zs[i - 1] < zs[i] and zs[i] > zs[i + 1]:
            r_maxes.append(R)
            z_maxes.append(zs[i])
        # Min
        elif zs[i - 1] > zs[i] and zs[i] < zs[i + 1]:
            r_mins.append(R)
            z_mins.append(zs[i])

    # Use final values as next innitial values
    xs[0], ys[0], zs[0] = xs[i], ys[i], zs[i]
# Plot

plt.scatter(r_maxes, z_maxes, color="black", s=0.5, alpha=0.2)
plt.scatter(r_mins, z_mins, color="red", s=0.5, alpha=0.2)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.show()
