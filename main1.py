from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

###### Model
Dt = 0.001 # Size of time step
a12 = 1.05 # Chaotic behaviour if a12 = 1.
a13 = 2.5
a21 = 1.5
a31 = 0.2
r2 = 0.6
r3 = 4.5
k3 = 1
d3 = 0.5


def initialize(x0, y0, z0):
    global x, xresult, y, yresult, z, zresult, t, timesteps
    x = x0
    y = y0
    z = z0
    xresult = [x]
    yresult = [y]
    zresult = [z]
    t = 0.
    timesteps = [t]

def observe():
    global x, xresult, y, yresult, z, zresult, t, timesteps
    xresult.append(x)
    yresult.append(y)
    zresult.append(z)
    timesteps.append(t)

def update():
    global x, xresult, y, yresult, z, zresult, t, timesteps

    nextx = x + (x - x*x - a12*x*y - a13*x*z) * Dt
    nexty = y + (r2*y - r2*y*y - a21*x*y) * Dt
    nextz = z + ((r3*x*z)/(x+k3) - a31*x*z - d3*z) * Dt

    x, y, z = nextx, nexty, nextz
    t = t + Dt

##### Simulate
#initialize(0.5, 0.1, 0.1)

initialize(0.1, 0.1, 0.1)
while t < 2000:
  update()
  observe()


#### Plotting
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(2,2,1, projection='3d')
ax.plot(xresult, yresult, zresult, "b")
ax.view_init(30, 45)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
#ax.set_title("Phase Space")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_zlim(0, 1)
ax1 = fig.add_subplot(2,2,2)
ax1.plot(xresult, yresult, "b")
xlabel("x")
ylabel("y")
plt.xlim([0, 1])
plt.ylim([0, 1])
ax1.grid(True)
ax1 = fig.add_subplot(2,2,3)
ax1.plot(xresult, zresult, "b")
xlabel("x")
ylabel("z")
plt.xlim([0, 1])
plt.ylim([0, 1])
ax1.grid(True)
ax1 = fig.add_subplot(2,2,4)
ax1.plot(yresult, zresult, "b")
xlabel("y")
ylabel("z")
plt.xlim([0, 1])
plt.ylim([0, 1])
ax1.grid(True)
fig.tight_layout()
#plt.savefig('phases_c.eps', format='eps')
plt.show()
