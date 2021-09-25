from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

Dt = 0.001 # Size of time step
a12 = 1 # Chaotic behaviour if a12 = 1.
a13 = 2.5
a21 = 1.5
a31 = 0.2
r2 = 0.6
r3 = 4.5
k3 = 1
d3 = 0.5
p = 1.1

def initialize(x0, y0, z0):
    global x, xresult, y, yresult, z, zresult, t, timesteps, u, uresult
    x = x0
    y = y0
    z = z0
    xresult = [x]
    yresult = [y]
    zresult = [z]

    u = 0.2
    uresult = [u]

    t = 0.
    timesteps = [t]

def observe():
    global x, xresult, y, yresult, z, zresult, t, timesteps, u, uresult
    xresult.append(x)
    yresult.append(y)
    zresult.append(z)
    uresult.append(u)
    timesteps.append(t)

def update():
    global x, xresult, y, yresult, z, zresult, t, timesteps, u, uresult

    nextx = x + (x - x*x - a12*x*y - a13*x*z - x*(0.15*(1-exp(-u)))) * Dt
    nexty = y + (r2*y - r2*y*y - a21*x*y) * Dt
    nextz = z + ((r3*x*z)/(x+k3) - a31*x*z - d3*z) * Dt
    nextu = u + (0.01-1.5*u*x) * Dt

    x, y, z, u = nextx, nexty, nextz, nextu
    t = t + Dt

####### PLOT
initialize(0.1, 0.1, 0.1)
while t < 2000:
  update()
  observe()
fig=plt.figure(figsize=(12,6))
plt.subplot(411)
plt.plot(timesteps, xresult)
plt.xlim([0, 2000])
plt.ylim([0, 1])
ylabel("x")
plt.title("Time responses of the system states")
plt.subplot(412)
plt.plot(timesteps, yresult)
plt.xlim([0, 2000])
plt.ylim([0, 1])
ylabel("y")
plt.subplot(413)
plt.plot(timesteps, zresult)
plt.xlim([0, 2000])
plt.ylim([0, 1])
xlabel("Time")
ylabel("z")
plt.subplot(414)
plt.plot(timesteps, uresult)
plt.xlim([0, 2000])
plt.ylim([0, 1])
xlabel("Time")
ylabel("u")
plt.tight_layout()
#plt.savefig('more.eps', format='eps')
plt.show()
