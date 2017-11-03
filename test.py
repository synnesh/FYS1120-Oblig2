import numpy as np
from math import ceil, cos, sqrt, pi
import matplotlib.pyplot as plt

Bx = 0    #NC^-1
By = 0
Bz = 1.5

q = -1.6E-19
me = 9.11E-31
mp = 1.67E-27
g = 9.81
E0 =float(25e3)/90e-6    #V/m
d = 0.1                  #diameter i meter
gap= 90e-6
B = sqrt(Bx**2+By**2+Bz**2)
omega = float(q*B)/mp
T = 2*pi/omega

t0 = 0
x0 = 0
y0 = 0
z0 = 0

v0x = 0 #m/s
v0y = 0
v0z = 0

time = 300.0e-9 #s
dt = 100e-15

#Variables
n = ceil(float(time)/dt)
x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)
vx = np.zeros(n)
vy = np.zeros(n)
vz = np.zeros(n)
ax = np.zeros(n)
ay = np.zeros(n)
az = np.zeros(n)
t = np.zeros(n)

#Initial conditions
x[0] = x0
y[0] = y0
z[0] = z0
vx[0] = v0x
vy[0] = v0y
vz[0] = v0z

def E(xi):
    if (xi>= -gap and xi<= gap):
        Ex = E0*cos(omega*T)
    else:
        Ex = 0
    return Ex
Ey = 0
Ez = 0
#Integration loop:
for i in range(int(n-1)):
    r = sqrt(x[i]**2+y[i]**2)
    if(r == 0.00001):
        break
    if (E(x[i]) != 0):
        Fx = E(x[i])*q
        Fy = Ey*q
        Fz = Ez*q
    else:
        Fx = q*(vy[i]*Bz - vz[i]*By)
        Fy = -q*(vx[i]*Bz - vz[i]*Bx)
        Fz = q*(vx[i]*By - vy[i]*Bx)
    # x
    ax[i]    = Fx/me
    vx[i+1]  = vx[i] + ax[i]*dt
    x[i+1]  = x[i] + vx[i+1]*dt
    t[i+1]  = t[i] + dt

    # y
    ay[i]    = Fy/me
    vy[i+1]  = vy[i] + ay[i]*dt
    y[i+1]  = y[i] + vy[i+1]*dt

    # z
    az[i] = Fz/me
    vz[i+1] = vz[i] + az[i]*dt
    z[i+1] = z[i] + vz[i+1]*dt

#Plot
plt.figure(1)
plt.title("Partikkel i syklotron")
plt.ylabel("y [m/s]")
plt.xlabel("x [m/s]")
plt.plot(y,x)
#plt.legend(["y","x"])
plt.grid()
plt.show()
