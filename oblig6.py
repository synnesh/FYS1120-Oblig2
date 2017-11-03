import matplotlib.pyplot as plt
import numpy as np
from math import ceil, sqrt

Bx = 0    #NC^-1
By = 0
Bz = 2.0

q = -1.6E-19
me = 9.11E-31
g = 9.81

t0 = 0
x0 = 0
y0 = 0
z0 = 0

v0x = 10e3 #m/s
v0y = 0
v0z = 0

time = 30.0e-12 #s
dt = 1e-15

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

#Integration loop
for i in range(int(n-1)):
    # x
    Fx = q*(vy[i]*Bz - vz[i]*By)
    ax[i]    = Fx/me
    vx[i+1]  = vx[i] + ax[i]*dt
    x[i+1]  = x[i] + vx[i+1]*dt
    t[i+1]  = t[i] + dt
    # y
    Fy = -q*(vx[i]*Bz - vz[i]*Bx)
    ay[i]    = Fy/me
    vy[i+1]  = vy[i] + ay[i]*dt
    y[i+1]  = y[i] + vy[i+1]*dt
    # z
    Fz = q*(vx[i]*By - vy[i]*Bx)
    az[i] = Fz/me
    vz[i+1] = vz[i] + az[i]*dt
    z[i+1] = z[i] + vz[i+1]*dt

#a = Ex*q/me
#r = 0.5*a*t**2

#Plot
plt.figure(1)
plt.title("Elektron i magnetisk felt")
plt.ylabel("fart [m/s]")
plt.xlabel("tid [s]")
plt.plot(t,vx,t,vy,t,vz)
plt.legend(["vx","vy","vz"])
plt.grid()
plt.show()
