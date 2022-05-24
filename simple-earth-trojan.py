#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, pi
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 7:
    print("Usage: ", sys.argv[0], " trojan-mass/mEarth dx dy dvx dvy duration")
    quit()

spice.furnsh('kernels/gm_de431.tpc')

trojanMass = float(sys.argv[1])
dx = float(sys.argv[2])
dy = float(sys.argv[3])
dvx = float(sys.argv[4])
dvy = float(sys.argv[5])
days = float(sys.argv[6])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)
gmSun = data[1][0]
mSun = gmSun/G

planet = 'EARTH BARYCENTER'
data = spice.bodvrd(planet, 'GM', 1)
gmEarth = data[1][0]
mEarth = gmEarth/G

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

au = spice.convrt(1.0,'AU','KM')

nEarth = sqrt((gmSun + gmEarth)/au**3)

pEarth = np.array([au, 0.0, 0.0], dtype = float)
vEarth = np.array([0.0, au*nEarth, 0.0], dtype = float)

earth = rebound.Particle(m=mEarth, x=pEarth[0], y=pEarth[1], z=pEarth[2], vx=vEarth[0], vy=vEarth[1], vz=vEarth[2])

sim.add(earth)

c60 = 0.5
s60 = sqrt(3.0)/2.0

pxElt = c60 * au + dx
pyElt = s60 * au + dy

vxElt = -s60 * au * nEarth + dvx
vyElt =  c60 * au * nEarth + dvy

earthLeadingTrojan = rebound.Particle(m = trojanMass*mEarth, x=pxElt, y=pyElt, z=0.0, vx=vxElt, vy=vyElt, vz=0.0 )

sim.add(earthLeadingTrojan)

tdata = []
xdata = []
ydata = []

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    #print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    tjd = sim.t/86400.0

    sun = sim.particles[0]
    earth = sim.particles[1]
    trojan = sim.particles[2]

    pEarth = np.array([earth.x - sun.x, earth.y - sun.y, earth.z - sun.z], dtype = float)
    vEarth = np.array([earth.vx - sun.vx, earth.vy - sun.vy, earth.vz - sun.vz], dtype = float)

    u = spice.vhat(pEarth)
    w = spice.vhat(spice.vcrss(pEarth, vEarth))
    v = spice.vhat(spice.vcrss(w, u))

    pTrojan = np.array([trojan.x - sun.x, trojan.y - sun.y, trojan.z - sun.z], dtype = float)

    dxTrojan = spice.vdot(u, pTrojan) - c60 * au
    dyTrojan = spice.vdot(v, pTrojan) - s60 * au
    dzTrojan = spice.vdot(w, pTrojan)

    print('%14.6f %14.3f %14.3f %14.3f' % (tjd, dxTrojan, dyTrojan, dzTrojan))

    tdata.append(tjd)
    xdata.append(dxTrojan/au)
    ydata.append(dyTrojan/au)


sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0 * 0.05

sim.integrate(86400.0 * days)

fig, ax = plt.subplots(2)

ax[0].plot(tdata, xdata, label='x')
ax[0].plot(tdata, ydata, label='y')
ax[0].legend()

ax[1].plot([0], [0], "or")
ax[1].plot(xdata, ydata)
ax[1].plot([1.0-c60],[-s60], "og")
ax[1].set(aspect=1)

plt.show()
