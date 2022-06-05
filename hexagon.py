#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, cos, sin, pi
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

tdata = []
xdata = []
ydata = []
rdata = []

for theta in [60, 120, 180, 240, 300]:
    ctheta = cos(theta*pi/180.0)
    stheta = sin(theta*pi/180.0)

    pxElt = ctheta * au + dx
    pyElt = stheta * au + dy

    vxElt = -stheta * au * nEarth + dvx
    vyElt =  ctheta * au * nEarth + dvy

    body = rebound.Particle(m = trojanMass*mEarth, x=pxElt, y=pyElt, z=0.0, vx=vxElt, vy=vyElt, vz=0.0 )

    sim.add(body)

    xdata.append([])
    ydata.append([])
    rdata.append([])

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    #print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    tjd = sim.t/86400.0

    sun = sim.particles[0]
    earth = sim.particles[1]

    pEarth = np.array([earth.x - sun.x, earth.y - sun.y, earth.z - sun.z], dtype = float)
    vEarth = np.array([earth.vx - sun.vx, earth.vy - sun.vy, earth.vz - sun.vz], dtype = float)

    u = spice.vhat(pEarth)
    w = spice.vhat(spice.vcrss(pEarth, vEarth))
    v = spice.vhat(spice.vcrss(w, u))

    tdata.append(tjd)

    for k in [2, 3, 4, 5, 6]:
        trojan = sim.particles[k]

        pTrojan = np.array([trojan.x - sun.x, trojan.y - sun.y, trojan.z - sun.z], dtype = float)

        dxTrojan = spice.vdot(u, pTrojan)
        dyTrojan = spice.vdot(v, pTrojan)

        xdata[k-2].append(dxTrojan/au)
        ydata[k-2].append(dyTrojan/au)
        rdata[k-2].append(spice.vnorm(spice.vsub(pEarth, pTrojan))/au)


sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0 * 0.05

sim.integrate(86400.0 * days)

fig, ax = plt.subplots()

ax.plot([0], [0], "or")
ax.plot([1], [0], "ob")

for k in [0, 1, 2, 3, 4]:
    ax.plot(xdata[k], ydata[k], label=str(k))

ax.legend()
ax.set(aspect=1)

plt.show()
