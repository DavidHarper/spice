#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, pi, sin, cos, atan2
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], " a[Earth1] a[Earth2] M[Earth2]/M[Earth1] days")
    quit()

spice.furnsh('kernels/gm_de431.tpc')

aEarth1 = float(sys.argv[1])
aEarth2 = float(sys.argv[2])
muEarth2 = float(sys.argv[3])
days = float(sys.argv[4])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)
gmSun = data[1][0]
mSun = gmSun/G

data = spice.bodvrd('EARTH BARYCENTER','GM',1)
mEarth = data[1][0]/G

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

nEarth1 = sqrt(gmSun/aEarth1**3)
vEarth1 = nEarth1 * aEarth1

earth1 = rebound.Particle(m=mEarth, x=aEarth1, y=0.0, z=0.0, vx=0.0, vy=vEarth1, vz=0.0)

sim.add(earth1)

nEarth2 = sqrt(gmSun/aEarth2**3)
vEarth2 = nEarth2 * aEarth2

earth2 = rebound.Particle(m=mEarth*muEarth2, x=-aEarth2, y=0.0, z=0.0, vx=0.0, vy=-vEarth2, vz=0.0)

sim.add(earth2)

tdata = []

rEarth1 = []
rEarth2 = []
rEarth2Earth1 = []
xEarth2Earth1 = []
yEarth2Earth1 = []
qEarth2Earth1 = []

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    tjd = sim.t/86400.0

    sun = sim.particles[0]
    earth1 = sim.particles[1]
    earth2 = sim.particles[2]

    pEarth1 = np.array([earth1.x - sun.x, earth1.y - sun.y, earth1.z - sun.z], dtype = float)
    vEarth1 = np.array([earth1.vx - sun.vx, earth1.vy - sun.vy, earth1.vz - sun.vz], dtype = float)

    u = spice.vhat(pEarth1)
    w = spice.vhat(spice.vcrss(pEarth1, vEarth1))
    v = spice.vhat(spice.vcrss(w, u))

    pEarth2 = np.array([earth2.x - sun.x, earth2.y - sun.y, earth2.z - sun.z], dtype = float)

    tdata.append(tjd/365.25)

    rEarth1.append(spice.vnorm(pEarth1) - aEarth1)
    rEarth2.append(spice.vnorm(pEarth2) - aEarth2)

    pEarth2Earth1 = spice.vsub(pEarth2, pEarth1)
    rEarth2Earth1.append(spice.vnorm(pEarth2Earth1))

    xEarth2Earth1.append(spice.vdot(u, pEarth2Earth1))
    yEarth2Earth1.append(spice.vdot(v, pEarth2Earth1))

    dx = spice.vdot(u, pEarth2)
    dy = spice.vdot(v, pEarth2)
    qEarth2Earth1.append(180.0*atan2(dy, dx)/pi)

sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0 * 0.1

sim.integrate(86400.0 * days)

fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(4, 2)

axDistances = fig.add_subplot(gs[0, :])
axDistances.plot(tdata, rEarth1, label='r[Earth1]')
axDistances.plot(tdata, rEarth2, label='r[earth2]')
axDistances.legend()

axDelta = fig.add_subplot(gs[1, :])
axDelta.plot(tdata, rEarth2Earth1, label='r[earth2-Earth1]')
axDelta.legend()

axTheta = fig.add_subplot(gs[2, :])
axTheta.plot(tdata, qEarth2Earth1, label='q[earth2-Earth1]')
axTheta.legend()

axEarth1Frame = fig.add_subplot(gs[3, 0])
axEarth1Frame.plot([0], [0], "or")
axEarth1Frame.plot(xEarth2Earth1, yEarth2Earth1, label='Earth2')
axEarth1Frame.set(aspect=1)


plt.show()
