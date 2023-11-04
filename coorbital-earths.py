#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, pi, sin, cos, atan2, fabs, log10, asin
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time

if len(sys.argv) < 5:
    print("Usage: ", sys.argv[0], " a[Earth1] a[Earth2] mu[Earth1] mu[Earth2] days")
    quit()

spice.furnsh('kernels/gm_de431.tpc')

aEarth1 = float(sys.argv[1])
aEarth2 = float(sys.argv[2])
muEarth1 = float(sys.argv[3])
muEarth2 = float(sys.argv[4])
days = float(sys.argv[5])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)
gmSun = data[1][0]
mSun = gmSun/G

data = spice.bodvrd('EARTH BARYCENTER','GM',1)
gmEarth = data[1][0]
mEarth = gmEarth/G

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

gmSystem = gmSun + gmEarth * (muEarth1 + muEarth2)

nEarth1 = sqrt(gmSystem/aEarth1**3)
vEarth1 = nEarth1 * aEarth1

earth1 = rebound.Particle(m=mEarth*muEarth1, x=aEarth1, y=0.0, z=0.0, vx=0.0, vy=vEarth1, vz=0.0)

sim.add(earth1)

nEarth2 = sqrt(gmSystem/aEarth2**3)
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
elEarth2 = []
magEarth2 = []
diaEarth2 = []

LOG_AU = log10(149600000.0)

def magnitude(r, d, p):
    theta = p/100.0
    dm = 5.0 * (log10(r*d) - 2.0 * LOG_AU)
    return -4.40 + dm + 0.09 * theta + 2.39 * theta * theta - 0.65 * theta * theta * theta

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

    r2 = spice.vnorm(pEarth2)
    rEarth2.append(r2 - aEarth2)

    pEarth2Earth1 = spice.vsub(pEarth2, pEarth1)
    delta = spice.vnorm(pEarth2Earth1)
    rEarth2Earth1.append(delta)

    xEarth2Earth1.append(spice.vdot(u, pEarth2Earth1))
    yEarth2Earth1.append(spice.vdot(v, pEarth2Earth1))

    dx = spice.vdot(u, pEarth2)
    dy = spice.vdot(v, pEarth2)
    qEarth2 = 180.0*atan2(dy, dx)/pi
    qEarth2Earth1.append(qEarth2)
    phaseEarth2 = 90.0-fabs(0.5*qEarth2)

    elEarth2.append(phaseEarth2)

    mag = magnitude(r2, delta, phaseEarth2)
    magEarth2.append(mag)

    dia = 2.0 * asin(6378.0/delta) * 3600.0 * 180.0/pi
    diaEarth2.append(dia)

sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0

ns_start = time.time_ns()

sim.integrate(86400.0 * days)

ns_end = time.time_ns()

ns_run = ns_end - ns_start

steps = len(elEarth2)

ns_per_step = ns_run/steps

print(f"Ran {steps} steps in {ns_run:,} ns which is {ns_per_step:0f} ns per step")

fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(7, 2)

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

axElong = fig.add_subplot(gs[3, :])
axElong.plot(tdata, elEarth2, label='el[Earth2]')
axElong.legend()

axMag = fig.add_subplot(gs[4, :])
axMag.plot(tdata, magEarth2, label='mag[Earth2]')
axMag.legend()

axDia = fig.add_subplot(gs[5, :])
axDia.plot(tdata, diaEarth2, label='dia[Earth2]')
axDia.legend()

axEarth1Frame = fig.add_subplot(gs[6, 0])
axEarth1Frame.plot([0], [0], "or")
axEarth1Frame.plot(xEarth2Earth1, yEarth2Earth1, label='Earth2')
axEarth1Frame.set(aspect=1)


plt.show()
