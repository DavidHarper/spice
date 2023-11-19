#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, pi, sin, cos
import sys
import spiceypy as spice
import rebound
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], " altitude phase delta-V days")
    quit()

spice.furnsh('kernels/gm_de431.tpc')

altitude = float(sys.argv[1])
phase = float(sys.argv[2]) * pi/180.0
deltav = float(sys.argv[3])
days = float(sys.argv[4])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('EARTH','GM',1)
gmEarth = data[1][0]
mEarth = gmEarth/G

aeEarth = 6378

data = spice.bodvrd('MOON','GM',1)
gmMoon = data[1][0]
mMoon = gmEarth/G

aeMoon = 1738

aMoon = 384400

sim = rebound.Simulation()

sim.G = G

earth = rebound.Particle(m=mEarth)

sim.add(earth)

gmSystem = gmEarth + gmMoon

nMoon = sqrt(gmSystem/aMoon**3)
vMoon = nMoon * aMoon

moon = rebound.Particle(m=mMoon, x=aMoon, y=0.0, z=0.0, vx=0.0, vy=vMoon, vz=0.0)

sim.add(moon)

aApollo0 = aeEarth + altitude
nApollo0 = sqrt(gmEarth/aApollo0**3)
vApollo0 = nApollo0 * aApollo0

x0Apollo = aApollo0 * cos(phase)
y0Apollo = aApollo0 * sin(phase)

aApollo1 = 0.5 * ( aeEarth + altitude + aMoon )
vApollo1 = sqrt(gmEarth * (2.0/aApollo0 - 1.0/aApollo1))

print("At altitude %.2f km, circular velocity is %.3f km/s, Hohmann transfer velocity is %.3f km/s" % (altitude, vApollo0, vApollo1))

vApollo1 = vApollo1 + deltav

vx0Apollo = - vApollo1 * sin(phase)
vy0Apollo =   vApollo1 * cos(phase)

apollo = rebound.Particle(m=0.0, x=x0Apollo, y=y0Apollo, z=0.0, vx=vx0Apollo, vy=vy0Apollo, vz=0.0)

sim.add(apollo)

tdata = []

xApolloRotating = []
yApolloRotating = []
xApolloFixed = []
yApolloFixed = []
xApolloMoon = []
yApolloMoon = []

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    tjd = sim.t/86400.0

    earth = sim.particles[0]
    moon = sim.particles[1]
    apollo = sim.particles[2]

    pMoon = np.array([moon.x - earth.x, moon.y - earth.y, moon.z - earth.z], dtype = float)
    vMoon = np.array([moon.vx - earth.vx, moon.vy - earth.vy, moon.vz - earth.vz], dtype = float)

    u = spice.vhat(pMoon)
    w = spice.vhat(spice.vcrss(pMoon, vMoon))
    v = spice.vhat(spice.vcrss(w, u))

    pApollo= np.array([apollo.x - earth.x, apollo.y - earth.y, apollo.z - earth.z], dtype = float)
    pApolloMoon = np.array([apollo.x - moon.x, apollo.y - moon.y, apollo.z - moon.z], dtype = float)

    tdata.append(tjd)

    xApolloFixed.append(apollo.x - earth.x)
    yApolloFixed.append(apollo.y - earth.y)

    xa = spice.vdot(u, pApollo)
    ya = spice.vdot(v, pApollo)

    xApolloRotating.append(xa)
    yApolloRotating.append(ya)

    xApolloMoon.append(spice.vdot(u, pApolloMoon))
    yApolloMoon.append(spice.vdot(v, pApolloMoon))

sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0

ns_start = time.time_ns()

sim.integrate(86400.0 * days)

ns_end = time.time_ns()

ns_run = ns_end - ns_start

steps = len(tdata)

ns_per_step = ns_run/steps

print(f"Ran {steps} steps in {ns_run:,} ns which is {ns_per_step:0f} ns per step")

fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(3, 1)

axApollo1 = fig.add_subplot(gs[0, 0])
axApollo1.plot([0], [0], "or")
axApollo1.plot(xApolloFixed, yApolloFixed, label='Earth centred, Fixed frame')
axApollo1.set(aspect=1)

axApollo2 = fig.add_subplot(gs[1, 0])
axApollo2.plot([0], [0], "or")
axApollo2.plot(xApolloRotating, yApolloRotating, label='Earth centred, Rotating frame')
axApollo2.set(aspect=1)

axApollo3 = fig.add_subplot(gs[2, 0])
axApollo3.plot([0], [0], "og")
axApollo3.plot(xApolloMoon, yApolloMoon, label='Moon centred, rotating frame')
axApollo3.set(aspect=1)

plt.show()