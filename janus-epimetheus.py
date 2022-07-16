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

if len(sys.argv) < 5:
    print("Usage: ", sys.argv[0], " a[Janus] mu[Janus] a[Epimetheus] mu[Epimetheus] days")
    quit()

spice.furnsh('kernels/gm_de431.tpc')

aJanus = float(sys.argv[1])
muJanus = float(sys.argv[2])
aEpimetheus = float(sys.argv[3])
muEpimetheus = float(sys.argv[4])
days = float(sys.argv[5])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SATURN BARYCENTER','GM',1)
gmSaturn = data[1][0]
mSaturn = gmSaturn/G

data = spice.bodvrd('JANUS','GM',1)
mJanus = data[1][0]/G

data = spice.bodvrd('EPIMETHEUS','GM',1)
mEpimetheus = data[1][0]/G

sim = rebound.Simulation()

sim.G = G

saturn = rebound.Particle(m=mSaturn)

sim.add(saturn)

nJanus = sqrt(gmSaturn/aJanus**3)
vJanus = nJanus * aJanus

janus = rebound.Particle(m=mJanus, x=aJanus, y=0.0, z=0.0, vx=0.0, vy=vJanus, vz=0.0)

sim.add(janus)

nEpimetheus = sqrt(gmSaturn/aEpimetheus**3)
vEpimetheus = nEpimetheus * aEpimetheus

epimetheus = rebound.Particle(m=mEpimetheus, x=-aEpimetheus, y=0.0, z=0.0, vx=0.0, vy=-vEpimetheus, vz=0.0)

sim.add(epimetheus)

tdata = []

rJanus = []
rEpimetheus = []
rEpimetheusJanus = []
xEpimetheus = []
yEpimetheus = []
qEpimetheus = []

def heartbeat2(sim_pointer):
    sim = sim_pointer.contents

    saturn = sim.particles[0]
    janus = sim.particles[1]
    epimetheus = sim.particles[2]

    print('%14.6f %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e' % (
    sim.t, janus.x, janus.y, janus.z, epimetheus.x, epimetheus.y, epimetheus.z))

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    tjd = sim.t/86400.0

    saturn = sim.particles[0]
    janus = sim.particles[1]
    epimetheus = sim.particles[2]

    pJanus = np.array([janus.x - saturn.x, janus.y - saturn.y, janus.z - saturn.z], dtype = float)
    vJanus = np.array([janus.vx - saturn.vx, janus.vy - saturn.vy, janus.vz - saturn.vz], dtype = float)

    u = spice.vhat(pJanus)
    w = spice.vhat(spice.vcrss(pJanus, vJanus))
    v = spice.vhat(spice.vcrss(w, u))

    pEpimetheus = np.array([epimetheus.x - saturn.x, epimetheus.y - saturn.y, epimetheus.z - saturn.z], dtype = float)

    tdata.append(tjd)

    rJanus.append(spice.vnorm(pJanus) - aJanus)
    rEpimetheus.append(spice.vnorm(pEpimetheus) - aEpimetheus)

    pEpimetheusJanus = spice.vsub(pEpimetheus, pJanus)
    rEpimetheusJanus.append(spice.vnorm(pEpimetheusJanus))

    xEpimetheus.append(spice.vdot(u, pEpimetheusJanus))
    yEpimetheus.append(spice.vdot(v, pEpimetheusJanus))

    qEpimetheus.append(190.0*atan2(spice.vdot(v, pEpimetheus), spice.vdot(u, pEpimetheus))/pi)

sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0 * 0.1

sim.integrate(86400.0 * days)

fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(4, 2)

axDistances = fig.add_subplot(gs[0, :])

axDistances.plot(tdata, rJanus, label='r[Janus]')
axDistances.plot(tdata, rEpimetheus, label='r[Epimetheus]')
axDistances.legend()

axDelta = fig.add_subplot(gs[1, :])
axDelta.plot(tdata, rEpimetheusJanus, label='r[Epimetheus-Janus]')
axDelta.legend()

axTheta = fig.add_subplot(gs[2, :])
axTheta.plot(tdata, qEpimetheus, label='q[Epimetheus]')
axTheta.legend()

axJanusFrame = fig.add_subplot(gs[3, 0])
axJanusFrame.plot([0], [0], "or")
axJanusFrame.plot(xEpimetheus, yEpimetheus, label='Epimetheus')
axJanusFrame.set(aspect=1)

plt.show()
