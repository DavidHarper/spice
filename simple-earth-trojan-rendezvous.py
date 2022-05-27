#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, pi, sin, cos
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if len(sys.argv) < 7:
    print("Usage: ", sys.argv[0], " trojan-mass/mEarth r theta dvx dvy duration")
    quit()

spice.furnsh('kernels/gm_de431.tpc')

trojanMass = float(sys.argv[1])
r = float(sys.argv[2])
theta = float(sys.argv[3]) * pi/180.0
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

pxElt = c60 * au
pyElt = s60 * au

vxElt = -s60 * au * nEarth
vyElt =  c60 * au * nEarth

earthLeadingTrojan = rebound.Particle(m = trojanMass*mEarth, x=pxElt, y=pyElt, z=0.0, vx=vxElt, vy=vyElt, vz=0.0 )

sim.add(earthLeadingTrojan)

pxProbe = au + r * cos(theta)
pyProbe = r * sin(theta)

vxProbe = dvx
vyProbe = au*nEarth + dvy

probe = rebound.Particle(m=0.0, x=pxProbe, y=pyProbe, z=0.0, vx=vxProbe, vy=vyProbe, vz=0.0)

sim.add(probe)

tdata = []
xTrojan = []
yTrojan = []
rTrojan = []
xProbe = []
yProbe = []
rProbeSun = []
rProbeTrojan = []
vrProbeTrojan = []
rProbeEarth = []
xProbeAbs = []
yProbeAbs = []
eProbe = []
aProbe = []

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    #print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    tjd = sim.t/86400.0

    sun = sim.particles[0]
    earth = sim.particles[1]
    trojan = sim.particles[2]
    probe = sim.particles[3]

    pEarth = np.array([earth.x - sun.x, earth.y - sun.y, earth.z - sun.z], dtype = float)
    vEarth = np.array([earth.vx - sun.vx, earth.vy - sun.vy, earth.vz - sun.vz], dtype = float)

    u = spice.vhat(pEarth)
    w = spice.vhat(spice.vcrss(pEarth, vEarth))
    v = spice.vhat(spice.vcrss(w, u))

    pTrojan = np.array([trojan.x - sun.x, trojan.y - sun.y, trojan.z - sun.z], dtype = float)
    vTrojan = np.array([trojan.vx - sun.vx, trojan.vy - sun.vy, trojan.vz - sun.vz], dtype = float)

    dxTrojan = spice.vdot(u, pTrojan) - c60 * au
    dyTrojan = spice.vdot(v, pTrojan) - s60 * au
    dzTrojan = spice.vdot(w, pTrojan)

    pProbe = np.array([probe.x - sun.x, probe.y - sun.y, probe.z - sun.z], dtype = float)

    delta = spice.vnorm(spice.vsub(pProbe, pEarth))

    if delta < 250000.0:
        return

    dxProbe = spice.vdot(u, pProbe) - c60 * au
    dyProbe = spice.vdot(v, pProbe) - s60 * au

    tdata.append(tjd)

    xTrojan.append(dxTrojan/au)
    yTrojan.append(dyTrojan/au)
    rTrojan.append(spice.vnorm(spice.vsub(pEarth, pTrojan))/au)

    xProbe.append(dxProbe/au)
    yProbe.append(dyProbe/au)
    rProbeTrojan.append(spice.vnorm(spice.vsub(pProbe, pTrojan))/au)
    rProbeEarth.append(delta/au)
    rProbeSun.append(spice.vnorm(pProbe)/au)

    rrProbeTrojan = np.array([probe.x - trojan.x, probe.y - trojan.y, probe.z - trojan.z], dtype = float)
    vProbeTrojan = np.array([probe.vx - trojan.vx, probe.vy - trojan.vy, probe.vz - trojan.vz], dtype = float)
    vrProbeTrojan.append(spice.vdot(spice.vhat(rrProbeTrojan), vProbeTrojan))

    xProbeAbs.append(probe.x/au)
    yProbeAbs.append(probe.y/au)

    pvProbe = np.array([probe.x, probe.y, probe.z, probe.vx, probe.vy, probe.vz], dtype=float)
    pvSun = np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz], dtype=float)

    oscElt = spice.oscltx(spice.vsubg(pvProbe, pvSun), sim.t, gmSun)

    eProbe.append(oscElt[1])
    aProbe.append(oscElt[9]/au)

sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0 * 0.05

sim.integrate(86400.0 * days)

fig = plt.figure(tight_layout=True)
gs = gridspec.GridSpec(5, 2)

axDistances = fig.add_subplot(gs[0, :])

axDistances.plot(tdata, rProbeTrojan, label='r[Probe-Trojan]')
axDistances.plot(tdata, rProbeSun, label='r[Probe-Sun]')
axDistances.plot(tdata, rProbeEarth, label='r[Probe-Earth]')
axDistances.legend()

axRelVel = fig.add_subplot(gs[1, :])
axRelVel.plot(tdata, vrProbeTrojan, label='v[Probe-Trojan]')
axRelVel.legend()

axSMAProbe = fig.add_subplot(gs[2, :])
axSMAProbe.plot(tdata, aProbe, label='a[Probe]')
axSMAProbe.legend()

axEccProbe = fig.add_subplot(gs[3, :])
axEccProbe.plot(tdata, eProbe, label='e[Probe]')
axEccProbe.legend()

axTrojanFrame = fig.add_subplot(gs[4, 0])

axTrojanFrame.plot([0], [0], "or")
axTrojanFrame.plot(xTrojan, yTrojan, label='Trojan')
axTrojanFrame.plot(xProbe, yProbe, label='Probe')
axTrojanFrame.plot([1.0-c60],[-s60], "og")
axTrojanFrame.plot([-c60], [-s60], "oy")
axTrojanFrame.set(aspect=1)

axAbsoluteFrame = fig.add_subplot(gs[4, 1])

axAbsoluteFrame.plot([0], [0], "oy")
axAbsoluteFrame.plot(xProbeAbs, yProbeAbs, label="Probe")

N = 400
t = np.linspace(0, 2 * pi, N)
xEarth, yEarth = np.cos(t), np.sin(t)
axAbsoluteFrame.plot(xEarth, yEarth, "r", label="Earth orbit")
axAbsoluteFrame.set(aspect=1)

plt.show()
