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

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], " start-date duration trojan-mass/mEarth")
    quit()

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

t0 = spice.str2et(sys.argv[1])
days = float(sys.argv[2])
trojanMass = float(sys.argv[3])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)

gmSun = data[1][0]

mSun = gmSun/G

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

planets = ['MERCURY BARYCENTER', 'VENUS BARYCENTER', 'EARTH BARYCENTER',
'MARS BARYCENTER', 'JUPITER BARYCENTER', 'SATURN BARYCENTER',
'URANUS BARYCENTER', 'NEPTUNE BARYCENTER']

for planet in planets:
    data = spice.bodvrd(planet, 'GM', 1)
    mPlanet = data[1][0]/G
    state=spice.spkezr(planet, t0, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    planet = rebound.Particle(m=mPlanet, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
    sim.add(planet)

# Create the leading Trojan Earth
planet = 'EARTH BARYCENTER'
data = spice.bodvrd(planet, 'GM', 1)
mEarth = data[1][0]/G
state=spice.spkezr(planet, t0, 'J2000', 'NONE', 'SUN')
pv = state[0]

u = spice.vhat(pv[0:3])
w = spice.vhat(spice.vcrss(pv[0:3], pv[3:6]))
v = spice.vhat(spice.vcrss(w, u))

px = spice.vdot(u, pv[0:3])
py = spice.vdot(v, pv[0:3])

vx = spice.vdot(u, pv[3:6])
vy = spice.vdot(v, pv[3:6])

c60 = 0.5
s60 = sqrt(3.0)/2.0

pxElt = c60 * px + s60 * py
pyElt = -s60 * px + c60 * py

vxElt = c60 * vx + s60 * vy
vyElt = -s60 * vx + c60 * vy

pElt = spice.vlcom(pxElt, u, pyElt, v)
vElt = spice.vlcom(vxElt, u, vyElt, v)

earthLeadingTrojan = rebound.Particle(m = trojanMass*mEarth, x=pElt[0], y=pElt[1], z=pElt[2], vx=vElt[0], vy=vElt[1], vz=vElt[2] )
sim.add(earthLeadingTrojan)

sim.t = t0
sim.dt = 86400.0 * 0.05

def distance(body1, body2):
    dx = body1.x - body2.x
    dy = body1.y - body2.y
    dz = body1.z - body2.z
    return sqrt(dx*dx+dy*dy+dz*dz)

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    #print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    djd = sim.t/86400.0

    kSun = 0
    sun = sim.particles[kSun]

    kEarth = 3
    earth = sim.particles[kEarth]

    kElt = 9
    elt = sim.particles[kElt]

    cx = earth.x - sun.x
    cy = earth.y - sun.y
    cz = earth.z - sun.z

    state=spice.spkezr(planets[kEarth-1], sim.t, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    vr = spice.vhat(pv[0:3])
    vz = spice.vhat(spice.vcrss(pv[0:3], pv[3:6]))
    vt = spice.vhat(spice.vcrss(vz, vr))
    dx = cx - pv[0]
    dy = cy - pv[1]
    dz = cz - pv[2]
    dp = np.array([dx, dy, dz], dtype = float)
    print('A %14.6f %2d %14.3f %14.3f %14.3f' % (djd, kEarth, spice.vdot(vr, dp), spice.vdot(vt, dp), spice.vdot(vz, dp)))

    pvEarth = np.array([earth.x, earth.y, earth.z, earth.vx, earth.vy, earth.vz], dtype=float)
    pvElt = np.array([elt.x, elt.y, elt.z, elt.vx, elt.vy, elt.vz], dtype=float)
    pvSun = np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz], dtype=float)

    oscEarth = spice.oscelt(spice.vsubg(pvEarth, pvSun), sim.t, gmSun)
    oscElt = spice.oscelt(spice.vsubg(pvElt, pvSun), sim.t, gmSun)

    dtheta = spice.vsep(pvEarth[0:3], pvElt[0:3]) * 180.0/pi

    print('B %14.6f %14.8f %14.8f %14.6f' % (djd, oscEarth[1], oscElt[1], dtheta))

sim.heartbeat=heartbeat

sim.integrate(86400.0 * days)
