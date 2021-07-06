#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np

if len(sys.argv) < 3:
    print("Usage: ", sys.argv[0], " start-date duration")
    quit()

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

t0 = spice.str2et(sys.argv[1])
days = float(sys.argv[2])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)

mSun = data[1][0]/G

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

sim.t = t0
sim.dt = 86400.0 * 0.05

def distance(body1, body2):
    dx = body1.x - body2.x
    dy = body1.y - body2.y
    dz = body1.z - body2.z
    return sqrt(dx*dx+dy*dy+dz*dz)

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    print('%14.2f %10.2f' % (sim.t-t0,sim.dt))
    kBody = 1
    for planet in planets:
        cx = sim.particles[kBody].x - sim.particles[0].x
        cy = sim.particles[kBody].y - sim.particles[0].y
        cz = sim.particles[kBody].z - sim.particles[0].z
        state=spice.spkezr(planet, sim.t, 'J2000', 'NONE', 'SUN')
        pv = state[0]
        vr = spice.vhat(pv[0:3])
        vz = spice.vhat(spice.vcrss(pv[0:3], pv[3:6]))
        vt = spice.vhat(spice.vcrss(vz, vr))
        dx = cx - pv[0]
        dy = cy - pv[1]
        dz = cz - pv[2]
        dp = np.array([dx, dy, dz], dtype = np.float)
        print('%2d %14.3f %14.3f %14.3f' % (kBody, spice.vdot(vr, dp), spice.vdot(vt, dp), spice.vdot(vz, dp)))
        kBody = kBody + 1

sim.heartbeat=heartbeat

sim.integrate(86400.0 * days)
