#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, pi, sin, cos, atan2, fabs
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np

if len(sys.argv) < 2:
    print("Usage: ", sys.argv[0], " t0 days step")
    quit()

spice.furnsh('kernels/gm_de431.tpc')
spice.furnsh('kernels/de430.bsp')
spice.furnsh('kernels/naif0012.tls')

t0 = spice.str2et(sys.argv[1])
days = float(sys.argv[2])
step = float(sys.argv[3])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)
gmSun = data[1][0]
mSun = gmSun/G

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

planets = ['EARTH', 'MOON']

for planet in planets:
    data = spice.bodvrd(planet, 'GM', 1)
    mPlanet = data[1][0]/G
    state=spice.spkezr(planet, t0, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    planet = rebound.Particle(m=mPlanet, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
    sim.add(planet)

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    tjd = sim.t/86400.0

    sun = sim.particles[0]
    earth = sim.particles[1]
    moon = sim.particles[2]

    pMoon = np.array([moon.x - earth.x, moon.y - earth.y, moon.z - earth.z], dtype = float)
    vMoon = np.array([moon.vx - earth.vx, moon.vy - earth.vy, moon.vz - earth.vz], dtype = float)

    u = spice.vhat(pMoon)
    w = spice.vhat(spice.vcrss(pMoon, vMoon))
    v = spice.vhat(spice.vcrss(w, u))

    state = spice.spkezr('MOON', sim.t, 'J2000', 'NONE', 'EARTH')

    pvMoon = state[0]

    pMoon0 = np.array([pvMoon[0], pvMoon[1], pvMoon[2]])

    dp = spice.vsub(pMoon, pMoon0)

    du = spice.vdot(u, dp)
    dv = spice.vdot(v, dp)
    dw = spice.vdot(w, dp)

    print('%14.6f %14.3f %14.3f %14.3f' % (tjd, du, dv, dw))

sim.heartbeat=heartbeat
sim.t = 0.0
sim.dt = 86400.0 * step

sim.integrate(86400.0 * days)
