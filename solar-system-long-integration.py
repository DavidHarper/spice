#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, cos, sin, pi
import os
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np

if len(sys.argv) < 3:
    print("Usage: ", sys.argv[0], " start-date duration [output-step]")
    quit()

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

# Constants taken from Souami, D. and Souchay, J. A&A 543, A133 (2012)
# "The solar system's invariable plane"
# DOI: 10.1051/0004-6361/201219011
#
# Table 4, DE405/DE406 on ICRF
# Node = 3.85262753 degrees
# Inclination = 23.00888397 degrees

ssiNode = pi * 3.85262753/180.0
ssiIncl = pi * 23.00888397/180.0

cosNode = cos(ssiNode)
sinNode = sin(ssiNode)

cosIncl = cos(ssiIncl)
sinIncl = sin(ssiIncl)

# These are the basis vectors for the invariable plane of the Solar System
# referred to the ICRF.
ssiX = np.array([cosNode, sinNode, 0.0], dtype=float)
ssiY = np.array([-sinNode * cosIncl, cosNode * cosIncl, sinIncl], dtype=float)
ssiZ = np.array([sinNode * sinIncl, -cosNode * sinIncl, cosIncl], dtype=float)

t0 = spice.str2et(sys.argv[1])
days = float(sys.argv[2])
step = 0

if len(sys.argv) >3:
    step = float(sys.argv[3])

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

sim.t = t0
sim.dt = 86400.0 * 0.05

def distance(body1, body2):
    dx = body1.x - body2.x
    dy = body1.y - body2.y
    dz = body1.z - body2.z
    return sqrt(dx*dx+dy*dy+dz*dz)

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    output_orbital_elements(sim)

def output_state_vectors(sim):
    print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    kBody = 1
    djd = sim.t/86400.0
    for planet in planets:
        cx = sim.particles[kBody].x - sim.particles[0].x
        cy = sim.particles[kBody].y - sim.particles[0].y
        cz = sim.particles[kBody].z - sim.particles[0].z
        vcx = sim.particles[kBody].vx - sim.particles[0].vx
        vcy = sim.particles[kBody].vy - sim.particles[0].vy
        vcz = sim.particles[kBody].vz - sim.particles[0].vz
        print('%14.6f %2d %18.3f %18.3f %18.3f %14.3f %14.3f %14.3f' % (djd, kBody, cx, cy, cz, vcx, vcy, vcz))
        kBody = kBody + 1
    sys.stdout.flush()

def output_orbital_elements(sim):
    print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    kBody = 1
    djd = sim.t/86400.0
    for planet in planets:
        cx = sim.particles[kBody].x - sim.particles[0].x
        cy = sim.particles[kBody].y - sim.particles[0].y
        cz = sim.particles[kBody].z - sim.particles[0].z
        vx = sim.particles[kBody].vx - sim.particles[0].vx
        vy = sim.particles[kBody].vy - sim.particles[0].vy
        vz = sim.particles[kBody].vz - sim.particles[0].vz
        dp = np.array([cx, cy, cz], dtype=float)
        dv = np.array([vx, vy, vz], dtype =float)
        pv = np.array([spice.vdot(dp, ssiX), spice.vdot(dp, ssiY), spice.vdot(dp, ssiZ), spice.vdot(dv, ssiX), spice.vdot(dv, ssiY), spice.vdot(dv, ssiZ)], dtype=float)
        elts = spice.oscelt(pv, sim.t, gmSun)
        print('%14.6f %2d %16.3f %16.12f %16.12f %16.12f %16.12f %16.12f' % (djd, kBody, elts[0], elts[1], elts[2], elts[3], elts[4], elts[5]))
        kBody = kBody + 1
    sys.stdout.flush()

if os.getenv('USE_WHFAST') != None:
    sim.integrator = "whfast"
    sim.ri_whfast.corrector = 17
    sim.ri_whfast.safe_mode = 0
    sim.dt = 86400.0
else:
    sim.ri_ias15.min_dt = 60000.0

if step == 0:
    sim.heartbeat=heartbeat
    sim.integrate(t0 + 86400.0 * days)
else:
    tNext = 0
    while True:
        tNext += step
        if tNext > days:
            break
        sim.integrate(t0 + 86400.0 * tNext)
        output_orbital_elements(sim)
