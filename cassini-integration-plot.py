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
import matplotlib.pyplot as plt
import re

def parse_step_size(step_str):
    p=re.compile(r"(\d+\.?\d?)([a-zA-z]?)")
    g=p.match(step_str).groups()
    stepsize=float(g[0])
    units=g[1].lower()
    if units=='' or units=='d':
        stepsize=stepsize*86400.0
    elif units=='h':
        stepsize=stepsize*3600.0
    elif units=='m':
        stepsize=stepsize*60.0
    elif units!='s':
        print("ERROR: unknown units: '" + units + "' [use d for days, h for hours, m for minutes, s for seconds]")

    return stepsize

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], " start-date stepsize steps")
    quit()

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

t0 = spice.str2et(sys.argv[1])
dt = parse_step_size(sys.argv[2])
nsteps = int(sys.argv[3])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

data = spice.bodvrd('SUN','GM',1)

mSun = data[1][0]/G

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

targets = ['MERCURY BARYCENTER', 'VENUS BARYCENTER', 'EARTH BARYCENTER',
'MARS BARYCENTER', 'JUPITER BARYCENTER', 'SATURN BARYCENTER',
'URANUS BARYCENTER', 'NEPTUNE BARYCENTER']

for planet in targets:
    data = spice.bodvrd(planet, 'GM', 1)
    mPlanet = data[1][0]/G
    state=spice.spkezr(planet, t0, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    planet = rebound.Particle(m=mPlanet, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
    sim.add(planet)

sim.N_active = sim.N

state = spice.spkezr('CASSINI', t0, 'J2000', 'NONE', 'SUN')
pv = state[0]
cassini = rebound.Particle(m=0, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
sim.add(cassini)
targets.append('CASSINI')
kCassini=sim.N-1

sim.t = t0

ts=np.empty(nsteps + 1)
dxs=np.empty(nsteps + 1)
dys=np.empty(nsteps + 1)
dzs=np.empty(nsteps + 1)
dvxs=np.empty(nsteps + 1)
dvys=np.empty(nsteps + 1)
dvzs=np.empty(nsteps + 1)

def report(sim, kStep):
    kBody = kCassini
    djd = sim.t/86400.0
    cx = sim.particles[kBody].x - sim.particles[0].x
    cy = sim.particles[kBody].y - sim.particles[0].y
    cz = sim.particles[kBody].z - sim.particles[0].z
    cvx = sim.particles[kBody].vx - sim.particles[0].vx
    cvy = sim.particles[kBody].vy - sim.particles[0].vy
    cvz = sim.particles[kBody].vz - sim.particles[0].vz
    state=spice.spkezr('CASSINI', sim.t, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    ts[kStep] = (sim.t - t0)/86400.0
    dxs[kStep] = cx - pv[0]
    dys[kStep] = cy - pv[1]
    dzs[kStep] = cz - pv[2]
    dvxs[kStep] = cvx - pv[3]
    dvys[kStep] = cvy - pv[4]
    dvzs[kStep] = cvz - pv[5]

sim.exact_finish_time = 1

report(sim, 0)

for kStep in range(1, nsteps+1):
    t_next = sim.t + dt
    sim.integrate(t_next)
    report(sim, kStep)

fig, axs = plt.subplots(2, 1, sharex=True)
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0)

axs[0].plot(ts, dxs, label='x')
axs[0].plot(ts, dys, label='y')
axs[0].plot(ts, dzs, label='z')

axs[0].legend()

axs[1].plot(ts, dvxs, label='vx')
axs[1].plot(ts, dvys, label='vy')
axs[1].plot(ts, dvzs, label='vz')

axs[1].legend()

plt.show()
