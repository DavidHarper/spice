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
from numpy.random import default_rng
import matplotlib.pyplot as plt

if len(sys.argv) < 6:
    print("Usage: ", sys.argv[0], " start-date stepsize steps swarm-size dv")
    quit()

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

t0 = spice.str2et(sys.argv[1])
dt = float(sys.argv[2]) * 86400.0
nsteps = int(sys.argv[3])
nswarm = int(sys.argv[4])
dv0 = float(sys.argv[5])

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

xSwarm = []
ySwarm = []
zSwarm = []

state = spice.spkezr('CASSINI', t0, 'J2000', 'NONE', 'SUN')
pv = state[0]
cassini = rebound.Particle(m=0, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
sim.add(cassini)

xSwarm.append(np.empty(nsteps + 1))
ySwarm.append(np.empty(nsteps + 1))
zSwarm.append(np.empty(nsteps + 1))

rng = default_rng()

for k in range(nswarm):
    ps = state[0][0:3]
    vs = state[0][3:6] + dv0 * rng.standard_normal(3)
    swarmer = rebound.Particle(m=0, x=ps[0], y=ps[1], z=ps[2], vx=vs[0], vy=vs[1], vz=vs[2])
    sim.add(swarmer)
    xSwarm.append(np.empty(nsteps + 1))
    ySwarm.append(np.empty(nsteps + 1))
    zSwarm.append(np.empty(nsteps + 1))

sim.t = t0

def report(sim, kStep):
    print('# %14.2f %10.2f' % (sim.t-t0,sim.dt))
    kBody = 1
    djd = sim.t/86400.0
    state=spice.spkezr('CASSINI', sim.t, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    vr = spice.vhat(pv[0:3])
    vz = spice.vhat(spice.vcrss(pv[0:3], pv[3:6]))
    vt = spice.vhat(spice.vcrss(vz, vr))

    for kBody in range(sim.N_active, sim.N):
        cx = sim.particles[kBody].x - sim.particles[0].x
        cy = sim.particles[kBody].y - sim.particles[0].y
        cz = sim.particles[kBody].z - sim.particles[0].z
        cvx = sim.particles[kBody].vx - sim.particles[0].vx
        cvy = sim.particles[kBody].vy - sim.particles[0].vy
        cvz = sim.particles[kBody].vz - sim.particles[0].vz
        dx = cx - pv[0]
        dy = cy - pv[1]
        dz = cz - pv[2]
        dvx = cvx - pv[3]
        dvy = cvy - pv[4]
        dvz = cvz - pv[5]
        dp = np.array([dx, dy, dz], dtype = float)
        dv = np.array([dvx, dvy, dvz], dtype = float)
        xSwarm[kBody-sim.N_active][kStep] = dx
        ySwarm[kBody-sim.N_active][kStep] = dy
        zSwarm[kBody-sim.N_active][kStep] = dz
        #xSwarm[kBody-sim.N_active][kStep] = spice.vdot(vr, dp)
        #ySwarm[kBody-sim.N_active][kStep] = spice.vdot(vt, dp)
        #zSwarm[kBody-sim.N_active][kStep] = spice.vdot(vz, dp)

sim.exact_finish_time = 1

for kStep in range(nsteps+1):
    report(sim, kStep)
    t_next = sim.t + dt
    sim.integrate(t_next)

ax = plt.figure().add_subplot(projection='3d')

for kBody in range(nswarm+1):
    ax.plot(xSwarm[kBody], ySwarm[kBody], zSwarm[kBody], lw=0.5)

ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Cassini Swarm")

plt.show()
