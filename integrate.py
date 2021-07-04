#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/

import rebound
from math import sqrt
import sys

if len(sys.argv) < 2:
    print("Usage: ", sys.argv[0], " days")
    quit()

days = float(sys.argv[1])

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)

file = open('initialstate.out', 'r')

t0 = float(file.readline())

mSun = float(file.readline().replace('D','e'))/G

nMass = int(file.readline())

masses = []

for i in range(0, nMass):
    masses.append(float(file.readline().replace('D','e'))/G)

nBody = int(file.readline())

state = []

for i in range(0, nBody):
    str = file.readline().replace('D','e')
    state.append(list(map(float, str.split())))

file.close()

sim = rebound.Simulation()

sim.G = G

sun = rebound.Particle(m=mSun)

sim.add(sun)

for i in range(0, len(masses)):
    pv = state[i]
    planet = rebound.Particle(m=masses[i], x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
    sim.add(planet)

for i in range(len(masses), len(state)):
    pv = state[i]
    target = rebound.Particle(x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
    sim.add(target)

sim.t = 0
sim.dt = 86400.0 * 0.05

def distance(body1, body2):
    dx = body1.x - body2.x
    dy = body1.y - body2.y
    dz = body1.z - body2.z
    return sqrt(dx*dx+dy*dy+dz*dz)

def heartbeat(sim_pointer):
    sim = sim_pointer.contents
    t=t0+sim.t/86400.0
    for i in range(1, len(sim.particles)):
        x = sim.particles[i].x - sim.particles[0].x
        y = sim.particles[i].y - sim.particles[0].y
        z = sim.particles[i].z - sim.particles[0].z
        vx = sim.particles[i].vx - sim.particles[0].vx
        vy = sim.particles[i].vy - sim.particles[0].vy
        vz = sim.particles[i].vz - sim.particles[0].vz
        print('%14.6f %2d %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e' % (
        t, i, x, y, z, vx, vy, vz
        ))

sim.heartbeat=heartbeat

sim.integrate(86400.0 * days)
