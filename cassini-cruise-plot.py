#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

from math import sqrt, log10
import sys
import glob
import spiceypy as spice
import rebound
import numpy as np
import matplotlib.pyplot as plt
import re

def cassini_offset(sim):
    kCassini=sim.N-1
    cx = sim.particles[kCassini].x - sim.particles[0].x
    cy = sim.particles[kCassini].y - sim.particles[0].y
    cz = sim.particles[kCassini].z - sim.particles[0].z
    #cvx = sim.particles[kCassini].vx - sim.particles[0].vx
    #cvy = sim.particles[kCassini].vy - sim.particles[0].vy
    #cvz = sim.particles[kCassini].vz - sim.particles[0].vz
    state=spice.spkezr('CASSINI', sim.t, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    return [cx - pv[0], cy - pv[1], cz - pv[2]]

def seconds_to_year(ts):
    return 2000.0+ts/(86400.0*365.25)

# The following list of TCMs is from Table 3 of the paper
# "Cassini-Huygens Maneuver Experience: Cruise and Arrival at Saturn"
# by T.D. Goodson, B.B. Buffington, Y. Hahn, N.J. Strange, S.V. Wagner, M. Wong
# in Proceedings of the AAS/AIAA Astrodynamics Specialist Conference,
# 7-11 August 2005
#
# Flyby dates are from Table 3 of the paper
# "Cassini Maneuver Experience: Finishing the Inner Cruise"
# by T.D. Goodson, D.L.Gray, Y. Hahn, F. Peralta
#
# Start date is launch plus one day, finish date is Saturn arrival minus one day.
#
# Note that there are gaps in the TCM numbering sequence because some TCMs
# were cancelled.  DSM is the Deep Space Manoeuvre.
events = [
    ['Start',   '16 Oct 1997 12:00'],
    ['TCM-1',    '9 Nov 1997 20:00'],
    ['TCM-2',   '25 Feb 1998 20:00'],
    ['Venus1',  '26 Apr 1998 13:45'],
    ['DSM',      '3 Dec 1998 06:00'],
    ['TCM-6',    '4 Feb 1999 20:00'],
    ['TCM-7',   '18 May 1999 17:00'],
    ['Venus2',  '24 Jun 1999 18:25'],
    ['TCM-9',    '6 Jul 1999 16:00'],
    ['TCM-10',  '19 Jul 1999 16:00'],
    ['TCM-11',   '2 Aug 1999 21:30'],
    ['TCM-12',  '11 Aug 1999 15:30'],
    ['Earth',   '18 Aug 1999 03:28'],
    ['TCM-13',  '31 Aug 1999 16:00'],
    ['TCM-14',  '14 Jun 2000 17:00'],
    ['Jupiter', '30 Dec 2000 10:36'],
    ['TCM-17',  '28 Feb 2001 17:30'],
    ['TCM-18',   '3 Apr 2002 18:00'],
    ['TCM-19',   '1 May 2003 20:00'],
    ['TCM-19a', '10 Sep 2003 20:00'],
    ['TCM-19b',  '2 Oct 2003 02:00'],
    ['TCM-20',  '27 May 2004 22:26'],
    ['TCM-21',  '16 Jun 2004 21:07'],
    ['Finish',  '30 Jun 2004 12:00']
]

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

t_data=[]
dx_data=[]
dy_data=[]
dz_data=[]

delta_t=6*3600.0

G = 6.6743e-20 # km^3 kg^(-1) s^(-2)
data = spice.bodvrd('SUN','GM',1)
mSun = data[1][0]/G

planets = ['MERCURY BARYCENTER', 'VENUS BARYCENTER', 'EARTH BARYCENTER',
    'MARS BARYCENTER', 'JUPITER BARYCENTER', 'SATURN BARYCENTER',
    'URANUS BARYCENTER', 'NEPTUNE BARYCENTER']

dt=3600.0

for i in range(0, len(events)-1):
    print("Segment " + str(i) + " : from " + events[i][0] + " to " + events[i+1][0])
    t_start = spice.str2et(events[i][1])+delta_t
    t_end_event = spice.str2et(events[i+1][1])
    t_finish = t_end_event-delta_t

    sim = rebound.Simulation()

    sim.G = G

    sun = rebound.Particle(m=mSun)

    sim.add(sun)

    for planet in planets:
        data = spice.bodvrd(planet, 'GM', 1)
        mPlanet = data[1][0]/G
        state=spice.spkezr(planet, t_start, 'J2000', 'NONE', 'SUN')
        pv = state[0]
        planet = rebound.Particle(m=mPlanet, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
        sim.add(planet)

        sim.N_active = sim.N

    state = spice.spkezr('CASSINI', t_start, 'J2000', 'NONE', 'SUN')
    pv = state[0]
    cassini = rebound.Particle(m=0, x=pv[0], y=pv[1], z=pv[2], vx=pv[3], vy=pv[4], vz=pv[5])
    sim.add(cassini)

    sim.t = t_start
    sim.exact_finish_time = 1

    # Add a NaN to the data to force Matplotlib to draw disjoint line segments.
    dp=cassini_offset(sim)
    t_data.append(seconds_to_year(sim.t))
    dx_data.append(dp[0])
    dy_data.append(dp[1])
    dz_data.append(dp[2])

    last_step=False

    while not last_step:
        t_next = sim.t + dt
        if t_next > t_finish:
            t_next=t_finish
            last_step=True
        sim.integrate(t_next)
        dp=cassini_offset(sim)
        t_data.append(seconds_to_year(sim.t))
        dx_data.append(dp[0])
        dy_data.append(dp[1])
        dz_data.append(dp[2])

    t_data.append(seconds_to_year(t_end_event))
    dx_data.append(np.nan)
    dy_data.append(np.nan)
    dz_data.append(np.nan)

fig, ax = plt.subplots()

ax.plot(t_data, dx_data, label='x')

ax.plot(t_data, dy_data, label='y')

ax.plot(t_data, dz_data, label='z')

y_label_values=[2000.0, 2400.0, 2800.0, 3200.0, 3600.0, 4000.0, 4400.0, 4800.0]
kLabel=0

for event in events:
    label=event[0]
    if label!='Start' and label!='Finish':
        t_label = seconds_to_year(spice.str2et(event[1]))
        y_label=y_label_values[kLabel%len(y_label_values)]
        ax.annotate(label,
            xy=(t_label, y_label),  # theta, radius
            xytext=(t_label, y_label),    # fraction, fraction
            horizontalalignment='center'
            )
        kLabel=kLabel+1


ax.legend()

plt.show()
