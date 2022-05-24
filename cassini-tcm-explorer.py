#!/usr/bin/env python3
#
# This program uses the REBOUND integrator by Hanno Rein and colleagues.
# https://rebound.readthedocs.io/en/latest/
#
# It also uses the SpiceyPy wrapper for NAIF SPICE

import sys
import glob
import spiceypy as spice
import numpy as np
from math import sqrt

for file in glob.glob('kernels/*.*'):
    spice.furnsh(file)

def cassini_velocity(t):
    state=spice.spkezr('CASSINI', t, 'J2000', 'NONE', 'SUN')
    return np.array(state[0][3:6], dtype=float)

def analyse_tcm(tcm_name, tcm_time, dt):
    print("Analysing " + tcm_name)
    t0 = spice.str2et(tcm_time)
    print("t0 = " + tcm_time + " = " + str(t0))

    vm2=cassini_velocity(t0-2.0*dt)
    vm1=cassini_velocity(t0-dt)
    vp1=cassini_velocity(t0+dt)
    vp2=cassini_velocity(t0+2.0*dt)

    print("vm2 = " + str(vm2))
    print("vm1 = " + str(vm1))
    print("vp1 = " + str(vp1))
    print("vp2 = " + str(vp2))

    am=(vm1-vm2)/dt
    ap=(vp2-vp1)/dt
    a0=0.5*(am+ap)

    print("am = " + str(am))
    print("ap = " + str(ap))
    print("a0 = " + str(a0))

    deltav=1000.0*(vp1-vm1-2.0*dt*a0)

    dv=sqrt(deltav[0]*deltav[0]+deltav[1]*deltav[1]+deltav[2]*deltav[2])

    print("delta v = " + str(deltav))
    print("dv = " + str(dv))
    print()

# The following list of TCMs is from Table 3 of the paper
# "Cassini-Huygens Maneuver Experience: Cruise and Arrival at Saturn"
# by T.D. Goodson, B.B. Buffington, Y. Hahn, N.J. Strange, S.V. Wagner, M. Wong
# in Proceedings of the AAS/AIAA Astrodynamics Specialist Conference,
# 7-11 August 2005
#
# Note that there are gaps in the TCM numbering sequence because some TCMs
# were cancelled.  DSM is the Deep Space Manoeuvre.
tcms = [
    ['TCM-1',    '9 Nov 1997 20:00'],
    ['TCM-2',   '25 Feb 1998 20:00'],
    ['DSM',      '3 Dec 1998 06:00'],
    ['TCM-6',    '4 Feb 1999 20:00'],
    ['TCM-7',   '18 May 1999 17:00'],
    ['TCM-9',    '6 Jul 1999 16:00'],
    ['TCM-10',  '19 Jul 1999 16:00'],
    ['TCM-11',   '2 Aug 1999 21:30'],
    ['TCM-12',  '11 Aug 1999 15:30'],
    ['TCM-13',  '31 Aug 1999 16:00'],
    ['TCM-14',  '14 Jun 2000 17:00'],
    ['TCM-17',  '28 Feb 2001 17:30'],
    ['TCM-18',   '3 Apr 2002 18:00'],
    ['TCM-19',   '1 May 2003 20:00'],
    ['TCM-19a', '10 Sep 2003 20:00'],
    ['TCM-19b',  '2 Oct 2003 02:00'],
    ['TCM-20',  '27 May 2004 22:26'],
    ['TCM-21',  '16 Jun 2004 21:07']
]

for tcm in tcms:
    analyse_tcm(tcm[0], tcm[1], 7200.0)
