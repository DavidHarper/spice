#!/usr/bin/env python3

from math import sqrt
import spiceypy as spice

spice.furnsh('kernels/000331R_SK_LP0_V1P32.bsp')
spice.furnsh('kernels/naif0012.tls')
spice.furnsh('kernels/gm_de431.tpc')

t0 = spice.str2et('1997 OCT 15 09:26:09')

for i in range(0, 96):
    t = t0 + i*900
    stateE = spice.spkezr('CASSINI', t, 'J2000', 'NONE', 'EARTH')
    stateS = spice.spkezr('CASSINI', t, 'J2000', 'NONE', 'SUN')
    tstr = spice.timout(t, 'YYYY Mon DD HR:MN:SC ::UTC')
    print("%20s %16.1f %7.3f %16.1f %7.3f"% (tstr, spice.vnorm(stateE[0][0:3]), spice.vnorm(stateE[0][3:6]), spice.vnorm(stateS[0][0:3]), spice.vnorm(stateS[0][3:6])))
