#!/usr/bin/env python3

import spiceypy as spice
import sys
from pathlib import Path
import glob

if len(sys.argv) < 2:
    print("Usage: ", sys.argv[0], " start-date")
    quit()

kernels=str(Path.home()) + '/kernels/*.*'

for kernel in glob.glob(kernels):
    spice.furnsh(kernel)

t=spice.str2et(sys.argv[1])

for target in ['SUN','VENUS','EARTH','JUPITER','SATURN']:
    p=spice.bodvrd(target, 'GM', 1)
    gm=p[1][0]
    print(target,' ',gm)

for target in ['VENUS BARYCENTER','EARTH BARYCENTER','JUPITER BARYCENTER','SATURN BARYCENTER','CASSINI']:
    state=spice.spkezr(target,t,'J2000','NONE','SUN')
    print(target,' ',state[0])
