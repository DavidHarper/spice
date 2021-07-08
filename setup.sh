#!/bin/bash

sudo apt-get update

sudo apt-get install -y python3-pip python3-venv

mkdir -p ~/venv

VENV=~/venv/spice

python3 -m venv ${VENV}

source ${VENV}/bin/activate

for PKG in wheel numpy rebound spiceypy
do
    pip3 install $PKG
done

mkdir kernels

cd kernels

for URL in \
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls \
    https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/data/spk/co_1997288_97318_launch_v1.bsp \
    https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/data/spk/co_1997319_99311_i_cru_v1.bsp \
    https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/data/spk/co_1999312_01066_o_cru_v1.bsp \
    https://naif.jpl.nasa.gov/pub/naif/pds/data/co-s_j_e_v-spice-6-v1.0/cosp_1000/data/spk/041014R_SCPSE_01066_04199.bsp \
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp \
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc \
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc
do
    echo "Downloading ${URL}"
    curl -L -O "${URL}"
done

echo -e "Remember to run\n\nsource ${VENV}/bin/activate\n\nto initialise the Python environment."
