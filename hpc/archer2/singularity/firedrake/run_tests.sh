#!/bin/bash

. "/home/firedrake/firedrake/bin/activate"
cd "/home/firedrake/firedrake/src/firedrake/" || exit 2
python -m pytest -v tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"