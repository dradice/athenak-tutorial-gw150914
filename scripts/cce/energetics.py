#!/usr/bin/env python3

import numpy as np
from scri import create_abd_from_h5
from scri.flux import energy_flux, angular_momentum_flux, momentum_flux
from scipy.integrate import cumulative_trapezoid

# Read the data
abd = create_abd_from_h5(
    file_name="CharacteristicExtractReduction.h5",
    file_format="spectrecce_v1",
    t_interpolate=np.linspace(0, 2500, 16384),
)
t = np.array(abd.h.t)

# Fluxes at infinity
E_dot = energy_flux(abd.h)
J_dot = angular_momentum_flux(abd.h)
P_dot = momentum_flux(abd.h)

# Integrated values
E = cumulative_trapezoid(E_dot, x=t)
J = cumulative_trapezoid(J_dot, x=t, axis=0)
P = cumulative_trapezoid(P_dot, x=t, axis=0)

# Output data
with open("energetics.txt", "w") as f:
    f.write("# 1:time 2:E 3:Jx 4:Jy 5:Jz 6:Px 7:Py 8:Pz\n")
    for i in range(t.shape[0] - 1):
        f.write(
            f"{0.5*(t[i] + t[i+1])} {E[i]} {J[i,0]} {J[i,1]} {J[i,2]} {P[i,0]} {P[i,1]} {P[i,2]}\n"
        )
