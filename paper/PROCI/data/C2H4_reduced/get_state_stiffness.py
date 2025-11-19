#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

import chemgen_conservative as cgc
import chemgen_temperature as cgt
from parameter_space_numerics import *

def stiffness(J):
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size < 2 else np.max(nz) / np.min(nz)

pressures_atm = [1, 5, 10]
phis = [0.8, 1.0, 1.2]

for i, phi in enumerate(phis):
    for j, p_atm in enumerate(pressures_atm):
        y = np.load(f"phi_{phi}_p{p_atm}.npy")
        stiffness_ratio = np.zeros(np.shape(y[:, :, 0]))
        linear_iterations_ratio = np.zeros(np.shape(y[:, :, 0]))
        nonlinear_iterations_ratio = np.zeros(np.shape(y[:, :, 0]))
        data_shape =  np.shape(stiffness_ratio)
        for k in range(data_shape[0]):
            for l in range(data_shape[1]):
                ykl = y[k,l,:]
                T = ykl[0]
                C = ykl[1:-1]
                dt = ykl[-1]
                Jt = cgt.source_jacobian(C, T)
                Jc = cgc.source_jacobian(C, T)
                stiffness_ratio[k, l] = stiffness(Jt)/stiffness(Jc)
                _, linear_iterations_cons, running_newton_cons = backwards_euler(C, T, 1e-8, cgc, "cons", gmres_tolerance=1e-6, max_iter=10)
                _, linear_iterations_temp, running_newton_temp = backwards_euler(C, T, 1e-8, cgt, "temp", gmres_tolerance=1e-6, max_iter=10)
                linear_iterations_ratio[k, l] = linear_iterations_temp/linear_iterations_cons
                nonlinear_iterations_ratio[k, l] = running_newton_temp/running_newton_cons

        print(f"saving for phi={phi} and p={p_atm}")
        np.save(f"phi_{phi}_p{p_atm}_stiffness.npy", (stiffness_ratio))
        np.save(f"phi_{phi}_p{p_atm}_linear.npy", (linear_iterations_ratio))
        np.save(f"phi_{phi}_p{p_atm}_nonlinear.npy", (nonlinear_iterations_ratio))

