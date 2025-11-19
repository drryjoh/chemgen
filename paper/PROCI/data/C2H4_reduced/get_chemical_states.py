#!python3
#!/usr/bin/env python3
import numpy as np
import cantera as ct
from copy import deepcopy

# ----------------------------
# User controls / constants
# ----------------------------
MECH = "mechanism.yaml"         # mechanism file
FUEL = "C2H4:1"             # ethylene as fuel
OX   = "O2:1, N2:3.76"      # air
pressures_atm = [1.0, 5.0, 10.0]
phis = [0.8, 1.0, 1.2]
T_init_grid = np.linspace(1800.0, 2800.0, 20)  # 10 temperatures
T_init_grid = T_init_grid[::-1]
dt = 1e-9                  # s, small fixed step for induction-time search
tmax = 0.1                 # s, safety cap for induction-time search
n_time_samples = 100        # 10 points along time: [0, 5*tau_id] inclusive

# ----------------------------
# Helpers
# ----------------------------
def set_initial_state(gas, T0, p_atm, phi):
    gas.TP = T0, p_atm * ct.one_atm
    gas.set_equivalence_ratio(phi, fuel=FUEL, oxidizer=OX)

def induction_time_via_dTdt_peak(gas, dt, tmax):
    """
    Constant-volume homogeneous reactor, advance with fixed target times
    at dt increments. Return time at which dT/dt attains its maximum.
    """
    gas0 = ct.Solution(gas.source)  # fresh Solution pointing to same mech
    gas0.TPX = gas.T, gas.P, gas.X  # store initial TPX

    r = ct.IdealGasReactor(gas)
    net = ct.ReactorNet([r])

    t = 0.0
    T_prev = r.T
    t_prev = 0.0

    dTdt_max = -np.inf
    t_at_max = 0.0

    # Stop early once temperature jumps significantly (typical ignition criterion)
    T_jump_threshold = 400.0

    while t < tmax:
        t += dt
        net.advance(t)

        # finite difference for dT/dt
        dTdt = (r.T - T_prev) / (t - t_prev) if t > t_prev else 0.0
        if dTdt > dTdt_max:
            dTdt_max = dTdt
            t_at_max = t

        if r.T - gas0.T > T_jump_threshold:
            # past strong exothermic rise; assume peak found
            break

        T_prev = r.T
        t_prev = t

    # restore gas to initial state before returning
    gas.TPX = gas0.T, gas0.P, gas0.X
    return max(t_at_max, dt)  # avoid zero

def sample_trajectory(gas, total_time, n_samples):
    """
    Re-run a constant-volume homogeneous reactor from the same initial
    condition and record [T, concentrations...] at n_samples uniform times
    up to total_time.
    """
    # capture initial state to ensure fresh run
    T0, P0, X0 = gas.T, gas.P, gas.X

    r = ct.IdealGasReactor(gas)
    net = ct.ReactorNet([r])

    times = np.linspace(0.0, total_time, n_samples)
    out = []

    # Step to each target time
    t_prev = net.time
    for targ in times:
        net.advance(targ)
        dt = net.time - t_prev
        vec = np.concatenate(([r.T], r.thermo.concentrations, np.array([dt])))
        t_prev = net.time
        out.append(vec)


    # restore gas
    gas.TPX = T0, P0, X0
    return np.vstack(out)  # shape (n_samples, 1 + n_species)

# ----------------------------
# Main
# ----------------------------
def main():
    gas = ct.Solution(MECH)
    n_species = gas.n_species

    # For each (phi, p): build array of shape (10 T0s, 10 times, 1+n_species)
    for phi in phis:
        for p_atm in pressures_atm:
            # container for 10×10 states
            t_by_T = np.zeros((len(T_init_grid), n_time_samples, 1 + n_species + 1)) #now includes time

            for iT, T0 in enumerate(T_init_grid):
                # Set initial state for this case
                set_initial_state(gas, T0, p_atm, phi)

                # 1) find ignition delay via dT/dt peak using small dt
                tau_id = induction_time_via_dTdt_peak(gas, dt, tmax)

                # 2) sample again from initial state up to 5*tau_id at 10 points
                total_time = 5.0 * tau_id
                # Problem statement: “calculate a new homogeneous reactor with
                # the initial condition at dt = 5*ind/10” which is equivalent to
                # uniform sampling of 10 points to total_time.
                set_initial_state(gas, T0, p_atm, phi)
                traj = sample_trajectory(gas, total_time, n_time_samples)  # (10, 1+n_species)

                t_by_T[iT, :, :] = traj

            # Save one file per (phi, p)
            fname = f"phi_{phi}_p{int(p_atm)}.npy"
            np.save(fname, t_by_T)
            # Optional: save metadata with species names
            with open(f"phi_{phi}_p{int(p_atm)}_species.txt", "w") as f:
                f.write("Order: [T, " + ", ".join(gas.species_names) + "]\n")
            print(f"Saved {fname} with shape {t_by_T.shape} (T0 x time x state)")

if __name__ == "__main__":
    main()

