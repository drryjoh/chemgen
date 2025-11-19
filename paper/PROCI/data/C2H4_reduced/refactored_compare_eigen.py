#!python3
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# Existing chemgen backends
import chemgen_conservative as cgc
import chemgen_temperature as cgt
import chemgen_dtemperature_ignore as cgc_it
import chemgen_conservative_ignore_species as cgc_is
import chemgen_temperature_ignore_species as cgt_is
# ---------- Numerics-------
import math
from scipy.sparse.linalg import gmres, LinearOperator
import eigen_gmres
def departure(A):
    return np.sum(~A.any(axis=1))
    #return np.linalg.norm(A @ A.T.conj() - A.T.conj() @ A)

def error_norm(x, x_ref, abs_tolerance, rel_tolerance):
    denom = abs_tolerance + rel_tolerance * np.abs(x_ref)
    return np.linalg.norm(x / denom)

def full_jacobian(J, dt):
    n = np.shape(J)
    I = np.eye(n[0])
    #gamma = 1.0 - 1.0 / math.sqrt(2.0)
    #Backward Euler
    return  J

def jacobi_precondition(A, eps=1e-14):
    """Return D^{-1}A and D^{-1} itself."""
    d = np.diag(A)
    nz = np.abs(d) > eps
    invd = np.zeros_like(d)
    invd[nz] = 1.0 / d[nz]
    M_inv = np.diag(invd)
    return M_inv @ A

def gmres_with_iterations(A, b, tol=1e-4, maxiter=100):
    diag_A = np.diag(A)
    M_inv = 1.0 / diag_A

    # Define LinearOperator for preconditioner application: z = M @ r
    def Mv(v):
        return M_inv * v

    M = LinearOperator(A.shape, matvec=Mv)

    iters = [0]
    def counter(residual):
        iters[0] += 1
    x, info = gmres(A, b, M=M, rtol=tol, atol=1e-7, maxiter=maxiter, callback=counter)
    return x, info, iters[0]

def sdirk2(C, T, dt, be, gmres_tolerance = 1e-6, max_iter=10, gmres_method = "numpy"):
    abs_newton_tol = 1e-7
    rel_newton_tol = 1e-6
    
    mod = be["mod"] #backedn module
    y = np.zeros(len(C)+1)
    y[1:] = C
    if be["require_internal_energy"] == "yes":
        y[0] = np.array(mod.internal_energy_volume_specific(C, T))
    else:
        y[0] = T

    y_init = y
    gamma = 1.0 - 1.0 / math.sqrt(2.0)
    one_minus_gamma = 1.0 - gamma
    temperature_init = be["temp_from_state"](y) #be["temp_from_state"](y)
    k1 = np.array(mod.source(y_init[1:], temperature_init))
    I = np.eye(len(k1))

    # stage 1
    running_newton = 0
    linear_iterations = 0
    for _ in range(max_iter):
        y_stage = y_init + gamma * dt * k1
        T_ = be["temp_from_state"](y_stage)
        f_val = np.array(mod.source(y_stage[1:], T_))
        res = f_val - k1
        J = I - gamma * dt * np.array(mod.source_jacobian(y_stage[1:], T_))
        if gmres_method == "numpy":
            dk, info, linear_iteration_1 = gmres_with_iterations(J, res)
        else:
            res = eigen_gmres.gmres_dense(J, res, rtol=1e-6, maxiter=100, restart=10)
            dk = res.x
            info = res.info
            linear_iteration_1 =  res.iters
    
        linear_iterations += linear_iteration_1
        if info > 0:
            print(f"Did not converge after {info} iterations")
        k1 = k1 + dk
        running_newton+=1
        if error_norm(dk, k1, abs_newton_tol, rel_newton_tol) < 1.0:
            break

    # stage 2
    k2 = k1
    for _ in range(max_iter):
        y_stage = y_init + one_minus_gamma * dt * k1 + gamma * dt * k2
        T_ = be["temp_from_state"](y_stage)
        f_val = np.array(mod.source(y_stage[1:], T_))
        res = f_val - k2
        J = I - gamma * dt * np.array(mod.source_jacobian(y_stage[1:], T_))
        if gmres_method == "numpy":
            dk, info, linear_iteration_2 = gmres_with_iterations(J, res)
        else:
            res = eigen_gmres.gmres_dense(J, res, rtol=1e-10, maxiter=100, restart=10)
            dk = res.x
            info = res.info
            linear_iteration_2 =  res.iters
        linear_iterations += linear_iteration_2
        if info > 0:
            print(f"Did not converge after {info} iterations")
        k2 = k2 + dk
        running_newton+=1
        if error_norm(dk, k2, abs_newton_tol, rel_newton_tol) < 1.0:
            break
            
    y_next = y_init + dt * one_minus_gamma * k1 + gamma * dt * k2
    return y_next, linear_iterations, running_newton

import numpy as np
import math

def backwards_euler(C, T, dt, be, gmres_tolerance=1e-6, max_iter=10, gmres_method="numpy"):
    """
    Implicit backward-Euler step:
      Solve [(I/dt) - J(y_guess)] * dy = f(y_guess) - (1/dt) * (y_guess - y_init)
      y_{n+1} = y_guess + dy, iterate Newton until error_norm < 1.
    Returns:
      y_next, total_linear_iters, newton_iters
    """
    abs_newton_tol = 1e-7
    rel_newton_tol = 1e-6
    # Build initial state y = [energy_or_T, concentrations...]
    mod = be["mod"]
    y = np.zeros(len(C) + 1, dtype=float)
    y[1:] = C
    if be.get("require_internal_energy", "no") == "yes":
        y[0] = float(mod.internal_energy_volume_specific(C, T))
    else:
        y[0] = float(T)

    y_init = y.copy()
    y_guess = y.copy()

    nvars = y_init.size
    I = np.eye(nvars, dtype=float)

    running_newton = 0
    linear_iterations = 0

    for iter_ in range(max_iter):
        # Temperature from current state
        T_ = be["temp_from_state"](y_guess)
        f = np.array(mod.source(y_guess[1:], T_), dtype=float)  
        J = np.array(mod.source_jacobian(y_guess[1:], T_), dtype=float) 
        A = I/dt - J
        res = f - (1.0 / dt) * (y_guess - y_init)

        # Solve for dy
        if be["key"] == "cgc":
            A_reduced = A[1:, 1:]
            res_reduced = res[1:]

            # Solve reduced system
            if gmres_method == "numpy":
                dy_reduced, info, iters_lin = gmres_with_iterations(A_reduced, res_reduced)
            else:
                r = eigen_gmres.gmres_dense(A_reduced, res_reduced, rtol=gmres_tolerance, maxiter=100, restart=50)
                dy_reduced = r.x
                info = r.info
                iters_lin = r.iters

            # Add a zero at the beginning of the solution
            dy = np.insert(dy_reduced, 0, 0.0)
        else:
            if gmres_method == "numpy":
                dy, info, iters_lin = gmres_with_iterations(A, res)
            else:
                r = eigen_gmres.gmres_dense(A, res, rtol=gmres_tolerance, maxiter=100, restart=50)
                dy = r.x
                info = r.info
                iters_lin = r.iters

        linear_iterations += iters_lin
        if info > 0:
            print(f"GMRES did not converge after {info} iterations")

        # Newton update
        y_guess = y_guess + dy
        running_newton += 1

        if error_norm(dy, y_guess, abs_newton_tol, rel_newton_tol) < 1.0:
            break

    # Return final state and counters
    return y_guess, linear_iterations, running_newton


# ---------- Utilities ----------
def inv_largest_tau(J):
    # tau_min = 1 / max(|lambda|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size == 0 else 1.0 / np.max(nz)

def inv_smallest_tau(J):
    # tau_max = 1 / min(|lambda|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size == 0 else 1.0 / np.min(nz)

def eigs(J):
    # tau_max = 1 / min(|lambda|)
    eig = np.linalg.eigvals(J)
    return np.real(eig)
    '''
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.array([])
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.array(nz)
    '''
def svs(A):
    return np.linalg.svd(A, compute_uv=False)  
def stiffness(J):
    # kappa = tau_max / tau_min = max(|lambda|)/min(|lambda|)
    eig = np.linalg.eigvals(J)
    mag = np.abs(np.real(eig))
    if mag.size == 0:
        return np.nan
    tol = 1e-12 * np.max(mag)
    nz = mag[mag > tol]
    return np.nan if nz.size < 2 else np.max(nz) / np.min(nz)

# Colorblind-friendly (Okabe–Ito)
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#000000", "#CC79A7", "#F0E442"]

# ---------- Register backends ----------
# Each backend supplies:
# - key: identifier
# - label: for legends
# - mod: imported module object
# - temp_from_state(y): extract temperature from state vector y returned by backwards_euler
# - style: plotting style

BACKENDS = [
    dict(
        key="cgc",
        label="Conservative",
        require_internal_energy = "yes",
        mod=cgc,
        temp_from_state=lambda y: cgc.temperature(y[1:], y[0]),
    ),
    dict(
        key="cgt",
        label="Temperature",
        require_internal_energy = "no",
        mod=cgt,
        temp_from_state=lambda y: y[0],
    ),
    #dict(
    #    key="cgc_it",
    #    label="Conservative_nodT",
    #    mod=cgc_it,
    #    temp_from_state=lambda y: cgc_it.temperature(y[1:], y[0]),
    #),
    #dict(
    #    key="cgc_is",
    #    label="Conservative_nodS",
    #    require_internal_energy = "yes",
    #    mod=cgc_is,
    #    temp_from_state=lambda y: cgc_it.temperature(y[1:], y[0]),
    #),
    #dict(
    #    key="cgt_is",
    #    label="Temperature_nodS",
    #    require_internal_energy = "no",
    #    mod=cgt_is,
    #    temp_from_state=lambda y: y[0],
    #),
    # When you add more, copy one of these blocks and adjust:
    # dict(key="cga", label="Alt-1", mod=cga, temp_from_state=lambda y: ...),
    # dict(key="cgb", label="Alt-2", mod=cgb, temp_from_state=lambda y: ...),
]

# Assign styles deterministically
for i, be in enumerate(BACKENDS):
    be["style"] = dict(color=PALETTE[i % len(PALETTE)],
                       ls="--" if i % 2 else "-.",
                       lw = 2)

# Choose a baseline for ratio plots (defaults to first backend if not found)
BASELINE_KEY = "cgc" if any(b["key"] == "cgc" for b in BACKENDS) else BACKENDS[0]["key"]

# ---------- Cantera setup ----------
gas = ct.Solution("mechanism.yaml")
gas.set_equivalence_ratio(phi=1.0, fuel="C2H4", oxidizer="O2:1.0, N2:3.76")  # O2:N2 = 1:3.76
gas.TP = 2500, 101325
reactor = ct.IdealGasReactor(gas)
net = ct.ReactorNet([reactor])

n_steps = 200
dt = 2e-7
time = np.linspace(0.0, n_steps * dt, n_steps)

# Reference trajectory from Cantera
T_phys = []

# Per-backend storage
store = {}
for be in BACKENDS:
    store[be["key"]] = dict(
        T_pred=[],
        inv_tau_min=[],
        inv_tau_max=[],
        kappa=[],
        lin_iter=[],
        newt_iter=[],
        eigs=[],
        svs=[],
        departure=[],
    )

# ---------- Time integration loop ----------
for t in time:
    T = reactor.T
    C = reactor.thermo.concentrations

    # Run each backend at the current physical state
    for be in BACKENDS:
        mod = be["mod"]
        y, n_lin, n_newt =  mod.backwards_euler(C, T, dt) #sdirk2(C,T,dt,be, gmres_method = "numpy") backwards_euler(C,T,dt,be, gmres_method = "eigen") #
        T_hat = be["temp_from_state"](y)

        s = store[be["key"]]
        s["T_pred"].append(T_hat)
        s["lin_iter"].append(n_lin)
        s["newt_iter"].append(n_newt)

    # Advance the physical reactor to this time
    net.advance(t)
    T_phys.append(reactor.T)

    # Form Jacobians at the updated physical state
    T_now = reactor.T
    C_now = reactor.thermo.concentrations
    for be in BACKENDS:
        J = np.array(be["mod"].source_jacobian(C_now, T_now))
        s = store[be["key"]]
        s["inv_tau_min"].append(inv_largest_tau(J))
        s["inv_tau_max"].append(inv_smallest_tau(J))
        # --- CHANGED: flatten eigenvalue accumulation ---
        s["eigs"].extend(eigs(J))       # extend instead of append
        s["svs"].append(svs(jacobi_precondition(full_jacobian(J,dt)))) #jacobi_precondition(
        s["kappa"].append(stiffness(J))
        s["departure"].append(departure(J))
        # --- END CHANGED ---

# Convert to arrays
T_phys = np.array(T_phys)
for be in BACKENDS:
    s = store[be["key"]]
    for k in s:
        s[k] = np.array(s[k])

time = time * 1e6

import seaborn as sns

sns.set(style="whitegrid", context="talk")

plt.figure(figsize=(15, 10))

# collect all eigenvalues per backend
for be in BACKENDS:
    if be["label"] == "Conservative" or be["label"] == "Temperature":
        key = be["key"]
        eigs_all = np.array(store[key]["eigs"])
        eigs_all = eigs_all[eigs_all != 0]  # avoid log(0)
        
        sns.histplot(
            eigs_all,
            bins=100,
            log_scale=(True, False),   # log x-axis
            alpha=0.5,                 # semi-transparent
            element="step",            # outlined style
            fill=True,
            label=be["label"],
            kde=False
        )

plt.xlabel(r"$1/|Re(\lambda)|$  [s]")
plt.ylabel("Count")
plt.legend()
plt.tight_layout()

fig, ax = plt.subplots(2, 2, figsize=(15, 10), sharey=False)
axs = [ax[0,0],ax[0,1],ax[1,0], ax[1,1]]
conservative_linear = store["cgc"]["lin_iter"]
temperature_linear = store["cgt"]["lin_iter"]
indexes_exceed = [i for i, (a, b) in enumerate(zip(conservative_linear, temperature_linear)) if b > a]

def nn_distances(x):
    x = np.sort(np.asarray(x, dtype=float))
    d = np.empty_like(x)
    d[0]  = x[1] - x[0]
    d[-1] = x[-1] - x[-2]
    d[1:-1] = np.minimum(x[1:-1]-x[:-2], x[2:]-x[1:-1])
    return d, x

def clustering_scores(x):
    d, xs = nn_distances(x)
    mean_d = d.mean()
    cv_nn  = d.std(ddof=0) / mean_d
    norm   = (xs[-1]-xs[0]) / (len(xs)-1) if len(xs) > 1 else np.nan
    R_nn   = mean_d / norm
    return {"R_NN": R_nn, "CV_NN": cv_nn}
import matplotlib.colors as mcolors
for be in BACKENDS:
    if be["label"] == "Conservative" or be["label"] == "Temperature":
        if be["label"] == "Conservative":
            idx = 0
        else:
            idx = 2
        key = be["key"]
        svs_all = store[key]["svs"]
        for i, t in enumerate(time):
            svs_t = np.array(svs_all[i])
            st = be["style"]
            #for _ in range(8):
            #    idx_min = np.argmin(svs_t)
            #    idx_max = np.argmax(svs_t)
            #svs_t = svs_t[svs_t<10**9]#np.delete(svs_t, [idx_min])
            #svs_t = svs_t[svs_t>10**3]
            axs[idx].semilogy( svs_t*0+time[i], svs_t,'o',mfc=mcolors.to_rgba(st["color"], alpha=0.7),mec='black', markersize=4, label=be["label"])
            cluster = clustering_scores(np.log(svs_t))
            #axs[1].plot(time[i], cluster["R_NN"],'o',mfc=mcolors.to_rgba(st["color"], alpha=0.7),mec='black', markersize=4, label=be["label"])
            axs[1].plot(time[i], cluster["R_NN"],'o',mfc=mcolors.to_rgba(st["color"], alpha=0.7),mec='black', markersize=4, label=be["label"])
            axs[1].set_ylabel("R_NN")
            axs[3].plot(time[i], cluster["CV_NN"],'o',mfc=mcolors.to_rgba(st["color"], alpha=0.7),mec='black', markersize=4, label=be["label"])
            axs[3].set_ylabel("CV_NN")
            #for j in range(10):
            #    idx_min = np.argmin(svs_t)
            #    idx_max = np.argmax(svs_t)
            #    svs_t = np.delete(svs_t, [idx_min, idx_max])
            #    axs[idx+1].semilogy(time[i], np.std(svs_t),'o',mfc=mcolors.to_rgba(st["color"], alpha=0.7),mec='black', markersize=4, label=be["label"])
            if t>20:
                break

plt.show()


# ---------- Plots ----------
fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

# 1) Temperature
axs[0].plot(time, T_phys, color="#000000", label="Cantera")
for be in BACKENDS:
    st = be["style"]
    axs[0].plot(time, store[be["key"]]["T_pred"], label=f"{be['label']}", **st)
axs[0].set_ylabel("Temperature [K]")
axs[0].set_title("Temperature evolution")
axs[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 4, fontsize = 10)

# 2) Inverse eigenvalues: tau_min and tau_max per backend
for be in BACKENDS:
    st = be["style"]
    axs[1].semilogy(time, store[be["key"]]["inv_tau_min"],
                    label=rf"$\tau_{{min}}$ {be['label']}", **st)
    # derive variant style for tau_max: dotted of the same color
    axs[1].semilogy(time, store[be["key"]]["inv_tau_max"],
                    linestyle=":", color=st["color"],
                    label=rf"$\tau_{{max}}$ {be['label']}")
axs[1].set_ylabel(r"Characteristic time $\tau = 1/\lambda$")
axs[1].set_title("Inverse of largest and smallest nonzero eigenvalues")
axs[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 4, fontsize=10)

# 3) Stiffness per backend
for be in BACKENDS:
    st = be["style"]
    axs[2].semilogy(time, store[be["key"]]["kappa"],
                    label=f"{be['label']}", **st)
axs[2].set_ylabel(r"Stiffness, $\tau_{max}/\tau_{min}$")
axs[2].set_xlabel("Time [s]")
axs[2].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, ncol = 2)

# 4) Stiffness ratios vs baseline
#--new code Non-normality
baseline = store[BASELINE_KEY]["kappa"]
for be in BACKENDS:
    st = be["style"]
    #ratio = store[be["key"]]["kappa"] / baseline
    ratio = store[be["key"]]["departure"]
    #axs[3].plot(time, ratio, label=f"{be['label']} / {BASELINE_KEY}", **st)
    axs[3].plot(time, ratio, label=f"{be['label']}", **st)
axs[3].set_ylabel("Nonnormality, AA^T-A^tA")
axs[3].set_xlabel("Time [s]")
axs[3].legend()
#axs[3].set_ylim([0.1, np.nanmax([np.nanmax(store[b['key']]['kappa']/baseline) for b in BACKENDS]) * 10.0])

# Grids and reference line at 10^0 = 1
for ax in axs:
    ax.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
    ax.grid(True, which="minor", linestyle="--", linewidth=0.3, alpha=0.5)
axs[3].axhline(1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.8)

plt.tight_layout()
plt.subplots_adjust(right=0.78)
plt.savefig("stiffness_formulation.png", dpi=300, bbox_inches="tight")

# Iteration plots
fig, axs = plt.subplots(3, 1, figsize=(15, 10), sharex=True)

# Linear iterations per step
for be in BACKENDS:
    st = be["style"]
    np.save(f"time_{be['key']}.npy", time)
    np.save(f"linear_iterations_{be['key']}.npy", store[be["key"]]["lin_iter"])
    total = np.nansum(store[be["key"]]["lin_iter"])
    axs[0].plot(time, store[be["key"]]["lin_iter"],
                label=f"{be['label']} total = {int(total)}", **st)
axs[0].set_ylabel("Linear iterations\nper time step")
axs[0].set_xlabel("Time [s]")
axs[0].legend(ncol=2)

# Nonlinear iterations per step
for be in BACKENDS:
    st = be["style"]
    axs[1].plot(time, store[be["key"]]["newt_iter"],
                label=f"{be['label']}", **st)
axs[1].set_ylabel("Nonlinear iterations\nper time step")
axs[1].set_xlabel("Time [s]")
axs[1].legend(ncol=2)

# Linear per nonlinear
for be in BACKENDS:
    st = be["style"]
    num = store[be["key"]]["lin_iter"].astype(float)
    den = store[be["key"]]["newt_iter"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(den > 0, num / den, np.nan)
    axs[2].plot(time, ratio, label=f"{be['label']}", **st)
axs[2].set_ylabel("Linear iterations\nper nonlinear step")
axs[2].set_xlabel("Time [s]")
axs[2].legend(ncol=2)

plt.tight_layout()
plt.show()

