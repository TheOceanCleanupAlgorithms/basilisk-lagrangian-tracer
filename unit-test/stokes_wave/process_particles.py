#!/usr/bin/env python3
"""
Post-process Lagrangian tracer particles from stokes_drift.c and compare
their horizontal drift with the analytic Stokes drift u_S(z).

Reads the position files written by write-ml-particles.h:
    <prefix>_particle_position_x.dat   (row: t  x0 x1 ... xN-1)
    <prefix>_particle_position_y.dat
    <prefix>_particle_position_z.dat

For each particle it estimates the drift velocity as the least-squares
slope of x(t) (more robust than an endpoint difference and it averages
out the orbital oscillation), and plots it against the particle's
time-mean depth, together with u_S(z) = omega a^2 k exp(2 k z).

Usage:
    python3 process_particles.py [prefix] [--tmin T0units] [--tmax T0units]

Defaults match stokes_drift.c (ak=0.11, k=2*pi, h=0.5, g=1).
"""

import sys
import math
import argparse
import numpy as np

# ---------------------------------------------------------------- parameters
# (keep in sync with stokes_drift.c)
AK = 0.11
K  = 2.0 * math.pi
H  = 0.5
G  = 1.0

OMEGA = math.sqrt(G * K * math.tanh(K * H))
A     = AK / K
T0    = K * 1.0 / math.sqrt(G * K)     # wave period (L0 = 1)


def uS(z):
    """Analytic deep-water surface Stokes drift at depth z (z <= 0)."""
    return OMEGA * A * A * K * np.exp(2.0 * K * z)


def load(prefix, comp):
    """Load one position file as (t[nrow], vals[nrow, npart])."""
    data = np.loadtxt(f"{prefix}_particle_position_{comp}.dat")
    if data.ndim == 1:                 # single row
        data = data[None, :]
    return data[:, 0], data[:, 1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix", nargs="?", default="tracer")
    ap.add_argument("--tmin", type=float, default=0.0,
                    help="start of fit window, in units of T0")
    ap.add_argument("--tmax", type=float, default=None,
                    help="end of fit window, in units of T0")
    ap.add_argument("--out", default="particle_validation.png")
    args = ap.parse_args()

    t, X = load(args.prefix, "x")
    _, Z = load(args.prefix, "z")
    npart = X.shape[1]

    # fit window (in physical time)
    tmin = args.tmin * T0
    tmax = args.tmax * T0 if args.tmax is not None else t[-1]
    m = (t >= tmin) & (t <= tmax)
    if m.sum() < 2:
        sys.exit("fit window contains < 2 samples")

    tw = t[m]
    drift = np.empty(npart)
    zmean = np.empty(npart)
    for p in range(npart):
        # least-squares slope of x(t) over the window = mean drift velocity
        drift[p] = np.polyfit(tw, X[m, p], 1)[0]
        zmean[p] = Z[m, p].mean()      # time-mean depth

    ua = uS(zmean)
    print(f"# span {t[0]:.3f}..{t[-1]:.3f}  ({(t[-1]-t[0])/T0:.1f} T0), "
          f"fit {args.tmin:.1f}..{tmax/T0:.1f} T0, {npart} particles")
    print(f"# {'z_mean':>9} {'drift_sim':>12} {'uS_analytic':>12} {'ratio':>7}")
    for p in np.argsort(zmean):
        r = drift[p] / ua[p] if ua[p] else float("nan")
        print(f"  {zmean[p]:9.4f} {drift[p]:12.6e} {ua[p]:12.6e} {r:7.2f}")

    # ---------------------------------------------------------------- plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        zz = np.linspace(-H, -1e-4, 300)
        fig, ax = plt.subplots(figsize=(5.2, 5))
        ax.plot(uS(zz), zz, "k-", lw=2,
                label=r"analytic $u_S(z)=\omega a^2 k\,e^{2kz}$")
        ax.plot(drift, zmean, "o", c="C3", ms=6,
                label="Lagrangian particle drift (sim)")
        ax.set_xlabel("horizontal drift velocity")
        ax.set_ylabel(r"time-mean depth $\bar z$")
        ax.set_ylim(-H, 0.02)
        ax.legend(fontsize=9)
        ax.set_title("Particle drift vs analytic Stokes drift")
        fig.tight_layout()
        fig.savefig(args.out, dpi=130)
        print(f"# figure written to {args.out}")
    except ImportError:
        print("# matplotlib not available; skipped figure")


if __name__ == "__main__":
    main()
