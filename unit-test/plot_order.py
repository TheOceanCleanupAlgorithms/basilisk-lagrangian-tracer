"""Order in time of the tracer particles: runs tconv and tcentered and
plots the error against the time step, in two columns.

    cd unit-test
    python3 plot_order.py [--out FILE]

It also writes the data of the figure next to it: FILE without its
extension + "_tconv.dat" (case, dt, error) and + "_tcentered.dat" (dt, t,
error).

It builds tconv and tcentered here, in unit-test/, with qcc and the
headers of the folder above (when they are missing or older than their
source), then runs them. qcc must be in the PATH, or $BASILISK set.

(a) tconv: rotation with rate 1 + cos 2t, linear in space (exact
    interpolation), error at t = 2. Third order with a prescribed u(t)
    and when u is updated like a solver does (solver=1), with equal or
    unequal steps; second order with P_RK2.
(b) tcentered: navier-stokes/centered.h with u = G t, which third order
    integrates exactly: the error at t = 1 must be round-off.

Needs numpy and matplotlib.
"""
import argparse
import os
import subprocess

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

HERE = os.path.dirname(os.path.abspath(__file__))  # unit-test/
ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--out", default=os.path.join(HERE, "order_in_time.png"),
                help="output figure (default: unit-test/order_in_time.png)")
args = ap.parse_args()

BASILISK = os.environ.get("BASILISK", os.path.expanduser("~/basilisk/src"))
ENV = dict(os.environ, BASILISK=BASILISK,
           PATH=os.environ.get("PATH", "") + os.pathsep + BASILISK)
HEADERS = [os.path.join(HERE, "..", h) for h in
           ("particle.h", "particle_classic.h", "particle_multilayer.h", "tracer-particles.h")]


def build(prog):
    """builds unit-test/<prog> from <prog>.c if missing or out of date"""
    exe, src = os.path.join(HERE, prog), os.path.join(HERE, prog + ".c")
    if os.path.exists(exe) and all(os.path.getmtime(exe) > os.path.getmtime(f)
                                   for f in [src] + HEADERS):
        return
    print(f"building {prog}")
    subprocess.run(["qcc", "-O2", "-disable-dimensions", "-grid=multigrid", "-I..",
                    prog + ".c", "-o", prog, "-lm"], cwd=HERE, env=ENV, check=True)


for prog in ("tconv", "tcentered"):
    build(prog)

DT = [0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125]


def run(prog, *a):
    """stdout lines of one run, without Basilisk's '#' lines"""
    out = subprocess.run([os.path.join(HERE, prog), *a], capture_output=True,
                         text=True, check=True).stdout
    return [l.split() for l in out.splitlines() if l and l[0] not in "#"]


def tconv(*opts):
    """time steps and largest position errors at t = 2"""
    rows = [run("tconv", f"dt={dt}", *opts)[0] for dt in DT]
    return np.array([float(r[0]) for r in rows]), np.array([float(r[2]) for r in rows])


def order(dt, e):
    return np.log(e[:-1]/e[1:])/np.log(dt[:-1]/dt[1:])


INK, MUTED, GRID = "#222222", "#666666", "#e6e6e6"
RK3, RK2 = "#0072B2", "#009E73"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3))

# (a) tconv
cases = [  # options, label, colour, line, marker, filled
    ((), "prescribed u(t), equal steps", RK3, "none", "o", False),
    (("r=1.5",), "prescribed u(t), unequal steps", RK3, "none", "s", False),
    (("solver=1",), "solver-like, equal steps", RK3, "-", "o", True),
    (("solver=1", "r=1.5"), "solver-like, unequal steps", RK3, "--", "s", True),
    (("rk2=1",), "P_RK2, equal steps", RK2, "-", "^", True),
]
DATA = []  # (case, dt, error), written next to the figure
for opts, label, col, ls, mk, filled in cases:
    dt, e = tconv(*opts)
    DATA += [(" ".join(opts) or "default", d, x) for d, x in zip(dt, e)]
    q = order(dt, e)
    print(f"tconv {' '.join(opts) or 'default':16s} order {' '.join(f'{x:.2f}' for x in q)}")
    if filled:
        a.plot(dt, e, linestyle=ls, color=col, lw=2, marker=mk, ms=7,
               label=f"{label} (order {q[-1]:.2f})")
    else:  # on top of the solver-like run, which gives the same errors
        a.plot(dt, e, linestyle="none", color=col, marker=mk, ms=13, mfc="none",
               mew=1.2, zorder=5, label=label)
x = np.array([0.0015, 0.03])
for p, c, xt, f in ((2, 0.25, 0.0016, 1/5), (3, 0.25, 0.0035, 1/6)):
    a.plot(x, c*x**p, color=MUTED, lw=1, ls=(0, (4, 3)))
    a.text(xt, c*xt**p*f, f"slope {p}", color=MUTED, fontsize=9)
a.set(xscale="log", yscale="log", xlabel="time step dt",
      ylabel="largest position error at t = 2",
      title="(a) rotating flow, rate 1 + cos 2t (tconv.c)")
a.legend(frameon=False, fontsize=8.5, loc="lower right")

# (b) tcentered
DTC = [dt for dt in DT if dt <= 0.02]
rows = [run("tcentered", f"{dt}")[0] for dt in DTC]
dtc = np.array([float(r[1]) for r in rows])
ec = np.array([abs(float(r[3])) for r in rows])
print("tcentered errors", " ".join(f"{e:.1e}" for e in ec))
b.axhspan(1e-18, 1e-14, color="#f0f0f0", zorder=0)
b.text(0.0013, 3e-18, "round-off", color=MUTED, fontsize=9)
b.plot(dtc, np.maximum(ec, 1e-18), "-", color=RK3, lw=2, marker="o", ms=7,
       label="particle position")
b.set(xscale="log", yscale="log", xlabel="time step dt",
      ylabel="position error at t = 1", ylim=(1e-18, 1e-2),
      title="(b) navier-stokes/centered.h, u = G t (tcentered.c)")
b.legend(frameon=False, fontsize=9, loc="center right")

for ax, ticks in ((a, DT), (b, DTC)):
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(args.out, dpi=150)
print(args.out)

# the data of the figure, next to it
base = os.path.splitext(args.out)[0]
with open(base + "_tconv.dat", "w") as fp:
    fp.write("# tconv: largest position error at t = 2\n# case dt error\n")
    for case, d, x in DATA:
        fp.write(f"{case.replace(' ', ','):16s} {d:.10g} {x:.6e}\n")
with open(base + "_tcentered.dat", "w") as fp:
    fp.write("# tcentered: position error at t = 1 (signed)\n# dt t error\n")
    for r in rows:
        fp.write(f"{r[1]} {r[2]} {r[3]}\n")
print(base + "_tconv.dat", base + "_tcentered.dat")
