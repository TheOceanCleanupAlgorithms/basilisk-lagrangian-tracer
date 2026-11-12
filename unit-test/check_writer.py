#!/usr/bin/env python
"""Checks of the text tables written by write_class_particles.h, for
run-tests.sh: check_writer.py <run directory>. Prints PASS/FAIL lines."""
import glob
import sys
import numpy as np

d = sys.argv[1]
px = np.loadtxt(f"{d}/w_particle_position_x.dat")
py = np.loadtxt(f"{d}/w_particle_position_y.dat")
t = np.loadtxt(f"{d}/tags.txt")  # list tag id x y: positions at the end
tag = t[:, 1].astype(int)
cols = {len(l.split()) for l in open(f"{d}/w_particle_position_x.dat")}
slip = glob.glob(f"{d}/w_particle_slip_velocity_*.dat")


def check(name, ok):
    print(("PASS " if ok else "FAIL ") + name)


check("writer: one column per tag in every row", cols == {1 + px.shape[1] - 1})
check("writer: last row = positions at the end (%.6e: 7 digits)",
      np.allclose(px[-1, 1 + tag], t[:, 3], rtol=1e-6, atol=1e-7) and
      np.allclose(py[-1, 1 + tag], t[:, 4], rtol=1e-6, atol=1e-7))
released = px[:, 1 + int(tag.max())]  # the last tag, released last
check("writer: released particles are nan before their release",
      np.isnan(released[0]) and not np.isnan(released[-1]))
check("writer: no slip-velocity files for tracer particles", not slip)
