# Lagrangian tracer particles for Basilisk, with MPI

Passive tracer particles for [Basilisk](http://basilisk.fr): each
particle moves with the fluid velocity interpolated at its position,
with a third-order Runge–Kutta scheme over two solver steps. They work
with the usual solvers (`navier-stokes/centered.h`, `stream.h`, ...),
with the multilayer solver (`layered/hydro.h`, `layered/nh.h`), on the
multigrid and on trees (quadtree, octree), serial or with MPI.

On the multigrid, particles share the fluid's domain decomposition:
each rank advances the particles in its own subdomain, and a particle
that crosses into another subdomain is sent to that neighbour only (at
most 8 messages per rank in 2D, whatever the number of ranks). Results
are the same, bit for bit, at any number of ranks. On trees, the
original all-gather exchange is used.

## Files

| file | needed | what it is |
|---|---|---|
| `particle.h` | always | particle lists, iterators, MPI exchange, dump and restore |
| `particle_classic.h` | always | functions for the non-layered solvers |
| `particle_multilayer.h` | always | functions for the multilayer solver |
| `tracer-particles.h` | always | tracer particles: time stepping, interpolation, initialisation, tags, releases, binary output |
| `write_class_particles.h` | optional | positions and velocities as text tables |
| `geography.h` | optional | metres/degrees conversions and longitude conventions |
| `unit-test/` | no | the tests (see `unit-test/README.md`) |

Keep both `particle_classic.h` and `particle_multilayer.h`, whatever the
solver: `particle.h` includes one or the other, and qcc needs both files
to exist.

## Using it

These headers replace the particle headers that Basilisk used to
provide. Either copy the four needed headers (plus the optional ones)
next to your `.c` file, or point qcc to this folder with
`-I/path/to/this/repo`.

```c
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "write_class_particles.h"       // optional: text output

Particles P;
TextWriter out;

event init (t = 0) {
  P = init_tp_square (n = 100, l = 0.8);  // 100 x 100 particles
  out = init_particle_files ("w", P, "run");
}

event output (t += 1., last) {           // last: after the particle events
  particles_kinematics (out, t);
}

event end (t = 10, last) {
  close_particle_files (&out);
}
```

Build:

```bash
qcc -O2 -Wall run.c -o run -lm
```

```bash
CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 run.c -o run-mpi -lm
```

```bash
mpirun -np 4 ./run-mpi
```

On macOS add `-D_DARWIN_C_SOURCE` to a hand-written `CC99` (the tree
grid needs it). Programs that include `view.h` also need
`-L$BASILISK/gl -lglutils -lfb_tiny`, in that order.

### Creating particles

All these functions are collective: every rank must call them, even a
rank that holds no particles.

| function | particles |
|---|---|
| `init_tp_square (n, xm, ym, l, zp)` | `n × n` on a square of side `l` centred on `(xm, ym)` |
| `init_tp_circle (n, xm, ym, l, zp)` | `n` on a circle |
| `init_tp_cells ()` | one at the centre of each cell (of each layer with the multilayer solver) |
| `release_tracer_particles (P, pos, n)` | adds `n` particles at `pos[0..n-1]` into list `P` during the run; every rank passes the same positions |
| `new_tracer_particles (n)` + `set_particle_attributes (P)` | your own placement |

Each particle has a **tag** that identifies it for its whole life,
whatever rank holds it. Tags are given at creation and are the same at
any number of ranks (on the multigrid; on trees, `init_tp_cells()` tags
in local order). `tp_ntags (P)` is the number of tags issued.

### Output

* **Text tables** (`write_class_particles.h`): one file per quantity
  (`<prefix>_particle_position_x.dat`, `..._velocity_x.dat`, ...), one
  row per call, one column per tag (`nan` before a release or after a
  loss), values in `%.6e`. Call `particles_kinematics()` from a `last`
  event at even steps to write the positions at time `t` exactly;
  otherwise it writes Runge–Kutta stage estimates and warns once.
  There are no slip-velocity files for tracer particles: a tracer has no
  slip.
* **Binary** (`tracer-particles.h`): `output_tracer_particles (P, fname,
  time, ntot)` appends one row of positions, `xyz[ntot][3]` doubles in
  tag order, after a header of two `int64` (`ntot`, 3).
* **Restart**: `pdump (fname, Dmem = true)` and `prestore (fname)`, with
  MPI-IO. A dump can be restored on any number of ranks. Use
  `Dmem = true` to keep tags and the Runge–Kutta state.

With MPI, all outputs are written in parallel (MPI-IO): nothing is
gathered on one rank.

### Solvers and switches

| case | what to do |
|---|---|
| usual solvers (`centered.h`, `stream.h`, ...) | nothing |
| multilayer, `layered/hydro.h` | nothing: each particle takes the velocity of the layer holding its height `z`; `z` stays fixed (this solver has no vertical velocity) |
| multilayer, `layered/nh.h` | include `layered/nh.h` **before** `tracer-particles.h`: particles also move vertically with `w` |
| particles that stay at the free surface | `#define PARTICLE_SURFACE 1` before `tracer-particles.h`: top layer always, `z` is the surface height `eta` in the outputs |
| `spherical.h` (longitude, latitude in degrees, `u` in m/s) | `#define PARTICLE_GEOGRAPHY 1` before `tracer-particles.h`, and set `Radius` in metres |
| text output with longitudes in [-180:180) | `#define WCP_LON_WRAPPED 1` before `write_class_particles.h` |

With `PARTICLE_GEOGRAPHY`, positions stay in the longitude convention of
the domain (e.g. [105:255] for a domain across the date line). For
positions given in another convention, `lon_to_domain()` in
`geography.h` converts them.

### MPI notes

* **Rank counts on the multigrid**: `n × 4^i` in 2D (`n × 8^i` in 3D),
  with `n` the product of the `dimensions()` you set: 1, 4, 16, 64 for
  a square domain, 2, 8, 32, 128 for `dimensions (2, 1)`.
* **Several lists**: they are exchanged together, in one message per
  neighbour.
* **Particles placed by hand on one rank**: call
  `mpi_redistribute_particles (P)` once after placing them.
* A message "forwarded ... beyond its neighbours" means a particle
  moved more than one subdomain in a step: the time step is far too
  large.

## Tests

```bash
cd unit-test && ./run-tests.sh
```

76 checks on the multigrid and quadtree, serial and with MPI, against
exact solutions and serial runs. See `unit-test/README.md`, which also
explains the Stokes-wave test (`unit-test/stokes_wave/`).

## Limitations

* Particles are passive tracers: no inertia, no rise or settling
  velocity, no diffusion.
* With `layered/hydro.h` alone, particles below the surface do not move
  vertically.
* Positions are wrapped by `L0` in every direction, as in Basilisk's
  original particle code, even in a non-periodic direction: a particle
  that reaches a wall jumps to the opposite side.
* The MPI exchange uses Basilisk internals of the multigrid (its
  Cartesian communicator); build with a recent Basilisk.

## License

GNU General Public License, version 3 (GPL-3.0): see `LICENSE`. The
headers are derived from Basilisk's particle code, which is distributed
under the same license.
