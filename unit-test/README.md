# Unit tests of the particle headers

The headers to use are in the folder above (`particle.h`,
`tracer-particles.h`, `particle_classic.h`, `particle_multilayer.h`,
and optionally `write_class_particles.h`, `geography.h`). These tests
check them on the multigrid and quadtree grids, serial and with MPI,
against exact solutions and serial runs: they test the current code
only, with the current Basilisk.

## The test suite

```bash
./run-tests.sh
```

It builds every test program, runs them on 1 to 64 ranks and prints
each check and the totals (76 checks; the script fails if one fails).
Builds, raw results and the list of checks (`results/checks.txt`) go to
`$WORK` (default `~/Library/Caches/particle-tests`). `./run-tests.sh build` only builds.
It needs Basilisk (`$BASILISK`, default `~/basilisk/src`), MPI and
Python 3 with numpy (`PYTHON=...` to choose it).

| program | what it tests |
|---|---|
| `tmig.c` | migration, tags, releases, output, dump/restore and the text writer, on the multigrid (`qcc -grid=multigrid`) or a quadtree (`-grid=quadtree`) |
| `tml.c` | the multilayer solver: layer selection, surface tracers (`PARTICLE_SURFACE`) |
| `tsph.c` | positions in degrees with `spherical.h` (`PARTICLE_GEOGRAPHY`) |
| `tnh.c` | vertical motion with `layered/nh.h` in a Stokes wave, serial and MPI |

## The Stokes wave (`stokes_wave/`)

A forced deep-water Stokes wave (`ak = 0.11`, N = 256, 20 layers, 30
periods) with a column of tracer particles: their drift must follow the
Stokes drift `u_S(z) = ω a² k e^{2kz}` near the surface. It takes about
30 minutes; 4 ranks are enough. From `stokes_wave/`:

```bash
export BASILISK=$HOME/basilisk/src PATH=$PATH:$HOME/basilisk/src
```

Serial:

```bash
qcc -I../.. -D_DARWIN_C_SOURCE -O2 -Wall stokes_drift.c -o stokes_drift -L$BASILISK/gl -lglutils -lfb_tiny -lm
```

With MPI, in two steps (`qcc -source` returns status 2 even when it
works):

```bash
CC99='mpicc -std=c99 -D_DARWIN_C_SOURCE' qcc -I../.. -source -D_MPI=1 stokes_drift.c
```

```bash
mpicc -std=c99 -D_DARWIN_C_SOURCE -O2 -D_MPI=1 -I$BASILISK _stokes_drift.c -o stokes_drift-mpi -L$BASILISK/gl -lglutils -lfb_tiny -lm
```

```bash
mpirun -np 4 ./stokes_drift-mpi
```

`process_particles.py` then compares the drift of each particle with
the theory. `-D_DARWIN_C_SOURCE` is for macOS; the graphics libraries
(`-lglutils -lfb_tiny`) are for the movie.
