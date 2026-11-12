/**
# Particle migration on the multigrid

Tracer particles are advected by a steady, analytic velocity field, so
that an MPI run must reproduce the serial one particle by particle.
Each particle carries an identifier `id` computed from its list and
initial lattice position, which (unlike `tag`) does not depend on the
domain decomposition.

Usage: `./tmig [key=value ...]` with

| key | default | meaning |
|---|---|---|
| `dt` | 0.005 | time step |
| `np` | 40 | `np x np` particles per list |
| `nl` | 1 | number of lists, each shifted by a fraction of a lattice step |
| `tend` | 2 | end time |
| `n` | 64 | grid size |
| `snap` | 0 | positions written every `snap` steps (0: final only) |
| `out` | 0 | rows written at the end, by gathering and by MPI-IO |
| `dump` | | file for `pdump()` at the end |
| `restore` | | file for `prestore()`, instead of placing particles |
| `release` | 0 | particles released into list 0 every 20 steps from step 10 |
| `cells` | 0 | 1: one particle per cell (`init_tp_cells()`) instead of the square |
| `write` | 0 | with `-DWRITER=1`: `particles_kinematics()` of list 0 every `write` steps (even) |
| `woff` | 0 | with `write`: write at steps `i % write == woff`; an odd `woff` writes Runge–Kutta stage positions |

With `-DBOX=1` the domain is the closed `dimensions(2,1)` box of the
North Pacific setup instead of the periodic unit square.

The grid is chosen on the command line, `qcc -grid=multigrid` or
`qcc -grid=quadtree` (uniform, at level `log2(n)`): the same flow on the
same cells, so the trajectories must be the same on both. There is no
grid `#include` here on purpose: qcc follows an `#include` even inside a
false `#if`, and two grid includes in one file gave a quadtree build that
crashed in its first `boundary()` call.

Output, per rank: `out-<pid>` (final `id x y`), `tags-<pid>` (final
`list tag id x y`), `traj-<pid>` (`i t id x y`) and, when restoring,
`restored-<pid>`. Rank 0 prints

~~~
ranks particles N steps ms_per_step recv_mean recv_max rounds msgs halo unlocated
~~~

per step, from step 2 on: `ms_per_step` is the time in the particle
events (slowest rank), `recv_*` the particles a rank receives from
others (mean and max over ranks), `rounds` the exchanges, `msgs` the
point-to-point messages a rank sends (mean), `halo` the interpolations
from the ghost layer (all ranks) and `unlocated` the particles that
could not be interpolated (all ranks, whole run). With `out`, a second
line

~~~
output ranks particles ms_gather ms_mpiio bytes_gather bytes_mpiio
~~~

gives the time per row (slowest rank) and the largest buffer on one
rank for each method. See [run-tests.sh](). */

#ifndef BOX
# define BOX 0
#endif
// no grid #include: qcc -grid=multigrid or -grid=quadtree (see above)
#include "run.h"
vector u[];

/**
The particle time is measured between `tic`, defined before the
particle events, and `toc`, defined after them. The particle events
are `last` events, which run after all others in the order they are
defined, so `tic`, `toc` and the output events are `last` too. */

static double wtime (void) {
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}

double tpart = 0., t0 = 0.;
long nsteps = 0;
bool counting = false;
double nrecv = 0., nmsgs = 0., ncoll = 0.;

event tic (i++, last) {
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  counting = (i >= 2);
  t0 = wtime();
}

#define ADD_PART_MEM coord u; coord u2; coord locn; long unsigned int tag; long id;
#include "tracer-particles.h"

/**
With MPI, the messages sent and received by the particle code are
counted through the MPI profiling interface. The original code
exchanges with `MPI_Allgatherv()` over all ranks, the neighbour
exchange with `MPI_Isend()` and `MPI_Recv()` on the multigrid's
Cartesian communicator (Basilisk's own ghost-cell messages use
`MPI_COMM_WORLD`). The current code also keeps its own counters
(`particle_stats`). */

#if _MPI
int MPI_Recv (void * buf, int count, MPI_Datatype type, int source,
	      int tag, MPI_Comm comm, MPI_Status * status)
{
#if defined(PNB) && !defined(PARTICLE_STATS) // the first neighbour exchange
  if (counting && type == particle_mpi_type())
    nrecv += count;
#endif
  return PMPI_Recv (buf, count, type, source, tag, comm, status);
}

int MPI_Isend (const void * buf, int count, MPI_Datatype type, int dest,
	       int tag, MPI_Comm comm, MPI_Request * request)
{
  if (counting && comm != MPI_COMM_WORLD)
    nmsgs++;
  return PMPI_Isend (buf, count, type, dest, tag, comm, request);
}

int MPI_Allgatherv (const void * sbuf, int scount, MPI_Datatype stype,
		    void * rbuf, const int rcounts[], const int displs[],
		    MPI_Datatype rtype, MPI_Comm comm)
{
  if (counting) {
    double n = 0.;
    for (int r = 0; r < npe(); r++)
      if (r != pid())
	n += rcounts[r];
    nrecv += n/sizeof(particle); // the original code sends bytes
    ncoll++;
  }
  return PMPI_Allgatherv (sbuf, scount, stype, rbuf, rcounts, displs,
			  rtype, comm);
}
#endif

#if WRITER
# include "write_class_particles.h"
#endif

#define MAXL 16
Particles P[MAXL];
int NP = 40, NL = 1, SNAP = 0, OUT = 0, RELEASE = 0, CELLS = 0, WRITE = 0, WOFF = 0;
double TEND = 2.;
char * DUMP = NULL, * RESTORE = NULL;

int main (int argc, char * argv[]) {
  DT = 0.005;
  N = 64;
  for (int a = 1; a < argc; a++) {
    char * v = strchr (argv[a], '=');
    if (!v) {
      fprintf (stderr, "usage: %s [key=value ...]\n", argv[0]);
      return 1;
    }
    *v++ = '\0';
    if (!strcmp (argv[a], "dt")) DT = atof (v);
    else if (!strcmp (argv[a], "np")) NP = atoi (v);
    else if (!strcmp (argv[a], "nl")) NL = atoi (v);
    else if (!strcmp (argv[a], "tend")) TEND = atof (v);
    else if (!strcmp (argv[a], "n")) N = atoi (v);
    else if (!strcmp (argv[a], "snap")) SNAP = atoi (v);
    else if (!strcmp (argv[a], "out")) OUT = atoi (v);
    else if (!strcmp (argv[a], "dump")) DUMP = v;
    else if (!strcmp (argv[a], "restore")) RESTORE = v;
    else if (!strcmp (argv[a], "release")) RELEASE = atoi (v);
    else if (!strcmp (argv[a], "cells")) CELLS = atoi (v);
    else if (!strcmp (argv[a], "write")) WRITE = atoi (v);
    else if (!strcmp (argv[a], "woff")) WOFF = atoi (v);
    else {
      fprintf (stderr, "unknown key '%s'\n", argv[a]);
      return 1;
    }
  }
  assert (NL >= 1 && NL <= MAXL);
#if BOX
  dimensions (2, 1);
  L0 = 2.; X0 = Y0 = 0.;
#else
  L0 = 1.; X0 = Y0 = 0.;
  periodic (right);
  periodic (top);
#endif
  run();
}

#if BOX
// cellular flow, closed streamlines inside [0:2]x[0:1]
double XM = 1., LP = 0.95;
#else
double XM = 0.5, LP = 0.9;
#endif

static void write_positions (const char * name, const char * mode,
			     bool step, int i, double t)
{
  FILE * fp = fopen (name, mode);
  for (int m = 0; m < NL; m++)
    foreach_particle_in (P[m])
      if (step)
	fprintf (fp, "%d %g %ld %.15g %.15g\n", i, t, p().id,
		 p().locn.x, p().locn.y);
      else
	fprintf (fp, "%ld %.15g %.15g\n", p().id, p().locn.x, p().locn.y);
  fclose (fp);
}

event init (t = 0) {
  foreach() {
#if BOX
    u.x[] = - sin(pi*x)*cos(pi*y);
    u.y[] =   cos(pi*x)*sin(pi*y);
#else
    u.x[] = 1. + 0.5*sin(2.*pi*y);
    u.y[] = 0.5 + 0.3*cos(2.*pi*x);
#endif
  }
  char name[80];
  if (RESTORE) {
    NL = prestore (RESTORE);
    tracer_particles = malloc (NL*sizeof(Particles));
    for (int m = 0; m < NL; m++)
      P[m] = tracer_particles[m] = m;
    sprintf (name, "restored-%d", pid());
    write_positions (name, "w", false, i, t);
  }
  else if (CELLS) {
    NL = 1;
    P[0] = init_tp_cells();
    foreach_particle_in (P[0])
      p().id = p().tag;
  }
  else {
    double d = LP/NP;
    for (int m = 0; m < NL; m++) {
      double xm = XM + m*d/NL;
      P[m] = init_tp_square (n = NP, xm = xm, ym = 0.5, l = LP);
      foreach_particle_in (P[m])
	p().id = m*NP*NP + lround ((x - (xm - LP/2. + d/2.))/d)*NP +
	  lround ((y - (0.5 - LP/2. + d/2.))/d);
    }
  }
  sprintf (name, "traj-%d", pid());
  FILE * fp = fopen (name, "w");
  fclose (fp);
}

event set_dtmax (i++) { dt = dtnext (DT); }

/**
Releases: `RELEASE` particles on a line across the domain, at a height
that changes from one release to the next. Every rank computes the same
positions. */

int nreleases = 0;
event releases (i = 10; i += 20) {
#ifdef TRACER_TAGS
  if (RELEASE > 0) {
    double Lx = L0, Ly = BOX ? L0/2. : L0;
    coord * pos = malloc (RELEASE*sizeof(coord));
    for (int k = 0; k < RELEASE; k++)
      pos[k] = (coord){X0 + Lx*(0.05 + 0.9*(k + 0.5)/RELEASE),
		       Y0 + Ly*(0.1 + 0.8*fmod (0.37*nreleases, 1.)), 0.};
    release_tracer_particles (P[0], pos, RELEASE);
    free (pos);
    nreleases++;
  }
#endif
}
event velocity (i++) {}

event toc (i++, last) {
  double el = wtime() - t0;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &el, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  if (counting) {
#ifdef PARTICLE_STATS
    if (nsteps == 0) { // count from step 2 on
      long u = particle_stats.unlocated;
      memset (&particle_stats, 0, sizeof(particle_stats));
      particle_stats.unlocated = u; // but report the whole run
    }
#endif
    tpart += el, nsteps++;
  }
  counting = false;
}

/**
## The particle writer

`particles_kinematics()` writes list 0 every `WRITE` steps, from a
`last` event after the particle events. The current writer is told the
final number of tags (the square plus the releases); the original one
writes whatever particles are present. */

#if WRITER
TextWriter writer;
double twrite = 0.;
long nwrites = 0;

event write_init (i = 0, last) {
  if (WRITE > 0) {
#ifdef WCP_W
    long steps = lround (TEND/DT), nrel = RELEASE > 0 && steps >= 10 ?
      (steps - 10)/20 + 1 : 0;
    writer = init_particle_files ("w", P[0], "w", ntot = NP*NP + nrel*RELEASE);
#else
    writer = init_particle_files ("w", P[0], "w");
#endif
  }
}

event write_out (i++, last) {
  if (WRITE > 0 && i % WRITE == WOFF) {
#if _MPI
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    double s = wtime();
    particles_kinematics (writer, t);
    double e = wtime() - s;
#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, &e, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    twrite += e, nwrites++;
  }
}
#endif // WRITER

/**
At even steps, after `tracer_particles_step1`, `p().locn` holds the
approximated position at time `t` (and `p().x` the first stage). */

event snapshot (i += 2, last) {
  if (SNAP > 0 && i % SNAP == 0) {
    char name[80];
    sprintf (name, "traj-%d", pid());
    write_positions (name, "a", true, i, t);
  }
}

#ifdef PARTICLE_STATS
/**
## Output: gathering on rank 0 and sorting, against MPI-IO

`output_gather_sort()` writes the same file format as
`output_tracer_particles()`, the way the output was done before: all
positions are gathered on rank 0, sorted by tag and written there. */

static double bytes_gather = 0., bytes_mpiio = 0.;

static void output_gather_sort (Particles P, char * fname, double time)
{
  long unsigned int n = pn[P];
  tagged_position * r = malloc ((n + 1)*sizeof(tagged_position)), * all = r;
  foreach_particle_in (P)
    r[_j_particle] = (tagged_position){p().tag, {x, y, z}};
  long unsigned int ntotal = n;
#if _MPI
  int nb = n*sizeof(tagged_position), cnt[npe()], dsp[npe()];
  MPI_Gather (&nb, 1, MPI_INT, cnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (pid() == 0) {
    dsp[0] = 0;
    for (int k = 1; k < npe(); k++)
      dsp[k] = dsp[k - 1] + cnt[k - 1];
    ntotal = (dsp[npe() - 1] + cnt[npe() - 1])/sizeof(tagged_position);
    all = malloc ((ntotal + 1)*sizeof(tagged_position));
  }
  MPI_Gatherv (r, nb, MPI_BYTE, all, cnt, dsp, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  if (pid() == 0) {
    qsort (all, ntotal, sizeof(tagged_position), tagged_position_compare);
    long head[2] = {0, 3};
    FILE * fp = fopen (fname, "r+b");
    if (fp)
      fread (head, sizeof(long), 2, fp);
    else {
      fp = fopen (fname, "w+b");
      head[0] = ntotal ? all[ntotal - 1].tag + 1 : 0;
      fwrite (head, sizeof(long), 2, fp);
    }
    long N = head[0];
    double * row = malloc ((1 + 3*N)*sizeof(double));
    row[0] = time;
    for (long j = 0; j < 3*N; j++)
      row[1 + j] = NAN;
    for (long unsigned int j = 0; j < ntotal && all[j].tag < N; j++)
      memcpy (row + 1 + 3*all[j].tag, all[j].c, 3*sizeof(double));
    fseek (fp, 0, SEEK_END);
    fwrite (row, sizeof(double), 1 + 3*N, fp);
    fclose (fp);
    free (row);
    bytes_gather = ntotal*sizeof(tagged_position) + (1 + 3*N)*sizeof(double);
  }
  if (all != r)
    free (all);
  free (r);
}
#endif // PARTICLE_STATS

event end (t = TEND, last) {
  char name[80];
  sprintf (name, "out-%d", pid());
  write_positions (name, "w", false, i, t);
  long unsigned int n = 0;
  for (int m = 0; m < NL; m++)
    n += pn[P[m]];
  double rmean = nrecv, rmax = nrecv, msgs = nmsgs;
#ifdef PARTICLE_STATS
  rmean = rmax = particle_stats.received;
  double halo = particle_stats.halo, unloc = particle_stats.unlocated;
  double rounds = particle_stats.rounds;
#else
  double halo = 0., unloc = 0., rounds = 0.;
#endif
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &n, 1, MPI_UNSIGNED_LONG, MPI_SUM,
		 MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &rmean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &rmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &msgs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &halo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &unloc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rmean /= npe(), msgs /= npe();
#ifndef PARTICLE_STATS
  /* Exchange rounds: an all-gather each in the original code; for the
     first neighbour exchange, the messages over the valid neighbours. */
# ifdef PNB
  int nbr[PNB], nvalid = 0;
  particle_neighbours (particle_cart_comm(), nbr);
  for (int k = 0; k < PNB; k++)
    nvalid += nbr[k] != MPI_PROC_NULL;
  rounds = nvalid ? nmsgs/nvalid : 0.;
# else
  rounds = ncoll;
# endif
#endif
#endif // _MPI
  if (pid() == 0)
    printf ("%d %lu %d %ld %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", npe(), n,
	    N, nsteps, nsteps ? 1e3*tpart/nsteps : 0.,
	    nsteps ? rmean/nsteps : 0., nsteps ? rmax/nsteps : 0.,
	    nsteps ? rounds/nsteps : 0., nsteps ? msgs/nsteps : 0.,
	    nsteps ? halo/nsteps : 0., unloc);

#ifdef PARTICLE_STATS
  /**
  The output comparison uses list 0. Its tags do not depend on the
  decomposition, so neither do the files. */

  if (OUT > 0) {
    double tg = 0., tm = 0.;
    if (pid() == 0)
      remove ("rows-gather.bin"), remove ("rows-mpiio.bin");
    for (int k = 0; k < OUT; k++) {
#if _MPI
      MPI_Barrier (MPI_COMM_WORLD);
#endif
      double s = wtime();
      output_gather_sort (P[0], "rows-gather.bin", t + k);
      double e = wtime() - s;
#if _MPI
      MPI_Allreduce (MPI_IN_PLACE, &e, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
#endif
      tg += e;
      s = wtime();
      output_tracer_particles (P[0], "rows-mpiio.bin", t + k);
      e = wtime() - s;
#if _MPI
      MPI_Allreduce (MPI_IN_PLACE, &e, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
      tm += e;
    }
    // buffers of output_tracer_particles(): sorted copy, packed positions, offsets
    bytes_mpiio = pn[P[0]]*(sizeof(tagged_position) + 3*sizeof(double)
#if _MPI
			    + sizeof(MPI_Aint)
#endif
			    );
    long unsigned int n0 = pn[P[0]];
#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, &bytes_mpiio, 1, MPI_DOUBLE, MPI_MAX,
		   MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, &n0, 1, MPI_UNSIGNED_LONG, MPI_SUM,
		   MPI_COMM_WORLD);
#endif
    if (pid() == 0)
      printf ("output %d %lu %.6g %.6g %.6g %.6g\n", npe(), n0,
	      1e3*tg/OUT, 1e3*tm/OUT, bytes_gather, bytes_mpiio);
  }
#endif // PARTICLE_STATS

  sprintf (name, "tags-%d", pid());
  FILE * fp = fopen (name, "w");
  for (int m = 0; m < NL; m++)
    foreach_particle_in (P[m])
      fprintf (fp, "%d %lu %ld %.15g %.15g\n", m, p().tag, p().id,
	       p().locn.x, p().locn.y);
  fclose (fp);
#ifdef TRACER_TAGS
  long unsigned int ntags = tp_ntags (P[0]);
  if (pid() == 0)
    printf ("ntags %lu\n", ntags);
#endif

#if WRITER
  if (WRITE > 0) {
    // largest buffer on one rank: the gathered particles for the original
    // writer (on every rank), the local particles for the current one
    long unsigned int nloc = pn[P[0]], ntotal = nloc;
#if _MPI
    MPI_Allreduce (&nloc, &ntotal, 1, MPI_UNSIGNED_LONG, MPI_SUM,
		   MPI_COMM_WORLD);
#endif
#ifdef WCP_W
    double bytes = nloc*(sizeof(wcp_values) + WCP_W + 8.);
    close_particle_files (&writer);
#else
    double bytes = (ntotal + nloc)*(double) sizeof(particle);
#endif
#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, &bytes, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    if (pid() == 0)
      printf ("writer %d %lu %ld %.6g %.6g\n", npe(), ntotal, nwrites,
	      nwrites ? 1e3*twrite/nwrites : 0., bytes);
  }
#endif // WRITER

  if (DUMP)
    pdump (DUMP, Dmem = true);
}
