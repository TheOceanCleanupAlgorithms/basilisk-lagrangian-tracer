/**
# Particle positions and velocities as text tables

One file per quantity, as before:

~~~
<prefix>_particle_position_x.dat        <prefix>_particle_velocity_x.dat
<prefix>_particle_position_y.dat        <prefix>_particle_velocity_y.dat
<prefix>_particle_slip_velocity_x.dat   <prefix>_particle_slip_velocity_y.dat
(and _z with dimension 3 or the multilayer solver)
~~~

The slip-velocity files exist only for particles with a velocity of
their own (inertial particles). A tracer moves with the fluid by
definition, so for a list of tracer particles (one of
`tracer_particles`) they are neither opened nor written.

Each call of `particles_kinematics()` adds one row: the time, then one
column per particle tag, so that **column `n + 1` is always the particle
tagged `n`**, whichever rank holds it and whether or not other
particles were lost or released meanwhile. Missing particles (lost, not
released yet) are `nan`. Values are tab-separated, in `%.6e` (7
significant digits), the time with 10 digits so that it stays exact to
the second over long runs, and `numpy.loadtxt()` reads the files as
before.

## How it is written

Every value is a fixed-width field of 13 characters followed by a tab
(a newline after the last column), and the time a field of 17
characters and a tab, so the field of tag `n` starts at byte
`18 + 14*n` of its row. To fit, values below 1e-99 in magnitude are
written as 0 and values of 1e100 or more as `inf`. With MPI each rank formats its own
particles and writes them in place, through collective MPI-IO: nothing
is gathered on rank 0 and nothing is sorted globally, so time and
memory per rank follow the particles it holds, and the files are the
same at any number of ranks. The files stay open between calls.

## What is written

* **Position.** For tracer particles called after stage 1 (a `last`
  event defined after including [tracer-particles.h](), at an even
  step), the position at time `t`
  (`p().locn`). Otherwise `p().x`, which for tracer particles is a
  Runge–Kutta stage estimate: a warning says so, once. A stage
  position just outside a periodic domain is written wrapped into it,
  as `particle_boundary()` would, so that the files are the same at any
  number of ranks. Longitudes are in
  the convention of the domain (e.g. [105:255]); with `WCP_LON_WRAPPED`
  defined to 1 (before including this file, or `-DWCP_LON_WRAPPED=1`),
  they are written in [-180:180) (as `lon_to_wrapped()` in
  [geography.h](), which is not needed).
  This changes the files only: the particles keep the domain's
  convention. With `PARTICLE_SURFACE` (surface tracers, see
  [particle.h]()), `z` is the surface height `eta` at the particle.
* **Velocity.** Inertial particles (`p().inertia`): their own velocity
  `p().u`. Tracer particles: the fluid velocity at the written
  position, in the units of `u` (m/s with `PARTICLE_GEOGRAPHY` too, not
  degrees per second).
* **Slip velocity** (inertial particles only, see above). Fluid minus
  particle velocity at the written position. With the multilayer solver
  the vertical fluid velocity is only known with `layered/nh.h`: tracer
  particles then move with `w` (see [tracer-particles.h]()), and their
  vertical position and velocity are written. Without `nh.h`, define
  `WCP_VERTICAL` as a vertical velocity field if you have one, otherwise
  the vertical slip of inertial particles is `nan`.

The fluid velocity is interpolated by the rank holding the particle,
with `particle_locate()`, so it is correct under MPI. With the
multilayer solver, it is that of the layer holding the particle.

## Usage

~~~literatec
#include "tracer-particles.h"
#include "write_class_particles.h"

TextWriter out;
event init (t = 0) {
  ...                                  // create the particles
  out = init_particle_files ("w", P, "run");
}
event write (t += 86400., last) {      // after the particle events
  particles_kinematics (out, t);
}
event close (t = end)
  close_particle_files (&out);
~~~

The number of columns is fixed when the files are opened: `ntot`, or
the number of tags issued so far (`p_ntags`). With releases during the
run, pass the final number of tags. With `mode = "a"` and existing
files, it is read from their first line. All functions are collective. */

#define WCP_W 14  // bytes per value, separator included: "%13.6e" + tab
#define WCP_WT 18 // bytes of the time field: "%17.9e" + tab

enum {
  WCP_PX, WCP_PY, WCP_UX, WCP_UY, WCP_SX, WCP_SY,
#if dimension == 3 || LAYERS
  WCP_PZ, WCP_UZ, WCP_SZ,
#endif
  WCP_NF
};

static const char * wcp_names[] = {
  "position_x", "position_y", "velocity_x", "velocity_y",
  "slip_velocity_x", "slip_velocity_y",
  "position_z", "velocity_z", "slip_velocity_z"
};

typedef struct {
  Particles particles;
  long N;           // columns (tags)
  bool slip;        // slip-velocity files (not for tracer particles)
#if _MPI
  particle_file f[WCP_NF];
#else
  FILE * f[WCP_NF];
#endif
} TextWriter;

extern vector u;

/**
The number of tags of `P`: its counter, or the largest tag + 1. */

static long wcp_ntags (Particles P)
{
  if (p_ntags[P] && p_ntags[P] != P_NTAGS_UNKNOWN)
    return p_ntags[P];
  long unsigned int n = 0;
  foreach_particle_in (P)
    if (p().tag + 1 > n)
      n = p().tag + 1;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &n, 1, MPI_UNSIGNED_LONG, MPI_MAX,
		 MPI_COMM_WORLD);
#endif
  return n;
}

/**
The slip velocity only exists for particles with a velocity of their
own: `wcp_has_slip()` is false for a list of tracer particles. */

static bool wcp_is_slip (int k)
{
  return k == WCP_SX || k == WCP_SY
#if dimension == 3 || LAYERS
    || k == WCP_SZ
#endif
    ;
}

static bool wcp_has_slip (Particles Q)
{
#ifdef TRACER_TAGS
  if (tracer_particles)
    foreach_P_in_list (tracer_particles)
      if (P == Q)
	return false;
#endif
  return true;
}

/**
The columns of an existing file, from the length of its first line. */

static long wcp_columns (const char * name)
{
  FILE * fp = fopen (name, "r");
  if (!fp)
    return -1;
  long len = 0;
  int c;
  while ((c = fgetc (fp)) != EOF && c != '\n')
    len++;
  fclose (fp);
  return c == '\n' && len + 1 > WCP_WT ? (len + 1 - WCP_WT)/WCP_W : -1;
}

TextWriter init_particle_files (char * mode, Particles P, char * prefix,
				long ntot = 0)
{
  TextWriter w;
  w.particles = P;
  w.N = ntot ? ntot : wcp_ntags (P);
  w.slip = wcp_has_slip (P);
  bool append = (mode[0] == 'a');
  char name[WCP_NF][256];
  for (int k = 0; k < WCP_NF; k++)
    snprintf (name[k], 256, "%s_particle_%s.dat", prefix, wcp_names[k]);
  if (append) {
    long N = -1;
    if (pid() == 0)
      N = wcp_columns (name[0]);
#if _MPI
    MPI_Bcast (&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
    if (N > 0)
      w.N = N;
  }
  for (int k = 0; k < WCP_NF; k++) {
    if (!w.slip && wcp_is_slip (k))
      continue; // no file: see wcp_has_slip()
#if _MPI
    if (MPI_File_open (MPI_COMM_WORLD, name[k],
		       MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
		       &w.f[k]) != MPI_SUCCESS) {
      fprintf (stderr, "init_particle_files(): cannot open '%s'\n", name[k]);
      exit (1);
    }
    if (!append)
      MPI_File_set_size (w.f[k], 0);
#else
    if (!(w.f[k] = fopen (name[k], append ? "a" : "w"))) {
      perror (name[k]);
      exit (1);
    }
#endif
  }
  return w;
}

void close_particle_files (TextWriter * w)
{
  for (int k = 0; k < WCP_NF; k++) {
    if (!w->slip && wcp_is_slip (k))
      continue;
#if _MPI
    MPI_File_close (&w->f[k]);
#else
    fclose (w->f[k]);
#endif
  }
}

/**
A value as a field of exactly `WCP_W` bytes: 13 characters, right
justified (`nan` included), and the separator. `%.6e` needs 13
characters with a two-digit exponent, hence the clamping. The time
field is wider. */

static void wcp_field (char * dst, double v, bool last)
{
  char tmp[32];
  if (fabs (v) < 1e-99)
    v = 0.;
  else if (fabs (v) >= 1e100)
    v = v > 0. ? HUGE_VAL : - HUGE_VAL;
  snprintf (tmp, sizeof(tmp), "%13.6e", v);
  memcpy (dst, tmp, WCP_W - 1);
  dst[WCP_W - 1] = last ? '\n' : '\t';
}

static void wcp_time_field (char * dst, double t)
{
  char tmp[32];
  snprintf (tmp, sizeof(tmp), "%17.9e", t);
  memcpy (dst, tmp, WCP_WT - 1);
  dst[WCP_WT - 1] = '\t';
}

typedef struct {
  long unsigned int tag;
  double v[WCP_NF];
} wcp_values;

static int wcp_compare (const void * a, const void * b)
{
  long unsigned int ta = ((wcp_values *) a)->tag, tb = ((wcp_values *) b)->tag;
  return ta < tb ? -1 : ta > tb;
}

/**
The values of one particle, computed by the rank holding it. */

static void wcp_particle (particle * q, double t, double * v, bool * stage)
{
  coord pos = {q->x, q->y, q->z};
#ifdef TRACER_TAGS
  if (!q->inertia) {
    if (tp_position_time == t) {
      pos.x = q->locn.x, pos.y = q->locn.y; // position at t
#if TP_VERTICAL
      pos.z = q->locn.z;
#endif
    }
    else
      *stage = true;
  }
#endif
  /* The periodic wrapping of particle_boundary(): a stage position just
     outside the domain is written inside it, at any number of ranks
     (wcp_send_to_owners() wraps p().x the same way on trees with MPI). */
  coord o = {X0, Y0, Z0};
  foreach_dimension() {
    if (pos.x < o.x)
      pos.x += L0;
    else if (pos.x > o.x + L0)
      pos.x -= L0;
  }
  coord f = {nodata, nodata, nodata}, c = pos; // fluid velocity at pos
  Point point = particle_locate (&c);
  if (point.level >= 0) {
    double xq = c.x, yq = c.y, zq = c.z; // not rotated by foreach_dimension()
    foreach_dimension()
      f.x = interpolate_linear (point, u.x, xq, yq, zq);
#if LAYERS && defined(WCP_VERTICAL)
    f.z = interpolate_linear (point, WCP_VERTICAL, xq, yq, zq);
#elif TP_VERTICAL
    f.z = interpolate_linear (point, w, xq, yq, zq); // layered/nh.h
#endif
  }
  coord vel = q->inertia ? q->u : f;
#if LAYERS && dimension == 2 && !TP_VERTICAL
  if (!q->inertia)
    vel.z = 0.; // tracer particles do not move vertically
#endif
  v[WCP_PX] = pos.x, v[WCP_PY] = pos.y;
#if WCP_LON_WRAPPED // as lon_to_wrapped() in geography.h, not included
  v[WCP_PX] = ((pos.x + 180) - 360*floor ((pos.x + 180)/360)) - 180;
#endif
  v[WCP_UX] = vel.x, v[WCP_UY] = vel.y;
  v[WCP_SX] = q->inertia ? f.x - vel.x : 0.;
  v[WCP_SY] = q->inertia ? f.y - vel.y : 0.;
#if dimension == 3 || LAYERS
  v[WCP_PZ] = pos.z, v[WCP_UZ] = vel.z;
#if LAYERS && PARTICLE_SURFACE
  v[WCP_PZ] = particle_surface_height (pos.x, pos.y); // surface tracers
#endif
  v[WCP_SZ] = q->inertia ? f.z - vel.z : 0.;
#endif
  for (int k = 0; k < WCP_NF; k++)
    if (v[k] == nodata || fabs (v[k]) >= nodata/2.)
      v[k] = NAN;
}

/**
On a tree with MPI, a rank finds a particle's cell only among its own
cells (on the multigrid it also accepts the first layer of ghost cells,
see `particle_locate()` in [particle.h]()). A written position that has
crossed into a neighbouring subdomain since the last exchange then had
no cell, and its velocity was written as `nan`. That happens when the
written position is `p().x`: a Runge–Kutta stage estimate for tracer
particles (a write at an odd step), or the position of inertial
particles. The particles are then first sent to the rank owning that
position: the periodic wrapping of `particle_boundary()` (without the
multilayer reflection), then `update_mpi()`. Only which rank holds a
particle changes, so trajectories are unchanged. At even steps tracers
are written at `p().locn`, which step 1 has just exchanged, so nothing
is done. Collective. */

#if _MPI && !MULTIGRID_MPI
static void wcp_send_to_owners (Particles P)
{
  coord o = {X0, Y0, Z0};
  foreach_particle_in (P)
    foreach_dimension() {
      if (p().x < o.x)
	p().x += L0;
      else if (p().x > o.x + L0)
	p().x -= L0;
    }
  update_mpi (P);
}
#endif

void particles_kinematics (TextWriter wr, const double t)
{
  Particles P = wr.particles;
  long N = wr.N;
  if (N <= 0)
    return;
#if _MPI && !MULTIGRID_MPI // trees: see wcp_send_to_owners()
  bool at_x = true;
# ifdef TRACER_TAGS
  at_x = wr.slip || tp_position_time != t; // tracers at even steps: locn
# endif
  if (at_x)
    wcp_send_to_owners (P);
#endif
  boundary ((scalar *){u});
#if TP_VERTICAL
  boundary ({w});
#endif
  particle_locate_boundary();
#if LAYERS && PARTICLE_SURFACE
  boundary ({eta}); // for z, the surface height
#endif
  long unsigned int n = pn[P];
  wcp_values * r = malloc ((n + 1)*sizeof(wcp_values));
  bool stage = false;
  foreach_particle_in (P) {
    r[_j_particle].tag = p().tag;
    wcp_particle (&p(), t, r[_j_particle].v, &stage);
  }
  qsort (r, n, sizeof(wcp_values), wcp_compare);
  long unsigned int nl = 0; // our particles with a column
  while (nl < n && r[nl].tag < (long unsigned int) N)
    nl++;

  static bool warned = false;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &stage, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#endif
  if (stage && !warned && pid() == 0) {
    fprintf (stderr, "particles_kinematics(): at t = %g, tracer positions are"
	     " Runge-Kutta stage estimates: call it from a 'last' event at"
	     " even steps\n", t);
    warned = true;
  }

  long rowsize = WCP_WT + N*WCP_W;
  char * buf = malloc ((nl + 1)*WCP_W), tfield[WCP_WT];
  wcp_time_field (tfield, t);
#if _MPI
  // the row to write, from the size of the first file
  MPI_Offset size = 0;
  if (pid() == 0)
    MPI_File_get_size (wr.f[0], &size);
  MPI_Bcast (&size, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
  MPI_Offset offset = (size/rowsize)*rowsize;
  long unsigned int nglobal = nl;
  MPI_Allreduce (MPI_IN_PLACE, &nglobal, 1, MPI_UNSIGNED_LONG, MPI_SUM,
		 MPI_COMM_WORLD);
  MPI_Aint * displ = malloc ((nl + 1)*sizeof(MPI_Aint));
  for (long unsigned int j = 0; j < nl; j++)
    displ[j] = WCP_WT + r[j].tag*WCP_W;
  MPI_Datatype field, view = MPI_BYTE;
  MPI_Type_contiguous (WCP_W, MPI_BYTE, &field);
  MPI_Type_commit (&field);
  if (nl) {
    MPI_Type_create_hindexed_block (nl, 1, displ, field, &view);
    MPI_Type_commit (&view);
  }
  // missing particles: each rank fills its share of the row with nan
  long first = N*pid()/npe(), last = N*(pid() + 1)/npe();
  char * nanbuf = NULL;
  if (nglobal < (long unsigned int) N) {
    nanbuf = malloc ((last - first + 1)*WCP_W);
    for (long j = first; j < last; j++)
      wcp_field (nanbuf + (j - first)*WCP_W, NAN, j == N - 1);
  }
  for (int k = 0; k < WCP_NF; k++) {
    if (!wr.slip && wcp_is_slip (k))
      continue;
    if (nanbuf) {
      MPI_File_set_view (wr.f[k], 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
      MPI_File_write_at_all (wr.f[k], offset + WCP_WT + first*WCP_W, nanbuf,
			     (last - first)*WCP_W, MPI_BYTE, MPI_STATUS_IGNORE);
      MPI_File_sync (wr.f[k]);
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_File_sync (wr.f[k]);
    }
    MPI_File_set_view (wr.f[k], 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    if (pid() == 0)
      MPI_File_write_at (wr.f[k], offset, tfield, WCP_WT, MPI_BYTE,
			 MPI_STATUS_IGNORE);
    for (long unsigned int j = 0; j < nl; j++)
      wcp_field (buf + j*WCP_W, r[j].v[k], r[j].tag == (long unsigned int) N - 1);
    MPI_File_set_view (wr.f[k], offset, MPI_BYTE, view, "native", MPI_INFO_NULL);
    MPI_File_write_all (wr.f[k], buf, nl*WCP_W, MPI_BYTE, MPI_STATUS_IGNORE);
  }
  if (nl)
    MPI_Type_free (&view);
  MPI_Type_free (&field);
  free (displ), free (nanbuf);
#else // !_MPI
  char * row = malloc (rowsize);
  memcpy (row, tfield, WCP_WT);
  for (int k = 0; k < WCP_NF; k++) {
    if (!wr.slip && wcp_is_slip (k))
      continue;
    for (long j = 0; j < N; j++)
      wcp_field (row + WCP_WT + j*WCP_W, NAN, j == N - 1);
    for (long unsigned int j = 0; j < nl; j++)
      wcp_field (row + WCP_WT + r[j].tag*WCP_W, r[j].v[k], r[j].tag == N - 1);
    fwrite (row, 1, rowsize, wr.f[k]);
    fflush (wr.f[k]);
  }
  free (row);
#endif // !_MPI
  free (buf), free (r);
}
