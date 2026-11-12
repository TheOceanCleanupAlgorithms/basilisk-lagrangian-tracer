/**
# Langrangian particles

Tracer-particle coordinates $\mathbf{p_i}$ are advected by a vector field
$\mathbf{u}$: 

$$\frac{\mathrm{d}\mathbf{p_i}}{\mathrm{d}t} = \mathrm{u}|_{\mathbf{x} = \mathbf{p_i}}$$

## Definititions and types 

The generic `particle` class is used and tuned: For the
time-advancement scheme, additional storage for velocities and
locations is required. Furthermore, particles may be `tag`ged. The
tracer particles are stored in an array of groups of `particle`s. Note 
that `Particles` reference a list of `particle`s. 
 */
extern vector u;
#ifndef ADD_PART_MEM
#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; 
#endif
#include "particle.h"
#define TRACER_TAGS 1 // tags independent of the decomposition, releases

/**
## Positions in degrees (`PARTICLE_GEOGRAPHY`)

With [spherical.h](/src/spherical.h), positions are longitude and
latitude in degrees, while `u` is in metres per second (the units of
`Radius`). Defined to 1 before including this file, or with
`-DPARTICLE_GEOGRAPHY=1`, `PARTICLE_GEOGRAPHY` converts every
interpolated velocity into degrees per second, on the model's sphere,
at the position where it was interpolated:

$$\frac{d\lambda}{dt} = \frac{u}{R\cos\phi}\frac{180}{\pi}, \qquad
\frac{d\phi}{dt} = \frac{v}{R}\frac{180}{\pi}$$

`p().u` and `p().u2` are then in degrees per second, and the
Runge–Kutta stages are unchanged. Without it, positions and `u` are
taken to be in the same length unit, as before, and none of this is
compiled. The switch needs `spherical.h` ("`Radius` undeclared"
otherwise), but not [geography.h](). Positions stay in the convention
of the domain (e.g. [105:255]): for positions given in another one,
include `geography.h` and use its `lon_to_domain()`. (No header is
included here under the switch: qcc follows `#include` even inside a
false `#if`, which would make `geography.h` required everywhere.) */

/**
## Vertical motion (multilayer solver with `layered/nh.h`)

A passive tracer follows the fluid in all three directions. With the
non-hydrostatic multilayer solver ([layered/nh.h](/src/layered/nh.h),
which defines `NH`), the vertical velocity `w` exists: each particle
then takes `w` of its layer, interpolated like `u`, and `z` is advanced
in the three Runge–Kutta stages together with `x` and `y` (`p().u.z`,
`p().u2.z` and `p().locn.z` hold the vertical parts). This is what makes
the Stokes drift of a wave right: a particle held at a fixed height
crosses the layers as they move with the wave, instead of moving with
its own layer (PROJECT-LOG §3.17).

`layered/nh.h` must be included *before* this file, so that `NH` is
defined here; a run where `w` exists but particles do not use it
warns at the first step. With `layered/hydro.h` alone there is no `w`
and `z` stays fixed. Surface tracers (`PARTICLE_SURFACE`) do not use
`z` at all. */

#if LAYERS && NH && !PARTICLE_SURFACE
# define TP_VERTICAL 1
#else
# define TP_VERTICAL 0
#endif

#if PARTICLE_GEOGRAPHY

static void tp_rate (coord p, coord * v)
{
  static bool warned = false;
  if (Radius == 1. && !warned) {
    fprintf (stderr, "tracer-particles.h: PARTICLE_GEOGRAPHY with Radius = 1:"
	     " set Radius in the units of u (e.g. 6371220 m)\n");
    warned = true;
  }
  v->x *= 180./(pi*Radius*cos (p.y*pi/180.));
  v->y *= 180./(pi*Radius);
}
#endif
Particles * tracer_particles = NULL;
bool P_RK2 = false;

/**
## Time integration  
 
The used particle-advection scheme is the RK-3 variant of Sanderse and
Veldman (2019) for $\alpha \neq \frac{2}{3} \pm \mathtt{tol}$:

$$ \begin{array}{c|ccc} 
0 & 0 & & \\ 
\alpha & \alpha & & \\ 
1 &1+\frac{1- \alpha}{\alpha (3\alpha -2)} & -\frac{1- \alpha}{\alpha
(3\alpha -2)} &\\ 
\hline
 & \frac{1}{2}-\frac{1}{6\alpha} &
\frac{1}{6\alpha(1-\alpha)} & \frac{2-3\alpha}{6(1-\alpha)} \\
\end{array} $$

see:  

B Sanderse and AEP Veldman, *Constraint-consistent Runge-Kutta
methods for one-dimensional incompressible multiphase flow*,
J. Comp. Phys. (2019) [link to
JCP.](https://www.sciencedirect.com/science/article/pii/S0021999119300683).

This scheme is used to advance the solution in two solver iterations
between $t_{n-2}$ and $t_n$. In case $\alpha \approx \frac{2}{3}$,
*or* `P_RK2 == true`, we skip the last stage and switch to a generic
second-order-accurate scheme.
 */
double tol = 3.e-2;
double dtf[2]; //timesteps
/**
The user is supposed to add particles in some relevant event, perhaps
called `init`. 
 */
event defaults (t = 0);

event init (t = 0);

event set_dtmax (i++);

event velocity (i++);
/**
By default, second-order-in-space "linear" interpolation is used to
approximate $\mathrm{u}|_{\mathbf{x} = \mathbf{p_i}}$. This could be
overloaded with [higher-order](higher-order.h) methods.

e.g for a 3rd order accurate method:

~~~literatec
...
#include "higher-order.h" //This is in /sandbox/Antoonvh/
#define interpolate_linear interpolate_quadratic
#include "tracer-particles.h" 
...
~~~

Velocities may found and stored in `p().u` or `p().u2`.

Near the edge of a subdomain the interpolation stencil reads the ghost
cells of `u`. The interpolation runs outside of `foreach()` loops, so
Basilisk does not update them automatically: `boundary()` does it if
they are out of date and costs nothing otherwise.

With the multilayer solver, each particle takes the velocity of the
layer holding its height `z` (see `particle_locate()` in
[particle.h]()).
 */
void interpolate_vel (Particles p, int swtch) {
  boundary ((scalar *){u});
#if TP_VERTICAL
  boundary ({w});
#endif
  particle_locate_boundary();
  foreach_particle_in(p) {
    coord * v = swtch == 1 ? &p().u : &p().u2, q = {x, y, z};
    Point point = particle_locate (&q);
    if (point.level < 0) {//freeze
      fprintf (stderr, "#pid: %d could not locate p.tag: %ld.\n",
	       pid(), p().tag); //This should not happen
      particle_stats.unlocated++;
      foreach_dimension()
	v->x = 0;
#if TP_VERTICAL
      v->z = 0;
#endif
    } else {
      double xq = q.x, yq = q.y, zq = q.z; // not rotated by foreach_dimension()
      foreach_dimension()
	v->x = interpolate_linear (point, u.x, xq, yq, zq);
#if TP_VERTICAL
      v->z = interpolate_linear (point, w, xq, yq, zq); // w of the same layer
#endif
#if PARTICLE_GEOGRAPHY
      tp_rate (q, v);
#endif
    }
  }
}
/**
The three Runge-Kutta stages are defined below:
 */
void RK_step1 (Particles p, double dtf[2]) {
  interpolate_vel(p, 1);
  foreach_particle_in(p) { 
    foreach_dimension() {
      p().locn.x = p().x; //store the locations at t_n 
      p().x += dtf[0]*p().u.x;
    }
#if TP_VERTICAL
    p().locn.z = p().z;
    p().z += dtf[0]*p().u.z;
#endif
  }
}

void RK_step2 (Particles p, double dtf[2]) {
  interpolate_vel(p, 2);
  double a1 = -1, a2 = 2., h = dtf[1] + dtf[0];
  if (dtf[1] != dtf[0] || P_RK2) {
    double c = dtf[0]/h;
    if (fabs (c - 2./3.) > tol && !P_RK2)
      a2 = (c - 1.)/(c*(3.*c - 2.));
    else // Raltson's 2nd order method
      a2 = 1./(2*c);
    a1 = 1 - a2;
  }
  foreach_particle_in(p) {
    foreach_dimension() 
      p().x = p().locn.x + h*(a1*p().u.x + a2*p().u2.x);
#if TP_VERTICAL
    p().z = p().locn.z + h*(a1*p().u.z + a2*p().u2.z);
#endif
  }
}

void RK_step3 (Particles p, double dtf[2]) {
  double h = dtf[1] + dtf[0];
  double c = dtf[0]/h;
  if (fabs(c - 2./3.) > tol) {  // RK-3
    boundary ((scalar *){u});
#if TP_VERTICAL
    boundary ({w});
#endif
    particle_locate_boundary();
    double b1 = 0.5 - 1./(6.*c);
    double b2 = 1./(6.*c*(1. - c));
    double b3 = 1. - (b1 + b2);
    foreach_particle_in(p) {
      coord q = {x, y, z}, u3 = p().u2; // a particle in transit uses u2
      Point point = particle_locate (&q);
      if (point.level >= 0) {
	double xq = q.x, yq = q.y, zq = q.z; // not rotated by foreach_dimension()
	foreach_dimension()
	  u3.x = interpolate_linear (point, u.x, xq, yq, zq);
#if TP_VERTICAL
	u3.z = interpolate_linear (point, w, xq, yq, zq);
#endif
#if PARTICLE_GEOGRAPHY
	tp_rate (q, &u3);
#endif
      }
      foreach_dimension()
        p().x = p().locn.x + h*(b1*p().u.x + b2*p().u2.x + b3*u3.x);
#if TP_VERTICAL
      p().z = p().locn.z + h*(b1*p().u.z + b2*p().u2.z + b3*u3.z);
#endif
    }
  }
}
/**
The particles are advanced in time in various RK-stage events. Note
that the particle coordinates (`p().x`) only contain the proper
approximated location inbetween `step3` and `step1`. These locations
are otherwise stored in `p().locn`. 

All the lists are exchanged at once (`particles_boundary()`). Step 1
directly follows step 3, within the same iteration, and in between
the particles only move by the third-stage correction, a small
fraction of a cell. On the multigrid, step 1 therefore skips the MPI
exchange and interpolates particles that have just crossed into a
neighbour from the ghost layer, with the same result (see
`particle_locate()`): two exchanges per two steps instead of
three. Step 1 exchanges at `i = 0`, and whenever step 3 did not run. */

event tracer_particles_step3 (i += 2, last);
event tracer_particles_step3 (i = 2; i += 2) {
  if (!P_RK2 && tracer_particles) {
    particles_boundary (tracer_particles);
    foreach_P_in_list (tracer_particles)
      RK_step3 (P, dtf);
  }
}

void tp_insert_releases (void);

/**
Right after step 1, `p().locn` holds the position at time `t` (the
third-order one, at even steps) and `p().u` the velocity there.
`tp_position_time` records that time, for output functions such as
[write_class_particles.h](). */

double tp_position_time = - HUGE;

event tracer_particles_step1 (i += 2, last) {
#if LAYERS && !TP_VERTICAL && !PARTICLE_SURFACE
  if (i == 0 && lookup_field ("w").i >= 0 && pid() == 0)
    fprintf (stderr, "tracer-particles.h: a field 'w' exists but particles do"
	     " not move vertically: include layered/nh.h before"
	     " tracer-particles.h\n");
#endif
  dtf[0] = dt;
  tp_insert_releases();
  if (tracer_particles) {
    particles_boundary (tracer_particles, exchange = (i == 0 || P_RK2));
    foreach_P_in_list (tracer_particles)
      RK_step1 (P, dtf);
  }
  tp_position_time = t;
}

event tracer_particles_step2 (i += 2, last);
event tracer_particles_step2 (i = 1; i += 2) {
  dtf[1] = dt;
  if (tracer_particles) {
    particles_boundary (tracer_particles);
    foreach_P_in_list (tracer_particles)
      RK_step2 (P, dtf);
  }
}

event free_tracers (t = end) {
  free (tracer_particles);
  tracer_particles = NULL;
}

/**
## Initializing tracer particles 

`n` New tracer particles can be declared by calling
`new_tracer_particles (n)`. This is a building-block of
tracer-particle initialization.
*/
Particles new_tracer_particles (long unsigned int n) {
  Particles p = new_particles (n);
  int l = 0, t = 0;
  if (tracer_particles != NULL) {
    while (pn[l] != terminate_int) {
      if (l == tracer_particles[t]) 
	t++;
      l++;
    }
  }
  tracer_particles = realloc (tracer_particles, (t + 1)*sizeof(Particles));
  tracer_particles[t] = p;
  return p;
}

/**
## Tags

A tag identifies a particle for its whole life, whichever rank holds
it. To be the same at any number of ranks, it must be given when the
particle is created, from something all ranks share, and never from
where the particle ends up after distribution:

* `init_tp_square()` and `init_tp_circle()` tag on rank 0, in placement
  order, before distributing the particles;
* `init_tp_cells()` tags from the global cell index (on the multigrid);
* `release_tracer_particles()` tags in release order.

`p_ntags[P]` (in [particle.h]()) counts the tags issued for list `P`,
so new tags continue from there, and `output_tracer_particles()` knows
its number of columns.

`tag_particles()` numbers the particles of `p` in their current order,
rank after rank, and sets the counter. Called before the particles are
distributed (all on one rank, or each rank holding a share that does not
depend on the decomposition), the tags do not depend on the number of
ranks; called after, they do. */

void tag_particles (Particles p) {
  long unsigned int offset = 0, n = pn[p], ntot = n;
#if _MPI
  MPI_Exscan (&n, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (pid() == 0)
    offset = 0;
  MPI_Allreduce (&n, &ntot, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  foreach_particle_in(p)
    p().tag = _j_particle + offset;
  p_ntags[p] = ntot;
}

/**
`tp_ntags()` returns the tag counter of `P`. After restoring a dump
written without it, it is recomputed as the largest tag + 1
(collective). */

long unsigned int tp_ntags (Particles P)
{
  if (p_ntags[P] == P_NTAGS_UNKNOWN) {
    long unsigned int n = 0;
    foreach_particle_in (P)
      if (p().tag + 1 > n)
	n = p().tag + 1;
#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, &n, 1, MPI_UNSIGNED_LONG, MPI_MAX,
		   MPI_COMM_WORLD);
#endif
    p_ntags[P] = n;
  }
  return p_ntags[P];
}

/**
### Releases

`release_tracer_particles (P, pos, n)` releases `n` particles at the
positions `pos[0..n-1]` into the tracer-particle list `P`, tagged
`tp_ntags(P)`, ..., `tp_ntags(P) + n - 1` in that order. It is
collective, and every rank must pass the same positions in the same
order: each rank keeps only those it owns (`particle_is_mine()`), so a
release needs no communication and its tags do not depend on the
number of ranks. Positions outside the domain are dropped, leaving
their tags unused.

The particles join the list at the next `tracer_particles_step1`, the
start of a Runge–Kutta cycle: released at an even step, they start at
that time; at an odd step, one step later. Joining in the middle of a
cycle would advance them as if they had existed since its start. */

typedef struct {
  Particles P;
  particle q;
} tp_release;

static tp_release * tp_released = NULL;
static long tp_nreleased = 0, tp_areleased = 0;

void release_tracer_particles (Particles P, coord * pos, long unsigned int n)
{
  long unsigned int first = tp_ntags (P);
  for (long unsigned int k = 0; k < n; k++)
    if (particle_is_mine (pos[k])) {
      if (tp_nreleased == tp_areleased) {
	tp_areleased = tp_areleased ? 2*tp_areleased : 64;
	tp_released = realloc (tp_released, tp_areleased*sizeof(tp_release));
      }
      particle q = {0};
      q.x = pos[k].x, q.y = pos[k].y, q.z = pos[k].z;
      q.tag = first + k;
      tp_released[tp_nreleased].P = P;
      tp_released[tp_nreleased++].q = q;
    }
  p_ntags[P] = first + n;
}

void tp_insert_releases (void)
{
  for (long k = 0; k < tp_nreleased; k++) {
    Particles P = tp_released[k].P;
    change_plist_size (P, 1);
    pl[P][pn[P] - 1] = tp_released[k].q;
  }
  tp_nreleased = 0;
}

/**
## Output sorted by tag

`output_tracer_particles()` appends the positions of the list `P` at
`time` to the binary file `fname`, as one row per call with one column
per particle: the particle tagged `n` is always in column `n`,
whichever rank holds it. The file holds

~~~
long ntot, nfields = 3       header; ntot = tp_ntags(P), or `ntot`, at the first call
double time, xyz[ntot][3]    one row per call
~~~

so that, in Python, `np.fromfile (f, np.int64, 2)` reads the header,
then each row is `1 + 3*ntot` doubles. Missing particles (lost, not
released yet, or tags never used) are NaN, and tags of `ntot` and above
are ignored: with releases during the run, pass the final number of
tags as `ntot` on the first call.

With MPI, each rank writes its own particles in place through
collective MPI-IO: nothing is gathered on one rank and nothing is
sorted globally, so the cost per rank and the memory scale with the
particles it holds. The function is collective. With `PARTICLE_SURFACE`
the third coordinate is the surface height `eta` at the particle. */

typedef struct {
  long unsigned int tag;
  double c[3];
} tagged_position;

static int tagged_position_compare (const void * a, const void * b)
{
  long unsigned int ta = ((tagged_position *) a)->tag,
    tb = ((tagged_position *) b)->tag;
  return ta < tb ? -1 : ta > tb;
}

void output_tracer_particles (Particles P, char * fname = "particles.bin",
			      double time = t, long unsigned int ntot = 0)
{
  long unsigned int n = pn[P];
  if (!ntot)
    ntot = tp_ntags (P);
  tagged_position * r = malloc ((n + 1)*sizeof(tagged_position));
#if LAYERS && PARTICLE_SURFACE
  boundary ({eta});
#endif
  foreach_particle_in (P) {
    r[_j_particle] = (tagged_position){p().tag, {x, y, z}};
#if LAYERS && PARTICLE_SURFACE
    r[_j_particle].c[2] = particle_surface_height (x, y);
#endif
  }
  qsort (r, n, sizeof(tagged_position), tagged_position_compare);
  long head[2] = {0, 3};
#if _MPI
  particle_file fh;
  if (MPI_File_open (MPI_COMM_WORLD, fname, MPI_MODE_RDWR | MPI_MODE_CREATE,
		     MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
    fprintf (stderr, "output_tracer_particles(): cannot open '%s'\n", fname);
    exit (1);
  }
  MPI_Offset size = 0;
  if (pid() == 0) {
    MPI_File_get_size (fh, &size);
    if (size == 0) {
      head[0] = ntot;
      MPI_File_write_at (fh, 0, head, 2, MPI_LONG, MPI_STATUS_IGNORE);
    }
    else
      MPI_File_read_at (fh, 0, head, 2, MPI_LONG, MPI_STATUS_IGNORE);
  }
  MPI_Bcast (head, 2, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast (&size, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
  long N = head[0];
  MPI_Offset rowsize = sizeof(double)*(1 + 3*N),
    row = size ? (size - sizeof(head))/rowsize : 0,
    offset = sizeof(head) + row*rowsize;

  // Our particles with a valid tag, packed
  long unsigned int nl = 0;
  while (nl < n && r[nl].tag < N)
    nl++;
  double * buf = malloc ((3*nl + 1)*sizeof(double));
  MPI_Aint * displ = malloc ((nl + 1)*sizeof(MPI_Aint));
  for (long unsigned int j = 0; j < nl; j++) {
    memcpy (buf + 3*j, r[j].c, 3*sizeof(double));
    displ[j] = r[j].tag*3*sizeof(double);
  }
  long unsigned int nglobal = nl;
  MPI_Allreduce (MPI_IN_PLACE, &nglobal, 1, MPI_UNSIGNED_LONG, MPI_SUM,
		 MPI_COMM_WORLD);

  // Fill missing particles with NaN, each rank its share of the row
  if (nglobal < N) {
    long first = N*pid()/npe(), last = N*(pid() + 1)/npe();
    double * nan = malloc ((3*(last - first) + 1)*sizeof(double));
    for (long j = 0; j < 3*(last - first); j++)
      nan[j] = NAN;
    MPI_File_write_at_all (fh, offset + sizeof(double)*(1 + 3*first), nan,
			   3*(last - first), MPI_DOUBLE, MPI_STATUS_IGNORE);
    free (nan);
    MPI_File_sync (fh);
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_File_sync (fh);
  }
  if (pid() == 0)
    MPI_File_write_at (fh, offset, &time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

  // Each particle at its tag
  MPI_Datatype xyz, view = MPI_DOUBLE;
  MPI_Type_contiguous (3, MPI_DOUBLE, &xyz);
  MPI_Type_commit (&xyz);
  if (nl) {
    MPI_Type_create_hindexed_block (nl, 1, displ, xyz, &view);
    MPI_Type_commit (&view);
  }
  MPI_File_set_view (fh, offset + sizeof(double), MPI_DOUBLE, view, "native",
		     MPI_INFO_NULL);
  MPI_File_write_all (fh, buf, 3*nl, MPI_DOUBLE, MPI_STATUS_IGNORE);
  if (nl)
    MPI_Type_free (&view);
  MPI_Type_free (&xyz);
  MPI_File_close (&fh);
  free (buf), free (displ);
#else // !_MPI
  FILE * fp = fopen (fname, "r+b");
  if (fp)
    fread (head, sizeof(long), 2, fp);
  else {
    if (!(fp = fopen (fname, "w+b"))) {
      perror (fname);
      exit (1);
    }
    head[0] = ntot;
    fwrite (head, sizeof(long), 2, fp);
  }
  long N = head[0];
  double * row = malloc ((1 + 3*N)*sizeof(double));
  row[0] = time;
  for (long j = 0; j < 3*N; j++)
    row[1 + j] = NAN;
  for (long unsigned int j = 0; j < n && r[j].tag < N; j++)
    memcpy (row + 1 + 3*r[j].tag, r[j].c, 3*sizeof(double));
  fseek (fp, 0, SEEK_END);
  fwrite (row, sizeof(double), 1 + 3*N, fp);
  fclose (fp);
  free (row);
#endif // !_MPI
  free (r);
}
/**
We need a ugly hack to access `ADD_PART_MEMBERS`, as they only become
available after the #define macros are processed.
 */
@define part p 
void set_a_particle_attributes (particle * p) {
  double x = p->x, y = p->y, z = p->z;
  foreach_dimension() {
    part->locn.x = p->x;
    part->u.x = interpolate (u.x, x, y, z);
    part->u2.x = 0.;
  }
}

void set_particle_attributes (Particles p) {
  interpolate_vel (p, 1);
  foreach_particle_in(p) {
    foreach_dimension() {
      p().locn.x = p().x;
      //p().u2.x = p().u.x; //Do or dont?
    }
#if TP_VERTICAL
    p().locn.z = p().z;
#endif
  }
  if (!p_ntags[p]) // see "Tags" above
    tag_particles (p);
}

/**
### Tracer-particle (`tp`) initialization

Tracer-particles may be initialized by calling the following
functions.
 */

/**
`init_tp_cells()` puts one particle at the centre of each cell (of each
layer, with the multilayer solver). On the multigrid the particle of
global cell `(I, J)` (layer `k`) is tagged `(J*Nx + I)*nl + k`, which does
not depend on the decomposition and needs no communication. With
`PARTICLE_SURFACE` (surface tracers, see [particle.h]()), there is one
particle per cell, at the surface, tagged `J*Nx + I`. On a tree
(quadtree, octree) the cells are not a regular array, so the particles
are tagged in their local order, which depends on the decomposition.
Basilisk defines `MULTIGRID` for trees too, hence the test on `TREE`.

The global index of a cell is computed here from `point.i` and the
rank's place in the process grid (`mpi_coords`): Basilisk's own `_I`,
`_J`, `_K` are undefined at the end of `grid/multigrid.h` since
2026-10, so user code cannot use them. */

#if MULTIGRID && !TREE
# if _MPI
#  define TP_GLOBAL(d) (mpi_coords[d]*(1 << point.level))
# else
#  define TP_GLOBAL(d) 0
# endif
# define TP_I ((long unsigned int) (point.i - GHOSTS + TP_GLOBAL(0)))
# define TP_J ((long unsigned int) (point.j - GHOSTS + TP_GLOBAL(1)))
# define TP_K ((long unsigned int) (point.k - GHOSTS + TP_GLOBAL(2)))
#endif

Particles init_tp_cells(void) {
#if LAYERS && !PARTICLE_SURFACE
  int nlp = nl; // place_in_cells() puts one particle per layer
#else
  int nlp = 1;
#endif
  int np = 0;
  foreach(serial) // local cells: no reduction, which would sum over ranks
    np++;
  Particles p = new_tracer_particles (np*nlp);
  place_in_cells(p);
#if MULTIGRID && !TREE // the regular multigrid only: global cell indices
  particles pp = pl[p];
  long unsigned int Nx = N, ncells = N;
  int j = 0;
  for (int d = 1; d < dimension; d++)
    ncells *= N*(&Dimensions.x)[d]/Dimensions.x;
  foreach(serial, reduction(+:j)) {
#if dimension == 1
    long unsigned int g = TP_I;
#elif dimension == 2
    long unsigned int g = TP_J*Nx + TP_I;
#else
    long unsigned int g = (TP_K*(N*Dimensions.y/Dimensions.x) + TP_J)*Nx + TP_I;
#endif
    for (int k = 0; k < nlp; k++)
      pp[j++].tag = g*nlp + k;
  }
  p_ntags[p] = ncells*nlp;
#else // trees
  tag_particles (p); // depends on the decomposition
#endif
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}

/**
`init_tp_square()` and `init_tp_circle()` place all the particles on
rank 0, so they are tagged there, in placement order, before being
distributed. */

Particles init_tp_square (int n = 10,
			  double xm = 0, double ym = 0,
			  double l = L0, double zp = 0)
{
  Particles p;
  if (pid() == 0) {
    p = new_tracer_particles (sq(n));
    place_in_square (p, n, xm, ym, l, zp);
  } else { 
    p = new_tracer_particles (0);
  }
  tag_particles (p);
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}

Particles init_tp_circle (int n = 100,
			  double xm = 0, double ym = 0,
			  double l = L0/2., double zp = 0)
{
  Particles p;
  if (pid() == 0) {
    p = new_tracer_particles (n);
    place_in_circle (p, n, xm, ym, l, zp);
  }
  else
    p = new_tracer_particles (0);
  tag_particles (p);
  particle_boundary (p);
  set_particle_attributes (p);
  return p;
}

/**
## To do

* Implement functions for useful tracer-particle diagnostics
* Check if `tracer-particles` can exist alongside other `Particles`
* Box boundary conditions for particles (other than periodic)
* ~~~Tag particles (for MPI)~~~
* ~~~Bind tracer particles with vof interfaces~~~

## Tests

* [MPI test](tpp.c)  
* [Advection scheme test](tp18.c)

## Usage

* [A confetti cannon](confetti.c) 
* [The collision of two vortex rings](two_rings.c)

## See also

* [Visualize `Particles`](scatter2.h)
* [Tracer particles bound by interfaces](vof-tracer-particles.h)
*/