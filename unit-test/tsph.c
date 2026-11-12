/**
# Tracer particles in degrees on the sphere

`PARTICLE_GEOGRAPHY` with [spherical.h](/src/spherical.h), on the
North Pacific domain: closed, `dimensions(2,1)`, longitude [105:255]
(the domain's convention, across the date line), latitude [-10:65],
`Radius` in metres. The velocity is uniform, so interpolation is exact
and so is the trajectory:

* `flow=zonal`: u = U, v = 0. The latitude stays constant and the
  longitude is λ0 + U t/(R cos φ0) 180/π, a different rate at each
  latitude; particles starting east of 180 must stay in [105:255].
* `flow=meridional`: u = 0, v = V. The longitude stays constant and the
  latitude is φ0 + V t/R 180/π.

Particles start at least 8 cells from the walls, where the ghost values
of `u` are not uniform. Usage: `./tsph [key=value ...]` with `flow`
(`zonal`), `n` (256 cells in longitude), `dt` (4000 s), `tend` (1e6 s).

Rank 0 prints

~~~
geography <flow> <ok|wrong> <largest error in degrees> <particles> <wrapped>
~~~

with `wrapped` the number of particles starting east of 180. Each rank
writes `out-<pid>`: `tag lon lat` at the end. With `-DWRITER=1` (and
`-DWCP_LON_WRAPPED=1`), `particles_kinematics()` writes the initial
positions to `w_particle_position_x.dat`. */

#include "grid/multigrid.h"
#include "spherical.h"
#include "run.h"
vector u[];

#define PARTICLE_GEOGRAPHY 1
#include "tracer-particles.h"
#if WRITER
# include "write_class_particles.h"
#endif

Particles P;
int NX = 18, NY = 11;           // particles in longitude and latitude
double U = 1., V = 1., TEND = 1e6;
bool ZONAL = true;

#define LON0(k) (115. + 85.*((k)/NY + 0.5)/NX)
#define LAT0(k) ((ZONAL ? 0. : -5.) + (ZONAL ? 50. : 55.)*((k)%NY + 0.5)/NY)

int main (int argc, char * argv[]) {
  DT = 4000.; // divides tend: see the note at set_dtmax
  N = 256;
  for (int a = 1; a < argc; a++) {
    char * v = strchr (argv[a], '=');
    if (!v) {
      fprintf (stderr, "usage: %s [key=value ...]\n", argv[0]);
      return 1;
    }
    *v++ = '\0';
    if (!strcmp (argv[a], "flow")) ZONAL = !strcmp (v, "zonal");
    else if (!strcmp (argv[a], "n")) N = atoi (v);
    else if (!strcmp (argv[a], "dt")) DT = atof (v);
    else if (!strcmp (argv[a], "tend")) TEND = atof (v);
    else {
      fprintf (stderr, "unknown key '%s'\n", argv[a]);
      return 1;
    }
  }
  dimensions (2, 1);
  size (150);
  origin (105, -10);
  Radius = 6371220.;
  run();
}

event init (t = 0) {
  foreach()
    u.x[] = ZONAL ? U : 0., u.y[] = ZONAL ? 0. : V;
  // placed on rank 0 and tagged in placement order: tag k is particle k
  if (pid() == 0) {
    P = new_tracer_particles (NX*NY);
    foreach_particle_in (P)
      p().x = LON0(_j_particle), p().y = LAT0(_j_particle), p().z = 0.;
  }
  else
    P = new_tracer_particles (0);
  tag_particles (P);
  particle_boundary (P);
  set_particle_attributes (P);
}

/**
`dtnext()` only knows the events before the one calling it, and
`set_dtmax` was declared in [tracer-particles.h](), before `end`. If
`dt` did not divide `tend`, the event system would shorten the last
step to reach `tend` without changing `dt`, and the particles would
move by `dt` regardless. The solvers call `dtnext()` in a `last` event,
after all others, and do not have this problem. */

event set_dtmax (i++) { dt = dtnext (DT); }
event velocity (i++) {}

#if WRITER
event write_init (i = 0, last) {
  TextWriter w = init_particle_files ("w", P, "w");
  particles_kinematics (w, t);
  close_particle_files (&w);
}
#endif

event end (t = TEND) {
  char name[80];
  sprintf (name, "out-%d", pid());
  FILE * fp = fopen (name, "w");
  double err = 0.;
  long n = 0, wrapped = 0;
  foreach_particle_in (P) {
    long k = p().tag;
    double lon = LON0(k), lat = LAT0(k);
    if (ZONAL)
      lon += U*t/(Radius*cos (lat*pi/180.))*180./pi;
    else
      lat += V*t/Radius*180./pi;
    err = max (err, max (fabs (x - lon), fabs (y - lat)));
    n++, wrapped += LON0(k) >= 180.;
    fprintf (fp, "%ld %.15g %.15g\n", k, x, y);
  }
  fclose (fp);
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &n, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &wrapped, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (pid() == 0)
    printf ("geography %s %s %.3g %ld %ld\n", ZONAL ? "zonal" : "meridional",
	    err < 1e-9 && n == NX*NY ? "ok" : "wrong", err, n, wrapped);
}
