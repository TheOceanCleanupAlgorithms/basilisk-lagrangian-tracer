/* Multilayer smoke test: particles advected by layered/hydro.h,
   serial vs MPI on the multigrid. With -DCELLS=1, one particle per cell
   and layer (init_tp_cells()) instead of a square, written with its
   tag as id. With -DSHEAR=1 (and -DCELLS=1), flat layers each moving at
   its own uniform velocity U(k), V(k): each particle must move with the
   velocity of its own layer.

   With -DPARTICLE_SURFACE=1 (surface tracers), one particle per cell.
   With SHEAR, half of them are then moved to z = -0.9 (bottom layer) and
   half to z = 5 (above the surface): all must move with the top layer,
   and z must not change (no reflection). With -DWRITER=1 instead, the
   positions at t = 0 are written by write_class_particles.h, and each
   rank writes eta-<pid> (tag, eta in the particle's cell), to check
   that the z column is the surface height. */
#include "grid/multigrid.h"
#include "layered/hydro.h"
#define ADD_PART_MEM coord u; coord u2; coord locn; long unsigned int tag; long id;
#include "tracer-particles.h"
#if WRITER
# include "write_class_particles.h"
#endif

Particles P;
int NP = 30;

int main() {
  L0 = 10.; X0 = Y0 = -5.;
  G = 9.81;
  N = 64;
  nl = 3;
  periodic (right);
  periodic (top);
  run();
}

#define U(k) (0.2*((k) + 1))
#define V(k) (0.1*(nl - (k)))
#define Z0(tag) ((tag) % 2 ? 5. : -0.9) // above the surface, bottom layer

event init (i = 0) {
  foreach() {
    zb[] = -1.;
#if SHEAR
    int k = 0;
    foreach_layer()
      h[] = 1./nl, u.x[] = U(k), u.y[] = V(k), k++;
#else
    foreach_layer()
      h[] = (1. + 0.2*exp(-(sq(x - 1.) + sq(y))))/nl, u.x[] = 0.8, u.y[] = 0.4;
#endif
  }
#if CELLS
  P = init_tp_cells();
  foreach_particle_in (P)
    p().id = p().tag;
#if PARTICLE_SURFACE && SHEAR
  foreach_particle_in (P)
    p().z = Z0(p().tag);
#endif
#else
  P = init_tp_square (n = NP, l = 8., zp = -0.5);
  double l = 8., d = l/NP;
  foreach_particle_in (P)
    p().id = lround ((x + l/2. - d/2.)/d)*NP + lround ((y + l/2. - d/2.)/d);
#endif
}

event end (t = 3) {
  char name[80];
  sprintf (name, "out-%d", pid());
  FILE * fp = fopen (name, "w");
  foreach_particle_in (P)
    fprintf (fp, "%ld %.15g %.15g\n", p().id, x, y);
  fclose (fp);
  long unsigned int n = pn[P];
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &n, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (pid() == 0)
    fprintf (stderr, "i=%d t=%g particles=%lu\n", i, t, n);
#if SHEAR
  // largest distance to the exact position, x0 + U(k) t, y0 + V(k) t
  double err = 0.;
  foreach_particle_in (P) {
#if PARTICLE_SURFACE
    long g = p().tag, k = nl - 1; // one particle per cell, top layer
    err = max (err, fabs (z - Z0(p().tag))); // z unchanged
#else
    long g = p().tag/nl, k = p().tag % nl;
#endif
    double ex = X0 + (g % N + 0.5)*L0/N + U(k)*t;
    double ey = Y0 + (g / N + 0.5)*L0/N + V(k)*t;
    double dx = remainder (x - ex, L0), dy = remainder (y - ey, L0);
    err = max (err, max (fabs (dx), fabs (dy)));
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
#if PARTICLE_SURFACE
  bool good = err < 1e-9 && n == (long unsigned int) N*N;
#else
  bool good = err < 1e-9;
#endif
  if (pid() == 0)
    printf ("shear %s %.3g\n", good ? "ok" : "wrong", err);
#endif
}

#if WRITER
event write_init (i = 0, last) { // after step 1: positions at t = 0
  TextWriter w = init_particle_files ("w", P, "w");
  particles_kinematics (w, t);
  close_particle_files (&w);
  char name[80];
  sprintf (name, "eta-%d", pid());
  FILE * fp = fopen (name, "w");
  foreach (serial) // the tag of init_tp_cells() with PARTICLE_SURFACE
    fprintf (fp, "%ld %13.6e\n", (long) (TP_J*N + TP_I), eta[]);
  fclose (fp);
}
#endif
