/**
# Particles in a Stokes wave: vertical motion with `layered/nh.h`

A reduced copy of the user's reference case
[stokes_wave/stokes_drift.c](): the same third-order Stokes wave
(`ak = 0.11`, deep water), non-hydrostatic multilayer solver with
remapping, and a column of particles at $(x, y) = (-L_0/4, -L_0/4)$, one
at the middle of each layer. Smaller and shorter (`N = 64`, 10 layers,
3 wave periods) and without the wave forcing, so that it runs in
seconds.

Particles must move vertically with `w`. `run-tests.sh` checks that,
with the prescribed flow (below), the trajectories are the same on 4
and 16 ranks as in serial. On 16 ranks the column lies exactly on a
corner of four subdomains.

Each rank writes `traj-<pid>`: `t tag x y z` every tenth of a period,
and rank 0 prints `fluid t eta w` at the column (`w` of the top layer),
to compare the fluid itself between serial and MPI runs. */

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

double ak = 0.11, RE = 40000.;
double k_ = 2.*pi, h_ = 0.5, g_ = 1.;
#define T0  (k_*L0/sqrt(g_*k_))

#include "stokes.h"

/**
With `-DPRESCRIBED=1` the flow is not the solver's: the layers, `u`
and `w` are set from the same Stokes wave, moving at its phase speed, at
time `t`, both before the solver's events (an ordinary event, read by
the particles) and after them (a `last` event defined before including
the particle code, which undoes the solver's update). The flow is then exactly the same at
any number of ranks, unlike the solved one (the pressure solver of
`nh.h` is not bit-identical across decompositions), so the trajectories
must be bit-identical between serial and MPI runs. */

#if PRESCRIBED
void prescribe_flow (void)
{
  double c = sqrt (g_*tanh (k_*h_)/k_), xs = c*t; // phase speed
  foreach() {
    double H = wave (x - xs, 0) + h_, z = - h_;
    eta[] = zb[] + H; // read by the reflection at the surface
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x (x - xs, z), u.y[] = 0.;
      w[] = u_y (x - xs, z);
      z += h[]/2.;
    }
  }
}

event prescribe0 (i++) prescribe_flow();
event prescribe (i++, last) prescribe_flow();
#endif

#include "tracer-particles.h"

Particles prt;

int main()
{
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N  = 64;
  nl = 10;
  G  = g_;
  nu = 1./RE;
  theta_H = 0.505;
  run();
}

event init (i = 0)
{
  CFL = 0.8;
  geometric_beta (1./3., true);
  foreach() {
    zb[] = - h_;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x(x, z);
      u.y[] = 0.;
      w[]   = u_y(x, z);
      z += h[]/2.;
    }
  }
  prt = new_tracer_particles (pid() == 0 ? nl : 0);
  if (pid() == 0) {
    particle * pp = pl[prt];
    double xref = - L0/4., yref = - L0/4.;
    double H0 = wave (xref, 0.) + h_, zc = - h_;
    for (int l = 0; l < nl; l++) {
      double hl = H0*beta[l];
      zc += hl/2.;
      pp[l].x = xref, pp[l].y = yref, pp[l].z = zc;
      zc += hl/2.;
    }
  }
  particle_boundary (prt);
  set_particle_attributes (prt);
  char name[80];
  sprintf (name, "traj-%d", pid());
  fclose (fopen (name, "w"));
}

event viscous_term (i++)
{
  foreach()
    vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.);
  horizontal_diffusion ({u, w}, nu, dt);
}

event snap (t += T0/10., last) // after the particle events
{
  char name[80];
  sprintf (name, "traj-%d", pid());
  FILE * fp = fopen (name, "a");
  foreach_particle_in (prt)
    fprintf (fp, "%.10g %ld %.15g %.15g %.15g\n", t, (long) p().tag, x, y, z);
  fclose (fp);

  /* the fluid at the column, to tell differences of the fluid solver
     (serial against MPI) from differences of the particle code */
  double e = interpolate (eta, - L0/4., - L0/4.), wtop = nodata;
  foreach_point (- L0/4. + 1e-9, - L0/4. + 1e-9, reduction(min:wtop)) {
    point.l = nl - 1;
    wtop = w[];
  }
  if (pid() == 0)
    printf ("fluid %.10g %.15g %.15g\n", t, e, wtop);
}

event end (t = 3.*T0) {}
