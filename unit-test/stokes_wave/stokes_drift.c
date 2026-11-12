/**
# Forced pure-Stokes-drift reference (multilayer solver)

A reference solution for the **wave-induced drift of a single periodic
deep-water wave**, *with wave forcing but no current and no rotation*.
This is the forced counterpart of the unforced reference in `../unforced/`:
the wave amplitude (and hence the Stokes drift) is held constant against
viscous decay, so that the Lagrangian drift measured by the tracer
particles can be compared with a fixed-amplitude analytic Stokes drift
without the decay bias of the unforced run. Only the wave is forced
(see the `forcing` event); there is no current and no rotation.

The geometry, resolution and viscosity are identical to the unforced
reference and to the transverse-current study in `../../drift_transverse.c`,
so the runs differ only by the forcing (here) and the current (there).

Adapted from the 3D breaking Stokes-wave test case, with: the steepness
lowered to the ocean value $ak = 0.11$ (single-valued, non-breaking),
a flat bottom, the drift diagnostics, Lagrangian particles, and the
wave forcing. */

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
// #include "scatter2.h"

#if !_GPU
# include "view.h"
#endif

/**
Wave parameters (identical to the main case): steepness `ak`,
wavenumber `k_`, depth `h_`, gravity `g_`, and Reynolds number
$Re = c\lambda/\nu$. */

double ak = 0.11;
double RE = 40000.;
double k_ = 2.*pi, h_ = 0.5, g_ = 1.;

#define T0  (k_*L0/sqrt(g_*k_))
#define TEND (30.*T0)

/**
## Wave forcing (this is the *forced* wave-only case)

Unlike the unforced reference, energy is injected into the wave's
*fluctuating* (orbital) motion so that its amplitude --- and hence the
Stokes drift --- stays constant against viscous decay, instead of
decaying through the run. The mean/drift part of `u.x` is never forced,
so the drift we measure evolves freely. A deficit-restoring feedback
holds the wave kinetic energy at its initial value `KE0f` over `tau_w`.
There is no current and no rotation here, so (unlike the main case)
only the wave is forced. Set `gammaw_on = false` to recover the
unforced case. See the `forcing` event below. */

bool   gammaw_on = true;
double gammaw    = 0.;    // set in main() to 2*nu*k^2 (baseline decay rate)
double tau_w     = 0.;    // set in main() to T0 (energy restoring time)
double KE0f      = -1.;   // wave-energy set point (captured on first call)

#include "stokes.h"

#include "tracer-particles.h"
#include "write_class_particles.h"
Particles prt;
TextWriter pw;

/**
Analytic deep-water Stokes drift $\uS(z) = \omega a^2 k\, e^{2kz}$,
used as the reference the simulation is compared against ($kh=\pi$, so
the deep-water form is accurate). */

double uS_analytic (double z)
{
  double omega = sqrt(g_*k_*tanh(k_*h_));
  double a = ak/k_;
  return omega*a*a*k_*exp(2.*k_*z);
}

int main()
{
  // size (2.)
  origin (-L0/2., -L0/2.);

  /**
  Doubly periodic and flat-bottomed, matching the main case exactly so
  the two runs differ only by the presence of the current/forcing. */

  periodic (right);
  periodic (top);

  N  = 256;
  nl = 20;
  G  = g_;
  nu = 1./RE;

  gammaw = 2.*nu*sq(k_);   // baseline wave-forcing rate ~ viscous decay
  tau_w  = T0;             // energy restoring timescale

  theta_H = 0.505;   // implicit damping of fast modes (as in breaking.c)

  run();
}

/**
## Initial conditions

Third-order Stokes wave on a flat bottom; no current
($v\equiv 0$). Layers refined toward the surface. */

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
      u.x[] = u_x(x, z);     // wave orbital velocity (x)
      u.y[] = 0.;            // no current
      w[]   = u_y(x, z);     // wave orbital velocity (vertical)
      z += h[]/2.;
    }
  }
  /**
  ## Lagrangian floats along a vertical column

  We seed `nl` particles in a single vertical column at the reference
  location $(x,y)=(0,0)$, one at each layer mid-depth, so their
  horizontal drift measures the wave-induced drift as a function of
  depth. The mid-depths are reconstructed from the geometric layer
  distribution `beta[]` (which `geometric_beta()` has just set),
  evaluated at the reference column where the total depth is
  $H_0 = \eta(0) - z_b = \mathrm{wave}(0,0) + h$.

  Note we cannot use `h[]` here: it is a grid field requiring a cell
  (Point) context, which the particle iterator does not provide. We
  index the particle array directly (`pl[column][l]`). */

  /**
  MPI-safe creation (cf. `init_tp_square` in tracer-particles.h):
  the `nl` particles are created and placed on rank 0 only --- creating
  `nl` on *every* rank would give `npe * nl` duplicates. The reference
  column is offset from the domain centre to avoid sitting on a
  process boundary; `particle_boundary()` then migrates each particle
  to its owning rank, and `set_particle_attributes()` assigns the
  global tags that the MPI writer sorts by (it must be called on all
  ranks). */

  prt = new_tracer_particles (pid() == 0 ? nl : 0);
  if (pid() == 0) {
    particle * pp = pl[prt];
    double xref = - L0/4., yref = - L0/4.;   // interior reference column
    double H0 = wave (xref, 0.) + h_;        // local total depth
    double zc = - h_;                        // start at the bottom
    for (int l = 0; l < nl; l++) {
      double hl = H0*beta[l];
      zc += hl/2.;                           // mid-depth of layer l
      pp[l].x = xref;
      pp[l].y = yref;
      pp[l].z = zc;
      zc += hl/2.;
    }
  }
  particle_boundary (prt);       // locate/migrate particles to owning ranks
  set_particle_attributes (prt); // assign global tags (collective)

  pw = init_particle_files ("w", prt, "tracer");
}

/**
## Viscous terms (as in breaking.c) */

event viscous_term (i++)
{
  foreach()
    vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.);
  horizontal_diffusion ({u, w}, nu, dt);
}

/**
## Wave forcing

Applied after viscosity so it directly counters the decay that step.
Energy is injected into the fluctuating part of the motion only
(`u.x - <u.x>_l` and `w`); the horizontal (wave-averaged) mean `mx[l]`
is left untouched so the drift evolves freely. The rate `g` restores
the wave kinetic energy `kew` to its set point `KE0f` over `tau_w`. */

event forcing (i++)
{
  if (!gammaw_on)
    return 0;

  /* per-layer horizontal mean of u.x (the drift) */
  double mx[nl];
  for (int l = 0; l < nl; l++) mx[l] = 0.;
  foreach (reduction(+:mx[:nl]))
    foreach_layer()
      mx[point.l] += u.x[]*sq(Delta);
  for (int l = 0; l < nl; l++) mx[l] /= sq(L0);

  /* wave kinetic energy (fluctuation only) for the feedback */
  double kew = 0.;
  foreach (reduction(+:kew))
    foreach_layer() {
      double up = u.x[] - mx[point.l];
      kew += (sq(up) + sq(w[]))*h[]*dv();
    }
  if (KE0f < 0.) KE0f = kew;

  double g = 0.;
  if (kew > 0.)
    g = gammaw + (KE0f - kew)/(2.*kew*tau_w);

  foreach()
    foreach_layer() {
      u.x[] += dt*g*(u.x[] - mx[point.l]);   // sustain wave, not drift
      w[]   += dt*g*w[];
    }
}

/**
## Diagnostics

### Energy */

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	      norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
    }
    gpe += sq(eta[])*dv();
  }
  fprintf (stderr, "%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}

/**
### Drift profile and Stokes transport

Once per wave period we write the horizontally averaged along-wave
velocity $\langle u\rangle_l(z)$ to `profiles.dat`, alongside the
analytic Stokes drift $\uS(z)$ for comparison. Columns:
`t/T0  layer  z_mean  <u_x>  uS_analytic`.

We also write the depth-integrated Stokes transport
$Q = \int \langle u\rangle\,\mathrm{d}z$ to `transport.dat`, compared
with the analytic value $\uS$-integral $= \omega a^2/2$. Because the
depth-integrated mean flux is the (conserved) wave momentum, $Q$ is the
most robust single measure of the wave-induced transport in this
finite-depth periodic domain, and is the recommended reference number.
Columns: `t/T0  Q_sim  Q_analytic`. */

event profiles (t = 0; t += T0/10.)
{
  static FILE * fp = NULL;
  if (!fp) { fp = fopen ("profiles.dat", "w");
    fprintf (fp, "# t/T0 layer z_mean <u_x> uS_analytic\n"); }
  static FILE * fq = NULL;
  if (!fq) { fq = fopen ("transport.dat", "w");
    fprintf (fq, "# t/T0 Q_sim Q_analytic\n"); }

  double area = L0*L0;
  double Qsim = 0.;
  for (int l = 0; l < nl; l++) {
    double um = 0., zm = 0., hm = 0.;
    foreach (reduction(+:um) reduction(+:zm) reduction(+:hm)) {
      double z = zb[];
      for (int j = 0; j < l; j++)
	z += h[0,0,j];
      z += h[0,0,l]/2.;
      um += u.x[0,0,l]*sq(Delta);
      zm += z*sq(Delta);
      hm += h[0,0,l]*sq(Delta);
    }
    um /= area; zm /= area; hm /= area;   // <u_x>, mean z, mean layer thickness
    fprintf (fp, "%g %d %g %g %g\n", t/T0, l, zm, um, uS_analytic (zm));
    Qsim += um*hm;                        // depth-integrated transport
  }
  fprintf (fp, "\n"); fflush (fp);

  double omega = sqrt(g_*k_*tanh(k_*h_)), a = ak/k_;
  double Qana = omega*a*a/2.;             // deep-water Stokes transport
  fprintf (fq, "%g %g %g\n", t/T0, Qsim, Qana);
  fflush (fq);
}

event writes (t += T0/10.) {
  particles_kinematics (pw, t);
}
/**
### Free-surface movie (top-layer x-velocity) */

event movie (t += T0/10.)
{
  char utopx[20]; sprintf (utopx, "u%d.x", nl - 1);
#if !_GPU
  view (fov = 17.3106, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768);
  char s[80]; sprintf (s, "t = %.1f T0", t/T0);
  draw_string (s, size = 80);
  // for (double x = -1; x <= 1; x++)
  //   translate (x)
      squares (utopx, linear = true, z = "eta", min = -0.15, max = 0.6);
  box ();
  save ("movie.mp4");
  save ("movie.png");
  if (t == 0)
    save ("movie0.png");
#else
  char utop[20]; sprintf (utop, "u%d", nl - 1);
  vector u = lookup_vector (utop);
  output_ppm (u.x, file = "movie.mp4", linear = true, min = -0.15, max = 0.6, n = 512);
#endif
}

event end (t = TEND) {}
