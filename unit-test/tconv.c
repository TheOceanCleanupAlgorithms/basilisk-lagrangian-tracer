/**
# Order of convergence in time of the tracer particles

Solid-body rotation about the origin with a time-dependent rate
$\omega(t) = 1 + \cos(2t)$:
$$u = -\omega(t)\,y, \qquad v = \omega(t)\,x$$
The flow is linear in space, so the linear interpolation is exact and
only the time error of the Runge–Kutta scheme is left. The exact
solution is a rotation by $\Theta(t) = t + \sin(2t)/2$.

Usage: `./tconv dt=... [r=1] [lag=0] [rk2=0] [nstep=...]`

* `r`: odd steps are `r*dt` (unequal steps, $\alpha = 1/(1+r)$);
* `solver=1`: like a solver (`centered.h`, `layered/hydro.h`), a `last`
  event defined before the particle code advances `u` from $t$ to
  $t + \Delta t$, so the particle `last` events see $u(t+\Delta t)$;
* `lag=1`: $u(t+\Delta t)$ is set in the ordinary event `velocity`
  (wrong on purpose: the particles then read a field one step ahead);
* `rk2=1`: `P_RK2`.

It prints `dt t err_t err_tdt err_shift`: the largest position error at
the final time against the exact solution at $t$, at $t + \Delta t$, and
at $t + \Delta t$ for a particle released at $t_1 = \Delta t_0$ instead
of 0. */

#include "run.h"
vector u[];
int SOLVER = 0;
void set_u (double tu);

event advance (i++, last) { // the solver's update, before the particle events
  if (SOLVER)
    set_u (t + dt);
}

#include "tracer-particles.h"

double R = 1., TEND = 2.;
int LAG = 0, NSTEP = 0;
double DT1 = 0.; // the first time step
Particles P;
coord * x0 = NULL;
long np = 0;

double rot_rate (double t) { return 1. + cos (2.*t); }
double rot_angle (double t) { return t + sin (2.*t)/2.; }

void set_u (double tu) {
  double om = rot_rate (tu);
  foreach() {
    u.x[] = - om*y;
    u.y[] = om*x;
  }
}

event set_dtmax (i++) { dt = dtnext (i % 2 ? R*DT : DT); }

event velocity (i++) {
  if (i == 0) DT1 = dt;
  if (!SOLVER)
    set_u (LAG ? t + dt : t);
}

event init (i = 0) {
  set_u (LAG ? t + dt : t);
  P = init_tp_square (n = 5, l = 0.4);
  np = pn[P];
  x0 = malloc (np*sizeof(coord));
  foreach_particle_in(P)
    x0[p().tag] = (coord){x, y};
}


double err_at (double te, double t0) {
  double c = cos (rot_angle (te) - rot_angle (t0)), s = sin (rot_angle (te) - rot_angle (t0));
  double e = 0.;
  foreach_particle_in(P) {
    coord e0 = x0[p().tag];
    double xe = c*e0.x - s*e0.y, ye = s*e0.x + c*e0.y;
    double d = sqrt (sq(p().locn.x - xe) + sq(p().locn.y - ye));
    if (d > e) e = d;
  }
  return e;
}

/* last, after step 1: p().locn is the position at the end of the
   last full Runge-Kutta cycle */
event stop (i = NSTEP, last) {
  printf ("%g %.12g %.6e %.6e %.6e\n", DT, t, err_at (t, 0), err_at (t + dt, 0),
	  err_at (t + dt, DT1));
  free (x0);
}

int main (int argc, char * argv[]) {
  int rk2 = 0;
  DT = 0.01;
  for (int a = 1; a < argc; a++) {
    char * v = strchr (argv[a], '=');
    if (!v) return 1;
    *v++ = '\0';
    if (!strcmp (argv[a], "dt")) DT = atof (v);
    else if (!strcmp (argv[a], "r")) R = atof (v);
    else if (!strcmp (argv[a], "lag")) LAG = atoi (v);
    else if (!strcmp (argv[a], "solver")) SOLVER = atoi (v);
    else if (!strcmp (argv[a], "rk2")) rk2 = atoi (v);
    else if (!strcmp (argv[a], "nstep")) NSTEP = atoi (v);
    else return 1;
  }
  P_RK2 = rk2;
  if (!NSTEP) // an even number of steps, about TEND
    NSTEP = 2*(int) round (TEND/((1. + R)*DT));
  origin (-0.5, -0.5);
  N = 32;
  run();
}
