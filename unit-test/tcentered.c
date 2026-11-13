/**
# Particle time with `navier-stokes/centered.h`

A uniform flow accelerated by `G` in a periodic box: $u(t) = G\,t$
exactly, so a particle moves to $x_0 + G t^2/2$, which the Runge–Kutta
scheme gives exactly (third order is exact for a velocity linear in
time).

`centered.h` advances `u` in its own `last` events, before the
particles' `last` events. The particles must still use $u(t)$, not
$u(t + \Delta t)$: with the old code, which read `u` in `last` events,
the error was $G\,t\,\Delta t$ (positions one step late).

Usage: `./tcentered [dt]`. It prints `centered dt t err` at $t = 1$,
with `err` the position error at time `t`. The steps grow from a small
first step to `dt` (`dtnext()`), so they are unequal at first. */

#include "navier-stokes/centered.h"
#include "tracer-particles.h"

double G = 0.05, XP = -0.1;
Particles P;
face vector av[];

event init (i = 0) {
  a = av;
  P = init_tp_square (n = 1, xm = XP, ym = 0., l = 0.1);
}

event acceleration (i++) {
  foreach_face(x)
    av.x[] += G;
}

/* Basilisk sets `t` to exactly 1 at the event, a few 1e-6 away from
   the sum of the steps: the particles are at the sum */
double tsum = 0.;

event stop (t = 1, last) {
  foreach_particle_in(P)
    printf ("centered %g %.12g %.3e\n", dt, tsum,
	    p().locn.x - (XP + G*sq(tsum)/2.));
}

event elapsed (i++, last) { tsum += dt; }

int main (int argc, char * argv[]) {
  DT = argc > 1 ? atof (argv[1]) : 0.01;
  origin (-0.5, -0.5);
  periodic (right);
  periodic (top);
  N = 16;
  run();
}
