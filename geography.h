/**
# Geographic conversions

Functions only: including this file changes nothing else, so it can be
used in any simulation. Tracer particles in degrees are switched on
separately, with `PARTICLE_GEOGRAPHY` (see [tracer-particles.h]()).

## Metres and degrees on the Earth's ellipsoid

From the American Practical Navigator (Vol. II, 1975 edition, p. 5).
`y` is the latitude in degrees at which the conversion applies: for long
displacements, use the midpoint latitude.

These are for real distances, such as observations. The model grid of
[spherical.h](/src/spherical.h) is a sphere of radius `Radius`, with
`Radius*pi/180` metres per degree everywhere: within 0.6 % of these,
but not the same, so tracer particles use the sphere (`tp_rate()` in
[tracer-particles.h]()). */

#ifndef GEOGRAPHY_H
#define GEOGRAPHY_H

double meters_to_degrees_lon (double dx_meters, double y) {
  double rlat = y*M_PI/180;
  if (fabs (cos (rlat)) < 1e-12) // at a pole (cos(pi/2) is 6e-17, not 0)
    return 0;
  return dx_meters/(111415.13*cos (rlat) - 94.55*cos (3*rlat));
}

double meters_to_degrees_lat (double dy_meters, double y) {
  double rlat = y*M_PI/180;
  return dy_meters/(111132.09 - 556.05*cos (2*rlat) + 1.2*cos (4*rlat));
}

double degrees_lon_to_meters (double dx, double y) {
  double rlat = y*M_PI/180;
  return dx*(111415.13*cos (rlat) - 94.55*cos (3*rlat));
}

double degrees_lat_to_meters (double dy, double y) {
  double rlat = y*M_PI/180;
  return dy*(111132.09 - 556.05*cos (2*rlat) + 1.2*cos (4*rlat));
}

/**
## Longitude conventions

A domain that straddles the date line is written in a convention such
as [105:255] (`origin()`, `X0`, `L0`), and particle positions must stay
in it: `locate()` and the MPI decomposition use it. Wrapped to
[-180:180), 205 becomes -155, which lies outside such a domain, and a
particle moved there is dropped as lost. So:

* `lon_to_domain()` maps any longitude into [`X0`:`X0` + 360), the
  convention of positions: use it for positions given in another
  convention, e.g. before `release_tracer_particles()`;
* `lon_to_wrapped()` maps into [-180:180) (180 gives -180): use it on
  output only (see `WCP_LON_WRAPPED` in [write_class_particles.h]()),
  never on `p().x`.

`fmod(a, b)` is `a - b*trunc(a/b)`, wrong for negative `a` here, hence
`a - b*floor(a/b)`. */

double constrain_longitude_to_valid_domain (double lon) {
  return ((lon + 180) - 360*floor ((lon + 180)/360)) - 180;
}

double lon_to_wrapped (double lon) {
  return constrain_longitude_to_valid_domain (lon);
}

double lon_to_domain (double lon) {
  return X0 + ((lon - X0) - 360*floor ((lon - X0)/360));
}

#endif // GEOGRAPHY_H
