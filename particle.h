/**
# A Particle class

This headerfile defines some useful functions and types for generic
particles in terms of initialization and cleanup as well as iterators. 
It then includes the relevant sub-library based on solver type.

A particle is defined by a 3D position and possible additional members
that may be `#define`d via the hook (mass, size etc.).
*/
typedef struct {
  double x;
  double y;
  double z;
  bool inertia;
#ifdef ADD_PART_MEM
  ADD_PART_MEM
#endif
} particle;


/**
In practice, more than one particle is stored in a list of
`particle`s. The code maintains a list of all lists (`pl`) together
with a related `terminate_int`-terminated array (`pn`) that stores the
numbers of particles in each list and the (dynammically) allocated
space. Furthermore, a list can be referenced by an integer number ,
`Particles` (mind the case), indexing the how-manieth entry in `pn`
and `pl` a list of particles is.
 */

typedef particle * particles; //omit double pointer types and mistakes
typedef int Particles;        //An index reference
particles * pl = NULL;        
long unsigned int * pn = NULL, * pna =  NULL; // Particles and space
long unsigned int terminate_int = ULONG_MAX;
bool * p_global = NULL; // List needs a global MPI redistribution

/**
`p_ntags[l]` counts the identifiers (tags) issued for list `l`, the
same number on every rank. It is set by the tagging functions of
[tracer-particles.h](), carried through `pdump()` and `prestore()`, and
`P_NTAGS_UNKNOWN` after restoring an older dump. */

#define P_NTAGS_UNKNOWN ULONG_MAX
long unsigned int * p_ntags = NULL;

/**
`n` new particles can be declared using the `new_particles()`
function. All member-values are set to zero.
*/

Particles new_particles (long unsigned int n) {
  int j = 1;
  if (pn != NULL) 
    while (pn[j++] != terminate_int);
  pn = realloc (pn, (j + 1)*sizeof(long unsigned int));
  pna = realloc (pna, (j + 1)*sizeof(long unsigned int));
  pl = realloc (pl, j*sizeof(particle));
  p_global = realloc (p_global, j*sizeof(bool));
  p_ntags = realloc (p_ntags, j*sizeof(long unsigned int));
  pn[j] = terminate_int;
  pn[--j] = n;
  pna[j] = n;
  p_global[j] = true; // new particles may lie anywhere
  p_ntags[j] = 0;
  particle * pt = calloc (n, sizeof(particle));
  pl[j] = pt;
  return j;
}
/**
A list of particles can be added to a Particles list
*/

Particles * add_Particles_to_list (Particles p, Particles * plist) {
  int l = 0, t = 0;
  if (plist != NULL) {
    while (pn[l] != terminate_int) {
      if (l == plist[t]) 
	t++;
      l++;
    }
  }
  plist = realloc (plist, (t + 1)*sizeof(Particles));
  plist[t] = p;
  return plist;
}


/**
   A particle can be added to a list ...
*/
void change_plist_size (Particles Plist, int n) {
  long unsigned int n_new = n + pn[Plist];
  if (n_new > pna[Plist] || 2.5*(n_new + 1) < pna[Plist]) {
    pna[Plist] = 2*n_new;
    pl[Plist] = realloc (pl[Plist] , pna[Plist]*sizeof(particle));
  }
  pn[Plist] += n;
}

int add_particle (particle p, Particles Plist) {
#if _MPI
  if (!(locate (p.x, p.y, p.z).level > 0)) //Only place particle where it can be found
    return 0;
#endif
  
  change_plist_size (Plist, 1);
  pl[Plist][pn[Plist] - 1] = p; //append data
  return 1;
}
/**
... or removed (based on index or condition)
*/
void remove_particles_index (Particles Plist, int n, long unsigned int ind[n]) {
  int m = 0, j = 0;
  while (j < pn[Plist] - n) {
    while (m < n ? j + m == ind[m] : 0)
      m++;
    while (m < n ? j < pn[Plist] - n && j + m != ind[m] : j < pn[Plist] - n) {
      pl[Plist][j] = pl[Plist][j + m];
      j++;
    }
  }
  change_plist_size(Plist, -n);
}
/**
the "undo-everything" function can be used to prevent the
all-important memory leaks.
 */
void free_p (void) {
  int j = 0;
  while (pn[j++] != terminate_int) {
    free (pl[j - 1]);
    pl[j -1] = NULL;
  }
  free (pl);
  pl = NULL;
  free (pn);
  pn = NULL;
  free (pna);
  pna = NULL;
  free (p_global);
  p_global = NULL;
  free (p_ntags);
  p_ntags = NULL;
}

event free_particles (t = end)
  free_p();

/**
## Iterators

Two particle iterators are `@def`ined. Within these
`foreach_particle..` iterators, the coordinates `x, y, z` become
available and the `particle` data itself (`p()`).

* The `foreach_particle()` iterator loops over every particles.  

* The `foreach_particle_in(Particles P)` iterator loops over every
 particle in a single list referenced by `P`

Furthermore, the `foreach_P_in(Particles * Pl)` loops over lists of
`Particle`s. Inside the iterator, `Particles P` becomes available. 

For example, if you want to set the x-coordinate of all particles in
some list (say `tracer_particles`) to 1, you could do:

~~~literatec
...
Particles * tracer_particles;
...
foreach_P_in(tracer_particles) {
  foreach_particle_in(P) {
    p().x = 1.
  }
}
...
~~~

If all particles happen to be tracer_particles, this behaves identical to

~~~literatec
...
foreach_particle() {
  p().x = 1.
}
...
~~~

The Implementation makes use of `qcc`s excellent `foreach...`
functions. Meaning that these iterators can also do [simple
reductions](/Basilisk%20C#parallel-programming).
*/
#define p() pl[_l_particle][_j_particle] 

@def PARTICLE_VARIABLES
  double x = pl[_l_particle][_j_particle].x; NOT_UNUSED(x);
  double y = pl[_l_particle][_j_particle].y; NOT_UNUSED(y);
  double z = pl[_l_particle][_j_particle].z; NOT_UNUSED(z);
@

macro foreach_particle() {
  {
    int _l_particle = 0;
    while (pn[_l_particle] != terminate_int) {
      for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
	PARTICLE_VARIABLES;
	{...}
      }
      _l_particle++;
    }
  }
}

macro foreach_particle_in (int PARTICLES, Reduce reductions = None) {
  {
    int _l_particle = PARTICLES;
    for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
      PARTICLE_VARIABLES;
      {...}
    }
  }
}


macro foreach_P_in_list (int * PARTICLES_LIST) {
  {
    int _lll_particle = 0, _ll_particle = 0;
    while (pn[_lll_particle] != terminate_int) {
      if (_lll_particle == PARTICLES_LIST[_ll_particle]) {
	Particles P = _lll_particle; NOT_UNUSED(P);
	{...}
	_ll_particle++;
      }
      _lll_particle++;
    }
  }
}

#define remove_particles(p,CONDITION) do {		\
  int nr = 0;						\
  foreach_particle_in ((p)) {				\
    if ((CONDITION))					\
      nr++;						\
  }							\
  long unsigned int ind[nr];				\
  int i = 0;						\
  foreach_particle_in ((p)) {				\
    if ((CONDITION))					\
      ind[i++] = _j_particle;					\
  }							\
  remove_particles_index ((p), nr, ind);		\
  } while (0)


/**
## MPI particle exchange 

With `_MPI`, particles may find themselves outside the realm of their
thread's domain. They are then sent to the rank that owns their new
position. A particle outside a non-periodic domain is lost.

The data type below sends whole `particle`s, so that every member
(position, Runge-Kutta velocities and locations, tag and anything
added with `ADD_PART_MEM`) migrates with the particle. Counts are in
particles rather than bytes, which avoids `int` overflow for large
messages. */

/**
Some counters, for diagnostics and the [tests](test/run-tests.sh). */

#define PARTICLE_STATS 1
struct {
  long rounds;    // neighbour exchanges
  long messages;  // point-to-point messages sent
  long received;  // particles received from neighbours
  long halo;      // interpolations from the ghost layer (`particle_locate()`)
  long unlocated; // particles a rank could not interpolate
} particle_stats = {0};

#if _MPI
/**
qcc parses the code without the system's `mpi.h` and does not know
`MPI_File`. It sees the typedef below, while the C compiler, after
it, also sees the deferred `#define`, which makes every
`particle_file` an `MPI_File`. */

typedef void * particle_file;
@define particle_file MPI_File

static MPI_Datatype particle_mpi_type (void)
{
  static MPI_Datatype t = MPI_DATATYPE_NULL;
  if (t == MPI_DATATYPE_NULL) {
    MPI_Type_contiguous (sizeof(particle), MPI_BYTE, &t);
    MPI_Type_commit (&t);
  }
  return t;
}

#if MULTIGRID_MPI
/**
### Multigrid: neighbour-only migration

The [multigrid](/src/grid/multigrid-mpi.h) splits the domain into a
regular Cartesian array of `Dimensions.x` $\times$ `Dimensions.y`
($\times$ `Dimensions.z`) subdomains, one per rank, each a square
(cube) of side `L0/Dimensions.x` holding `1 << depth()` cells per
direction. `mpi_coords[]` is the position of this rank in the array.

The owner of any position therefore follows from its coordinates
alone. Velocities are interpolated from the rank's own data, using the
ghost cells that the fluid solver already exchanges, so interpolation
needs no extra communication either.

We use the Cartesian communicator of the multigrid itself. Its ranks
are those of `MPI_COMM_WORLD`, but Basilisk sends its own (ghost cell)
messages on `MPI_COMM_WORLD`, so the particle messages cannot be
mistaken for them. */

static MPI_Comm particle_cart_comm (void)
{
  for (Boundary ** b = boundaries; b && *b; b++)
    if ((*b)->level == mpi_boundary_level)
      return ((MpiBoundary *) *b)->cartcomm;
  assert (false); // the multigrid has no MPI boundary
  return MPI_COMM_NULL;
}

/**
`particle_owner()` sets `c[]` to the coordinates of the rank owning
particle `q`. In periodic directions, the position is first wrapped
into the domain. It returns `false` if `q` lies outside a non-periodic
domain. The cell index is computed as in `locate()`. */

static bool particle_owner (particle * q, int c[dimension])
{
  int n = 1 << depth();
  double * pos[3] = {&q->x, &q->y, &q->z};
  double o[3] = {X0, Y0, Z0};
  for (int d = 0; d < dimension; d++) {
    int nd = (&Dimensions.x)[d];
    double ext = L0*nd/Dimensions.x; // domain size in direction d
    if ((&Period.x)[d]) {
      *pos[d] -= ext*floor ((*pos[d] - o[d])/ext);
      if (*pos[d] >= o[d] + ext) // round-off
	*pos[d] = o[d];
    }
    double i = (*pos[d] - o[d])/L0*n*Dimensions.x;
    if (!(i > -1.) || i >= n*nd) // also catches NaN
      return false;
    c[d] = ((int) i)/n;
  }
  return true;
}

/**
Neighbours are indexed by their offset $o_d \in \{-1, 0, 1\}$ along
each direction, $k = \sum_d (o_d + 1)\,3^d$, so that $k = $ `PNB/2` is
this rank and the opposite neighbour of $k$ is `PNB - 1 - k`. In 2D
these are the eight neighbours, diagonals included. Ranks across a
non-periodic edge of the domain are `MPI_PROC_NULL`. */

#define PNB (dimension == 1 ? 3 : dimension == 2 ? 9 : 27)

static void particle_neighbours (MPI_Comm comm, int nbr[PNB])
{
  for (int k = 0; k < PNB; k++) {
    int c[dimension], m = k;
    bool valid = (k != PNB/2);
    for (int d = 0; d < dimension; d++, m /= 3) {
      int nd = (&Dimensions.x)[d];
      c[d] = mpi_coords[d] + m % 3 - 1;
      if (c[d] < 0 || c[d] >= nd) {
	if ((&Period.x)[d])
	  c[d] = (c[d] + nd) % nd;
	else
	  valid = false;
      }
    }
    nbr[k] = MPI_PROC_NULL;
    if (valid)
      MPI_Cart_rank (comm, c, &nbr[k]);
  }
}

/**
A particle moves a fraction of a cell per step, so on leaving its
subdomain it enters one of its neighbours. Each rank sends one
(possibly empty) message to every neighbour, tagged with the
neighbour's index, and receives one from each: communication scales
with the number of particles crossing subdomain boundaries, not with
the total, and does not grow with the number of ranks.

`update_mpi_lists()` exchanges the `nls` lists `ls` at once, with a
single message per neighbour: `nls` counts, then the particles of each
list in turn. The number of messages therefore does not grow with the
number of lists either. All ranks must pass the same lists.

A particle further away than a neighbour (time step far above the CFL
limit, or particles placed anywhere by the user) is forwarded one
neighbour towards its owner and reaches it after a few calls. Call
`mpi_redistribute_particles()` to place such particles at once. */

static void update_mpi_lists (Particles * ls, int nls)
{
  if (nls < 1)
    return;
  MPI_Comm comm = particle_cart_comm();
  int nbr[PNB];
  particle_neighbours (comm, nbr);

  // Keep local particles (in place), pack the leaving ones
  size_t hs = nls*sizeof(int), nbytes[PNB], nalloc[PNB] = {0};
  char * buf[PNB] = {NULL};
  int lost = 0, far = 0;
  for (int m = 0; m < nls; m++) {
    Particles l = ls[m];
    long unsigned int kept = 0;
    particle * pp = pl[l];
    for (long unsigned int j = 0; j < pn[l]; j++) {
      particle q = pp[j];
      int c[dimension];
      if (!particle_owner (&q, c)) {
	lost++;
	continue;
      }
      int k = 0, s = 1;
      bool hop = false;
      for (int d = 0; d < dimension; d++, s *= 3) {
	int nd = (&Dimensions.x)[d], o = c[d] - mpi_coords[d];
	if ((&Period.x)[d]) { // shortest way round
	  if (o > nd/2)
	    o -= nd;
	  else if (o < -(nd - 1)/2)
	    o += nd;
	}
	if (o > 1)
	  o = 1, hop = true;
	else if (o < -1)
	  o = -1, hop = true;
	k += (o + 1)*s;
      }
      if (k == PNB/2)
	pp[kept++] = q;
      else {
	assert (nbr[k] != MPI_PROC_NULL);
	far += hop;
	if (!buf[k]) {
	  nalloc[k] = hs + 16*sizeof(particle);
	  buf[k] = calloc (nalloc[k], 1);
	  nbytes[k] = hs;
	}
	if (nbytes[k] + sizeof(particle) > nalloc[k]) {
	  nalloc[k] *= 2;
	  buf[k] = realloc (buf[k], nalloc[k]);
	}
	memcpy (buf[k] + nbytes[k], &q, sizeof(particle));
	nbytes[k] += sizeof(particle);
	((int *) buf[k])[m]++;
      }
    }
    pn[l] = kept;
  }

  // Send to each neighbour, then receive from each
  MPI_Request req[PNB];
  int nr = 0, * zero = calloc (nls, sizeof(int)), * in = calloc (nls, sizeof(int));
  for (int k = 0; k < PNB; k++)
    if (nbr[k] != MPI_PROC_NULL)
      MPI_Isend (buf[k] ? buf[k] : (char *) zero, buf[k] ? nbytes[k] : hs,
		 MPI_BYTE, nbr[k], k, comm, &req[nr++]);
  char * rbuf[PNB] = {NULL};
  for (int k = 0; k < PNB; k++)
    if (nbr[PNB - 1 - k] != MPI_PROC_NULL) {
      MPI_Status st;
      int size;
      MPI_Probe (nbr[PNB - 1 - k], k, comm, &st);
      MPI_Get_count (&st, MPI_BYTE, &size);
      rbuf[k] = malloc (size);
      MPI_Recv (rbuf[k], size, MPI_BYTE, nbr[PNB - 1 - k], k, comm,
		MPI_STATUS_IGNORE);
      for (int m = 0; m < nls; m++)
	in[m] += ((int *) rbuf[k])[m];
    }

  // Append the arrivals to each list, in neighbour order
  for (int m = 0; m < nls; m++) {
    Particles l = ls[m];
    long unsigned int n0 = pn[l];
    change_plist_size (l, in[m]);
    particle * dst = pl[l] + n0;
    for (int k = 0; k < PNB; k++)
      if (rbuf[k]) {
	int * h = (int *) rbuf[k];
	size_t off = hs;
	for (int mm = 0; mm < m; mm++)
	  off += h[mm]*sizeof(particle);
	memcpy (dst, rbuf[k] + off, h[m]*sizeof(particle));
	dst += h[m];
	particle_stats.received += h[m];
      }
  }
  MPI_Waitall (nr, req, MPI_STATUSES_IGNORE);
  particle_stats.rounds++;
  particle_stats.messages += nr;
  for (int k = 0; k < PNB; k++)
    free (buf[k]), free (rbuf[k]);
  free (zero), free (in);

  if (lost)
    fprintf (stderr, "#pid: %d lost %d particles outside the domain\n",
	     pid(), lost);
  if (far)
    fprintf (stderr, "#pid: %d forwarded %d particles beyond its neighbours"
	     " (consider mpi_redistribute_particles())\n", pid(), far);
}

/**
### Global redistribution

Newly created lists (from `new_particles()`, which every rank calls,
with zero particles on most ranks for e.g. `init_tp_square()` or
`prestore()`) may hold particles anywhere. Their first update sends
every particle directly to its owner with a single all-to-all. This
function is collective and may also be called by the user, e.g. after
injecting particles on one rank. */

void mpi_redistribute_particles (Particles l)
{
  MPI_Comm comm = particle_cart_comm();
  MPI_Datatype ptype = particle_mpi_type();
  int np = npe(), lost = 0;
  int * scount = calloc (np, sizeof(int)), * rcount = malloc (np*sizeof(int));
  int * sdispl = malloc (np*sizeof(int)), * rdispl = malloc (np*sizeof(int));
  int * dest = malloc ((pn[l] + 1)*sizeof(int));
  for (long unsigned int j = 0; j < pn[l]; j++) {
    int c[dimension];
    if (particle_owner (&pl[l][j], c)) {
      MPI_Cart_rank (comm, c, &dest[j]);
      scount[dest[j]]++;
    }
    else
      dest[j] = -1, lost++;
  }
  sdispl[0] = 0;
  for (int r = 1; r < np; r++)
    sdispl[r] = sdispl[r - 1] + scount[r - 1];
  particle * sbuf = malloc ((pn[l] - lost + 1)*sizeof(particle));
  int * cur = malloc (np*sizeof(int));
  memcpy (cur, sdispl, np*sizeof(int));
  for (long unsigned int j = 0; j < pn[l]; j++)
    if (dest[j] >= 0)
      sbuf[cur[dest[j]]++] = pl[l][j];

  MPI_Alltoall (scount, 1, MPI_INT, rcount, 1, MPI_INT, comm);
  rdispl[0] = 0;
  for (int r = 1; r < np; r++)
    rdispl[r] = rdispl[r - 1] + rcount[r - 1];
  int in = rdispl[np - 1] + rcount[np - 1];
  pn[l] = 0;
  change_plist_size (l, in);
  particle dummy;
  MPI_Alltoallv (sbuf, scount, sdispl, ptype,
		 in ? pl[l] : &dummy, rcount, rdispl, ptype, comm);
  p_global[l] = false;

  free (scount); free (rcount); free (sdispl); free (rdispl);
  free (dest); free (sbuf); free (cur);
  if (lost)
    fprintf (stderr, "#pid: %d lost %d particles outside the domain\n",
	     pid(), lost);
}

#else // !MULTIGRID_MPI
/**
### Tree grids: global exchange

The decomposition of a tree is not a regular array, so the leaving
particles are gathered by all ranks, each keeping those it can
locate. */

static void update_mpi_global (Particles p) {
  int l = 0;
  while (pn[l] != terminate_int) {
    if (l == p) {
      int outt, in = 0, m = 0, out = 0;
      foreach_particle_in(p) 
	if (locate (p().x, p().y, p().z).level < 0) {
	  out++;
	}
      //get indices and outgoing data
      int * ind = malloc (sizeof(int)*out);
      particle * senddata = malloc (sizeof(particle)*out);
      foreach_particle_in(p) { 
	if (locate (p().x, p().y, p().z).level < 0) {
	  ind[m] = _j_particle;
	  senddata[m++] = p();
	}
      }
      //remove the senddata from arrays (shrink)
      m = 0;
      int j = 0;
      while (j < pn[l] - out) {
	while (m < out ? j + m == ind[m] : 0)
	  m++;
	while (m < out ? j < pn[l] - out && j + m != ind[m] : j < pn[l] - out) {
	  pl[l][j]   = pl[l][j + m];
	  j++;
	}
      }
      // Gather lost particles among threads:
      // First, count all of them
      int outa[npe()], outat[npe()];
      outat[0] = 0;
      
      MPI_Allgather (&out, 1, MPI_INT, &outa[0], 1, MPI_INT, MPI_COMM_WORLD);
      //Compute displacements
      for (int j = 1; j < npe(); j++) 
	outat[j] = outa[j - 1] + outat[j - 1];
      outt = outat[npe() - 1] + outa[npe() - 1]; 
      // Allocate recieve buffer and gather
      particle * recdata = malloc (sizeof(particle)*outt);
      for (int j = 0; j < npe(); j++) {
	outat[j] *= sizeof(particle);
	outa[j]  *= sizeof(particle);
      }
    //send and recieve data
      MPI_Allgatherv (&senddata[0], outa[pid()], MPI_BYTE,
		      &recdata[0], outa, outat, MPI_BYTE,
		      MPI_COMM_WORLD); 
    //count new particles
    for (int j = 0; j < outt ; j++) 
      if (locate (recdata[j].x, recdata[j].y, recdata[j].z).level > 0)
	in++;
    long unsigned int po = pn[l] - out;
    //Manage the memory if required...
    change_plist_size (l, in - out);
    //Collect new particles from `recdata`
    if (in > 0) {
      int indi[in];
      m = 0;
      for (int j = 0; j < outt; j++) 
	if (locate (recdata[j].x, recdata[j].y, recdata[j].z).level > 0)
	  indi[m++] = j;
      m = 0;
      for (j = po; j < pn[l]; j++) {
	pl[l][j] = recdata[indi[m]];
	m++;
      }
    }
    // clean the mess
    free (ind); free (senddata); free (recdata);
    }
    l++;
  }
  
}

void mpi_redistribute_particles (Particles p) {
  update_mpi_global (p);
  p_global[p] = false;
}

static void update_mpi_lists (Particles * ls, int nls) {
  for (int m = 0; m < nls; m++)
    update_mpi_global (ls[m]);
}
#endif // !MULTIGRID_MPI

/**
### Deferred exchanges

`particle_boundary()` (in [particle_classic.h]() and
[particle_multilayer.h]()) ends with `update_mpi()`. While `p_defer` is
set, `update_mpi()` only records the list, so that
`particles_boundary()` below can exchange all the lists at once. A new
list (`p_global`) is redistributed straight away, which is collective
and consistent since all ranks create the same lists. */

static int p_defer = 0;
static Particles * p_pending = NULL;
static int p_npending = 0;

void update_mpi (Particles p) {
  if (p_global[p])
    mpi_redistribute_particles (p);
  else if (p_defer) {
    for (int m = 0; m < p_npending; m++)
      if (p_pending[m] == p)
	return;
    p_pending = realloc (p_pending, (p_npending + 1)*sizeof(Particles));
    p_pending[p_npending++] = p;
  }
  else
    update_mpi_lists (&p, 1);
}
#endif // _MPI

/**
## Locating a particle for interpolation

`particle_locate()` is `locate()`, except that on the multigrid with
MPI it also accepts the first layer of ghost cells across an edge
shared with another rank, or across a periodic edge. The ghost values
there are exact copies of the neighbour's (the exchange sends both
ghost layers, corners included), and bilinear interpolation from such
a cell reads at most the second layer. A particle that has just left
its subdomain is then interpolated by its own rank with the same
result as by its owner, which lets `tracer_particles_step1` skip an
exchange. On a closed wall there is no copy, so the particle is not
found there. `q` is changed to the periodic image of the position that
lies in the ghost layer, which is the position to interpolate at.

With the multilayer solver, `particle_locate()` also selects the layer
holding the particle's height `q->z`: it sets `point.l`, so that fields
read at the point, as in `interpolate_linear()`, are those of that
layer. Layer `k` spans from $z_b + h_0 + \dots + h_{k-1}$ to $z_b + h_0 +
\dots + h_k$ in the particle's cell. A height below the bottom or above
the surface gives the bottom or the top layer. Within a layer the
velocity is constant in the vertical, as in the solver; horizontally it
is interpolated from the same layer of the neighbouring cells. This reads
`zb` and `h`, also in the ghost cells: call `particle_locate_boundary()`
before locating particles.

### Surface tracers (`PARTICLE_SURFACE`)

With the multilayer solver, `PARTICLE_SURFACE` defined to 1 (before
including [tracer-particles.h](), or `-DPARTICLE_SURFACE=1`) makes every
particle a passive tracer of the free surface. The surface moves with
the water, so a surface tracer stays on it, and its velocity is the
fluid velocity there: that of the top layer, since the solver's velocity
is constant within a layer. So, with the switch:

* `particle_locate()` sets `point.l` to the top layer, whatever `z`,
  and `zb` and `h` are not read (no layer search, no ghost-cell update
  for them);
* `particle_boundary()` in [particle_multilayer.h]() does not reflect
  `z` at the surface and the bottom;
* `z` is neither used nor updated: outputs write the surface height
  `eta` at the particle instead (`particle_surface_height()`);
* `init_tp_cells()` puts one particle per cell, at the surface.

Particles cannot be in the bulk of the flow in the same program: the
switch applies to all of them. Without it, none of this is compiled. */

#if LAYERS && !PARTICLE_SURFACE
static Point particle_layer (Point point, double zp) // z is reserved here
{
  double top = zb[];
  for (point.l = 0; point.l < nl - 1; point.l++)
    if (zp < (top += h[]))
      break;
  return point;
}
#endif

/**
`particle_locate_boundary()` updates the ghost cells that
`particle_locate()` reads (with the multilayer solver, those of `zb` and
`h`). It is collective, and costs nothing when they are up to date. */

void particle_locate_boundary (void)
{
#if LAYERS && !PARTICLE_SURFACE
  boundary ({zb, h});
#endif
}

static Point particle_locate_cell (coord * q);

Point particle_locate (coord * q)
{
  Point point = particle_locate_cell (q);
#if LAYERS && PARTICLE_SURFACE
  if (point.level >= 0)
    point.l = nl - 1; // surface tracers: always the top layer
#elif LAYERS
  if (point.level >= 0)
    point = particle_layer (point, q->z);
#endif
  return point;
}

#if LAYERS && PARTICLE_SURFACE
/**
`particle_surface_height (xp, yp)` is the height of the free surface
`eta` at a position, interpolated as the velocity is, or `NAN` if this
rank cannot locate it. It reads the ghost cells of `eta`: call
`boundary ({eta})` first, collectively. */

double particle_surface_height (double xp, double yp)
{
  coord q = {xp, yp, 0.};
  Point point = particle_locate (&q);
  if (point.level < 0)
    return NAN;
  double xq = q.x, yq = q.y, zq = q.z;
  return interpolate_linear (point, eta, xq, yq, zq);
}
#endif

static Point particle_locate_cell (coord * q)
{
  Point point = locate (q->x, q->y, q->z);
#if MULTIGRID_MPI
  if (point.level >= 0)
    return point;
  int n = 1 << depth(), ind[3] = {0, 0, 0};
  double pos[3] = {q->x, q->y, q->z}, o[3] = {X0, Y0, Z0};
  for (int d = 0; d < dimension; d++) {
    int nd = (&Dimensions.x)[d], c = mpi_coords[d];
    bool per = (&Period.x)[d], found = false;
    double ext = L0*nd/Dimensions.x;
    for (int s = 0; s < (per ? 3 : 1) && !found; s++) {
      double xs = pos[d] + (s == 1 ? - ext : s == 2 ? ext : 0.);
      double a = (xs - o[d])/L0*n*Dimensions.x + GHOSTS - c*n; // as locate()
      if (!(a >= GHOSTS - 1) || a >= n + GHOSTS + 1)
	continue;
      int i = a;
      if (!per && ((i < GHOSTS && c == 0) || (i >= n + GHOSTS && c == nd - 1)))
	continue; // a closed wall
      ind[d] = i, pos[d] = xs, found = true;
    }
    if (!found)
      return point; // level -1
  }
  point.level = depth();
  point.i = ind[0];
#if dimension > 1
  point.j = ind[1];
#endif
#if dimension > 2
  point.k = ind[2];
#endif
  q->x = pos[0], q->y = pos[1], q->z = pos[2];
  particle_stats.halo++;
#endif // MULTIGRID_MPI
  return point;
}

/**
`particle_is_mine()` tells whether this rank owns position `p`, so
that when every rank tests the same positions, each position is kept
by exactly one rank (none if it lies outside the domain). On the
multigrid, every rank computes the owner with the same arithmetic. On a
tree (quadtree, octree), the owner is the rank where `locate()` finds
the position. Basilisk defines `MULTIGRID` for trees too (they are built
on the multigrid code), hence the test on `TREE`. */

bool particle_is_mine (coord p)
{
#if MULTIGRID && !TREE // the regular multigrid only, not trees
  double pos[3] = {p.x, p.y, p.z}, o[3] = {X0, Y0, Z0};
  for (int d = 0; d < dimension; d++)
    if (!(pos[d] >= o[d]) ||
	pos[d] >= o[d] + L0*(&Dimensions.x)[d]/Dimensions.x)
      return false;
#if MULTIGRID_MPI
  particle q = {0};
  q.x = p.x, q.y = p.y, q.z = p.z;
  int c[dimension];
  if (!particle_owner (&q, c))
    return false;
  for (int d = 0; d < dimension; d++)
    if (c[d] != mpi_coords[d])
      return false;
#endif
  return true;
#else // trees and other grids
  return locate (p.x, p.y, p.z).level >= 0;
#endif
}




// Some shared placement functions

void place_in_square (Particles p, int n = 10,
		      double xm = 0, double ym = 0,
		      double l = L0, double zp = 0)
{
  long unsigned int j = 0;
  particles pp = pl[p];
  double dx = l/(double)n;
  double x0 = xm - l/2. + dx/2;
  double y0 = ym - l/2. + dx/2;
  for (double xp = x0; xp < xm + l/2.; xp += dx) {
    for (double yp = y0; yp < ym + l/2.; yp += dx) {
      pp[j].x = xp;
      pp[j].y = yp;
      pp[j++].z = zp;
    }
  }
}

void place_in_circle (Particles p, int n = 10,
		      double xm = 0, double ym = 0,
		      double l = L0, double zp = 0) {
  particles pp = pl[p];
  for (int j = 0; j < n; j++) {
    double angle = noise()*pi;
    double R = sqrt(fabs(noise()))*l;
    pp[j].x = R*sin(angle) + xm;
    pp[j].y = R*cos(angle) + ym;
    pp[j].z = zp;
  }
}

// A shared particle statistics struct
typedef struct {
  coord avg;
  coord min;
  coord max;
  coord stddev;
} pstats;

/**
### Dump and restore particles
The `pdump()` and `prestore()` functions are implemented below.
 */

/**
After the particles, `pdump()` writes a trailer: 8 bytes
`PDUMP_TRAILER` and the tag counter `p_ntags` of each list. Readers that
do not know it (such as the original `prestore()`) ignore it. */

#define PDUMP_TRAILER "pntags1"

bool pdump (char * fname = "pdump",  //File name
	    Particles * td = NULL,   //Particle lists  Default: all
	    FILE * fp = NULL,        //File pointer    Default: Not used
	    bool Dmem = false,       //Member data     Default: false, only x,y,z
	    bool print = false) {    //Show dump data  Default: false
  // Get particle lists data
  int j = 0, n = 0;
  bool all = !td;
  if (td) {
    foreach_P_in_list (td) 
      j++;
  } else  {// all
    while (pn[j++] != terminate_int);
    td = malloc (sizeof(int)*j);
    j = 0;
    while (pn[j] != terminate_int) { 
      td[j] = j;
      j++;
    }
  }

  // Nr of particles
  long unsigned int np[j];
  foreach_P_in_list (td) {
    np[n] = pn[P];
    n++;
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, np, j, MPI_UNSIGNED_LONG, MPI_SUM,
		 MPI_COMM_WORLD);
#endif
  int size = Dmem ? sizeof(particle) : 3*sizeof(double);
  int Headersize = sizeof(int) + j*sizeof (long unsigned int) + sizeof(bool); 
  if (print && pid() == 0) 
    fputs ("# P\tamount\n", stderr);

#if _MPI
  /**
  With MPI and no file pointer, the ranks write their particles in
  place through MPI-IO, one after the other in rank order (as
  before). The file is created and emptied once, collectively: with
  `fopen (file, "wb")` on every rank, a late rank could empty the
  file after others had written to it. */

  if (!fp) {
    particle_file fh;
    if (MPI_File_open (MPI_COMM_WORLD, fname,
		       MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
		       &fh) != MPI_SUCCESS) {
      fprintf (stderr, "pdump(): cannot open '%s'\n", fname);
      exit (1);
    }
    MPI_File_set_size (fh, 0);
    if (pid() == 0) {
      char head[Headersize];
      memcpy (head, &j, sizeof(int));
      memcpy (head + sizeof(int), np, j*sizeof(long unsigned int));
      memcpy (head + sizeof(int) + j*sizeof(long unsigned int), &Dmem,
	      sizeof(bool));
      MPI_File_write_at (fh, 0, head, Headersize, MPI_BYTE,
			 MPI_STATUS_IGNORE);
    }
    MPI_Datatype rec;
    MPI_Type_contiguous (size, MPI_BYTE, &rec);
    MPI_Type_commit (&rec);
    MPI_Offset offset = Headersize;
    for (int m = 0; m < j; m++) {
      Particles P = td[m];
      long unsigned int nl = pn[P], before = 0;
      MPI_Exscan (&nl, &before, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (pid() == 0)
	before = 0;
      char * buf = (char *) pl[P], dummy;
      if (!Dmem) {
	double * c = malloc ((3*nl + 1)*sizeof(double));
	foreach_particle_in (P)
	  c[3*_j_particle] = x, c[3*_j_particle + 1] = y,
	    c[3*_j_particle + 2] = z;
	buf = (char *) c;
      }
      MPI_File_write_at_all (fh, offset + (MPI_Offset) before*size,
			     nl ? buf : &dummy, nl, rec, MPI_STATUS_IGNORE);
      if (!Dmem)
	free (buf);
      if (print && pid() == 0) 
	fprintf (stderr, "# %d\t%ld\n", m, np[m]);
      offset += (MPI_Offset) np[m]*size;
    }
    if (pid() == 0) {
      long unsigned int nt[j];
      for (int m = 0; m < j; m++)
	nt[m] = p_ntags[td[m]];
      MPI_File_write_at (fh, offset, PDUMP_TRAILER, 8, MPI_BYTE,
			 MPI_STATUS_IGNORE);
      MPI_File_write_at (fh, offset + 8, nt, j*sizeof(long unsigned int),
			 MPI_BYTE, MPI_STATUS_IGNORE);
    }
    MPI_Type_free (&rec);
    MPI_File_close (&fh);
    if (all)
      free (td);
    return true;
  }
#endif // _MPI

  // The file:
  char * file = fp ? NULL : fname;
  if (file) {
    if ((fp = fopen (file, "wb")) == NULL) {
      perror (file);
      exit (1);
    }
  }
  assert (fp);
  
  // Print header
  if (pid() == 0) {
    fwrite (&j, sizeof(int), 1, fp);
    fwrite (np, sizeof(long unsigned int), j, fp);
    fwrite (&Dmem, sizeof(bool), 1, fp);
  }
  fseek (fp, Headersize, SEEK_SET); //offset

  // Print particle data
  for (int m = 0; m < j; m++) { //Nesting within foreach_P... did not work
    Particles P = td[m];
#if _MPI
    long unsigned int outa[npe()], outat[npe()  + 1];
    outat[0] = 0;
    MPI_Allgather (&pn[td[m]], 1, MPI_UNSIGNED_LONG,
		   &outa[0], 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    //Compute and set displacements
    for (int j = 1; j <= npe(); j++) 
      outat[j] = outa[j - 1] + outat[j - 1];
    fseek (fp, outat[pid()]*size, SEEK_CUR);
#endif
    foreach_particle_in (P) {
      if (Dmem) {
	fwrite (&p(), sizeof(particle), 1, fp);
      } else {
	double c[3] = {x, y, z};
	fwrite (c, sizeof(double), 3, fp);
      }
    }
    if (print && pid() == 0) 
      fprintf (stderr, "# %d\t%ld\n", m, np[m]);
#if _MPI
    fseek (fp, (np[m] - outat[pid() + 1])*size, SEEK_CUR);
#endif
  }
  if (pid() == 0) {
    fwrite (PDUMP_TRAILER, 1, 8, fp);
    for (int m = 0; m < j; m++)
      fwrite (&p_ntags[td[m]], sizeof(long unsigned int), 1, fp);
  }
  fclose (fp);
  if (all)
    free (td);
  return true; 
}

/**
These include statements separate out all functions which behave differently in 
the multilayer solver compared to the other solvers in basilisk C. The multilayer patch was developed by [larswd](basilisk.fr/sandbox/larswd/README). The [multilayer code](particle_multilayer.h) and the [code for the other solvers](particle_classic.h) are separated for ease of comparison.
*/

#if LAYERS
#include "particle_multilayer.h"
#else
#include "particle_classic.h"
#endif

/**
## Boundary conditions for several lists at once

`particles_boundary()` applies `particle_boundary()` to each list of
`plist` and then exchanges them all at once (one message per
neighbour, see `update_mpi_lists()`).

With `exchange = false` on the multigrid, the exchange is skipped and
the particles stay on their rank. This is only correct if they cannot
have moved further than the first ghost layer since the last exchange:
they are then interpolated from it (see `particle_locate()`). Other
grids always exchange. */

void particles_boundary (Particles * plist, bool exchange = true)
{
  if (!plist || !pn)
    return;
#if _MPI
  p_defer++;
#endif
  foreach_P_in_list (plist)
    particle_boundary (P);
#if _MPI
  p_defer--;
#if MULTIGRID_MPI
  if (exchange)
#endif
    update_mpi_lists (p_pending, p_npending);
  p_npending = 0;
#endif
}
/**
#### Particle Restoration 

The restore function forgets all existing particles and does not
assign particles to any list, nor does it know about the names of the
particle lists references. See also [this test](test_prestore.c). */
int prestore (char * fname = "pdump",  //File name
	      FILE * fp = NULL,        //File pointer    Default: Not used
	      bool print = false)      //Show dump data  Default: false
{
  // Open file
  if (fp) fname = NULL;
  else if (fname && (fp = fopen (fname, "r")) == NULL)
    return 0;
  assert (fp);

  //read nr. of lists and nr. of particles
  long base = ftell (fp);
  int j;
  bool Dmem;
  fread (&j, sizeof(int), 1, fp);
  long unsigned int np[j];
  fread (np, sizeof(long unsigned int), j, fp);
  fread (&Dmem, sizeof(bool), 1, fp);
  
  // Print some reference data
  if (print) {
    if (pid() == 0) {
      if (Dmem)
	fputs ("# Restoring members...\n", stderr);
      fputs ("# P\tamount\n", stderr);
      for (int m = 0; m < j; m++)
	fprintf (stderr, "# %d\t%ld\n", m, np[m]);
    }
  }
  
  // Remove existing particles
  if (pl != NULL)
    free_p();

  /**
  Each rank reads a contiguous slice of every list (with MPI; all of it
  otherwise), and the first boundary condition sends the particles to
  their owners (see `mpi_redistribute_particles()`). */

  int size = Dmem ? sizeof(particle) : 3*sizeof(double);
  long offset = base + sizeof(int) + j*sizeof(long unsigned int) + sizeof(bool);
  Particles P;
  for (int m = 0; m < j; m++) {
    long unsigned int first = 0, count = np[m];
#if _MPI
    first = np[m]*pid()/npe();
    count = np[m]*(pid() + 1)/npe() - first;
#endif
    P = new_particles (count);
    fseek (fp, offset + first*size, SEEK_SET);
    if (Dmem) 
      fread (pl[P], sizeof(particle), count, fp);
    else {
      foreach_particle_in (P) {
	double c[3];
	fread (&c, sizeof(double), 3, fp);
	p().x = c[0]; p().y = c[1]; p().z = c[2];
      }
    }
    offset += np[m]*size;
    // Apply boundary conditions.
    particle_boundary (P);
  }

  // Tag counters, from the trailer when the dump has one
  char trailer[8] = {0};
  fseek (fp, offset, SEEK_SET);
  bool tagged = fread (trailer, 1, 8, fp) == 8 &&
    !strncmp (trailer, PDUMP_TRAILER, 8);
  for (int m = 0; m < j; m++)
    if (!tagged || fread (&p_ntags[m], sizeof(long unsigned int), 1, fp) != 1)
      p_ntags[m] = P_NTAGS_UNKNOWN;
  if (!Dmem && pid() == 0)
    fputs ("prestore(): positions only: the other members (tags, velocities)"
	   " are zero; dump with Dmem = true to keep them\n", stderr);
  if (fname)
    fclose (fp);
  return j;
}


/**
## todo

* ~~~Dynamically add and delete single `particle`s and `Particles` lists~~~
* Sort particles along grid-iterator curve ([radix.h]())
* More flexible attributes/members for `particle`s. 
* ~~~[Tie particles to Basilisk fields](tie_coord_to_cell.c)~~~
* ~~~dump and restore particles~~~


## Tests

* [Dump and restore particles with MPI](test_prestore.c)
* [Length of a loop](tlengt.c)
* [Adding and removing particles](inject_particles.c)

## Usage

* [flow tracer particles](tracer-particles.h)
* [Inertial particles](inertial-particles.h)
* [Brownian motion](brownian.c)

## See also

* [Bview draw function for particles](scatter2.h)
* [Reference particles using scalar data](particle_reference.h)
 */
