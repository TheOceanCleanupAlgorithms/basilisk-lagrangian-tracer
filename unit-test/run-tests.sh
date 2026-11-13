#!/bin/bash
# Unit tests of the particle headers in the folder above.
#
#   ./run-tests.sh            build and run everything
#   ./run-tests.sh build      build only
#
# Builds, raw results and the pass/fail list ($WORK/results/checks.txt)
# go to $WORK. Ranks above the core count are run with --oversubscribe.
set -e
cd "$(dirname "$0")"
BASILISK=${BASILISK:-$HOME/basilisk/src}
export BASILISK PATH=$PATH:$BASILISK
PYTHON=${PYTHON:-python3}     # with numpy
MPIRUN="mpirun --oversubscribe"
WORK=${WORK:-$HOME/Library/Caches/particle-tests}
mkdir -p $WORK

if [ "$(ps aux | grep -c '[n]orth_pacific-mpi')" != 0 ]; then
  echo "a north_pacific run is going: not competing with it" >&2
  exit 1
fi

build() { # variant headers...
  local v=$1; shift
  rm -rf $WORK/build/$v; mkdir -p $WORK/build/$v
  cp "$@" tmig.c tml.c $WORK/build/$v/
  (cd $WORK/build/$v
   for b in "" "-DBOX=1" "-DWRITER=1"; do
     case "$b" in "") n=tmig;; -DBOX=1) n=tmig-box;; *) n=tmigw;; esac
     qcc -O2 -disable-dimensions -grid=multigrid $b tmig.c -o $n-serial -lm 2> $n.log
     CC99='mpicc -std=c99 -D_DARWIN_C_SOURCE' qcc -O2 -disable-dimensions -grid=multigrid $b -D_MPI=1 \
	 tmig.c -o $n-mpi -lm 2>> $n.log
   done
   for b in "" "-DCELLS=1" "-DCELLS=1 -DSHEAR=1"; do
     case "$b" in "") n=tml;; -DCELLS=1) n=tmlc;; *) n=tmls;; esac
     qcc -O2 -disable-dimensions $b tml.c -o $n-serial -lm 2> $n.log
     CC99='mpicc -std=c99 -D_DARWIN_C_SOURCE' qcc -O2 -disable-dimensions $b -D_MPI=1 \
	 tml.c -o $n-mpi -lm 2>> $n.log
   done
   for n in tmig tmig-box tmigw tml tmlc tmls; do
     for m in serial mpi; do
       [ -x $n-$m ] || { echo "build of $v/$n-$m failed, see $WORK/build/$v/$n.log" >&2; exit 1; }
     done
   done)
  echo "built $v"
}

build after ../particle.h ../particle_classic.h ../particle_multilayer.h ../tracer-particles.h ../write_class_particles.h
# the switches of the current code only: positions in degrees
# (PARTICLE_GEOGRAPHY, tsph) and surface tracers (PARTICLE_SURFACE, tmlsurf)
cp tsph.c $WORK/build/after/
(cd $WORK/build/after
 for b in "" "-DWRITER=1 -DWCP_LON_WRAPPED=1" \
	  "-DCELLS=1 -DSHEAR=1 -DPARTICLE_SURFACE=1" \
	  "-DCELLS=1 -DWRITER=1 -DPARTICLE_SURFACE=1" "-grid=quadtree" \
	  "-grid=quadtree -DWRITER=1"; do
   case "$b" in
     "") n=tsph; s=tsph.c;;
     *WCP_LON*) n=tsphw; s=tsph.c;;
     *SHEAR*) n=tmlsurf; s=tml.c;;
     -grid=quadtree) n=tmigq; s=tmig.c;;
     "-grid=quadtree -DWRITER=1") n=tmigqw; s=tmig.c;;
     *) n=tmlsurfw; s=tml.c;;
   esac
   qcc -O2 -disable-dimensions $b $s -o $n-serial -lm 2> $n.log
   CC99='mpicc -std=c99 -D_DARWIN_C_SOURCE' qcc -O2 -disable-dimensions $b -D_MPI=1 \
       $s -o $n-mpi -lm 2>> $n.log
   for m in serial mpi; do
     [ -x $n-$m ] || { echo "build of after/$n-$m failed, see $WORK/build/after/$n.log" >&2; exit 1; }
   done
 done)
# Order in time (tconv, run.h) and the time of u with a real solver
# (tcentered, navier-stokes/centered.h)
cp tconv.c tcentered.c $WORK/build/after/
(cd $WORK/build/after
 qcc -O2 -disable-dimensions -grid=multigrid tconv.c -o tconv -lm 2> tconv.log
 qcc -O2 -disable-dimensions -grid=multigrid tcentered.c -o tcentered -lm 2> tcentered.log
 for n in tconv tcentered; do
   [ -x $n ] || { echo "build of after/$n failed, see $WORK/build/after/$n.log" >&2; exit 1; }
 done)
# Particles in a Stokes wave with layered/nh.h (vertical motion with w):
# tnh (solved flow) and tnhp (prescribed flow)
cp tnh.c stokes_wave/stokes.h $WORK/build/after/
(cd $WORK/build/after
 for p in "" "-DPRESCRIBED=1"; do
   n=tnh${p:+p}
   qcc -O2 -disable-dimensions $p tnh.c -o $n-serial -lm 2> $n.log
   CC99='mpicc -std=c99 -D_DARWIN_C_SOURCE' qcc -O2 -disable-dimensions $p -D_MPI=1 \
       tnh.c -o $n-mpi -lm 2>> $n.log
   for m in serial mpi; do
     [ -x $n-$m ] || { echo "build of after/$n-$m failed, see $WORK/build/after/$n.log" >&2; exit 1; }
   done
 done)
[ "$1" = build ] && exit 0

B=$WORK/build
R=$WORK/results
rm -rf $R; mkdir -p $R

# run <dir> <ranks> <program> args...: merged, sorted outputs in <dir>
# <program> is <variant>/<name>, e.g. after/tmig; -serial or -mpi is added
run() {
  local d=$1 np=$2 bin=$B/$3; shift 3
  rm -rf $d; mkdir -p $d
  (cd $d
   set +e
   if [ $np = 1 ]; then $bin-serial "$@" > stdout.txt 2> log.err
   else $MPIRUN -np $np $bin-mpi "$@" > stdout.txt 2> log.err; fi
   status=$?
   [ $status = 0 ] || echo $status > crashed
   grep -v '^#' stdout.txt | grep -v '^output' | grep -v '^ntags' > stats.txt
   grep '^output' stdout.txt > output.txt
   grep '^ntags' stdout.txt > ntags.txt
   for f in out restored; do
     cat $f-* 2>/dev/null | sort -n > $f.txt; rm -f $f-*
   done
   cat tags-* 2>/dev/null | sort -k1,1n -k2,2n > tags.txt; rm -f tags-*
   cat traj-* 2>/dev/null | sort -k1,1n -k3,3n > traj.txt; rm -f traj-*
   [ -s restored.txt ] || rm -f restored.txt
   true)
  echo "$d: $(cat $d/stats.txt)$([ -f $d/crashed ] && echo CRASHED) | warnings: $(grep -c 'lost\|forwarded\|could not' $d/log.err)"
}
check() { # <description> <file> <file>
  if cmp -s "$2" "$3"; then echo "PASS $1"; else echo "FAIL $1"; fi | tee -a $R/checks.txt
}

# Accuracy: trajectories every 10 steps, periodic square and closed box
for v in after; do
  for np in 1 4 16; do
    run $R/accuracy/periodic/$v-np$np $np $v/tmig snap=10
    run $R/multilayer/$v-np$np $np $v/tml
  done
  for np in 1 2 8; do
    run $R/accuracy/box/$v-np$np $np $v/tmig-box snap=10
  done
done
# Four lists (M2): 4 x 400 particles
for v in after; do
  for np in 1 4 16; do
    run $R/accuracy/lists/$v-np$np $np $v/tmig nl=4 np=20 snap=10
  done
done
for np in 4 16; do
  check "periodic, $np ranks = serial" $R/accuracy/periodic/after-np{1,$np}/traj.txt
  check "four lists, $np ranks = serial" $R/accuracy/lists/after-np{1,$np}/traj.txt
  check "multilayer, $np ranks = serial" $R/multilayer/after-np{1,$np}/out.txt
done
# Multilayer, sheared layers: each particle moves with its own layer
echo ok > $R/shear-ok.txt
for np in 1 4 16; do
  run $R/multilayer/shear-np$np $np after/tmls
  awk '/^shear/ {print $2}' $R/multilayer/shear-np$np/stdout.txt > $R/shear-np$np.txt
  check "multilayer, sheared layers on $np ranks: each layer's velocity" $R/shear-np$np.txt $R/shear-ok.txt
done
for np in 4 16; do
  check "multilayer, sheared layers: $np ranks = serial" $R/multilayer/shear-np{1,$np}/out.txt
done
# Stokes wave with layered/nh.h: particles move vertically with w; with
# the prescribed flow, identical at 4 and 16 ranks (on 16 the column is on
# a subdomain corner). The solved flow is not bit-identical across rank
# counts (the nh.h pressure solver), so it is only run.
for np in 1 4 16; do
  run $R/nh/prescribed-after-np$np $np after/tnhp
  run $R/nh/solved-after-np$np $np after/tnh
done
for np in 4 16; do
  check "Stokes wave, prescribed flow: $np ranks = serial" $R/nh/prescribed-after-np{1,$np}/traj.txt
done

# Order in time: rotation with a time-dependent rate, linear in space
# (exact interpolation), error at t = 2 for dt = 0.02, 0.01, 0.005
mkdir -p $R/order
order() { # <description> <expected order> tconv args...
  local name=$1 o=$2; shift 2
  for dt in 0.02 0.01 0.005; do $B/after/tconv dt=$dt "$@" 2>/dev/null; done |
    grep '^[0-9]' > $R/order/$name.txt || true
  awk -v o=$o 'NR > 1 { q = log(e/$3)/log(2); if (q < o - 0.1 || q > o + 0.1) bad++ }
    { e = $3 } END { print (NR == 3 && !bad ? "ok" : "wrong") }' $R/order/$name.txt \
    > $R/order/$name.ok
  echo ok > $R/order/expected.ok
  check "order in time: $name, order $o" $R/order/$name.ok $R/order/expected.ok
}
order equal-steps 3
order unequal-steps 3 r=1.5
order solver-equal-steps 3 solver=1
order solver-unequal-steps 3 solver=1 r=1.5
order rk2 2 rk2=1
# centered.h: u(t) = G t, the position must be exact at t (old code: G t dt)
for dt in 0.01 0.005; do $B/after/tcentered $dt 2>/dev/null; done | grep '^centered' > $R/order/centered.txt || true
awk '{ if ($4 > 1e-12 || $4 < -1e-12) bad++ } END { print (NR == 2 && !bad ? "ok" : "wrong") }' \
  $R/order/centered.txt > $R/order/centered.ok
check "centered.h: particles use u(t), not u(t+dt)" $R/order/centered.ok $R/order/expected.ok

# Surface tracers (PARTICLE_SURFACE): sheared layers, particles placed in
# the bottom layer or above the surface must follow the top layer with z
# unchanged; the writer's z column must be the surface height eta
for np in 1 4 16; do
  run $R/multilayer/surface-np$np $np after/tmlsurf
  awk '/^shear/ {print $2}' $R/multilayer/surface-np$np/stdout.txt > $R/surface-np$np.txt
  check "surface tracers on $np ranks: top layer, z unchanged" $R/surface-np$np.txt $R/shear-ok.txt
done
for np in 4 16; do
  check "surface tracers: $np ranks = serial" $R/multilayer/surface-np{1,$np}/out.txt
done
for np in 1 4 16; do
  run $R/multilayer/surfacew-np$np $np after/tmlsurfw
done
d=$R/multilayer/surfacew-np1
cat $d/eta-* | sort -n > $d/eta.txt
awk 'NR == FNR {e[$1] = $2; next}
  FNR == 1 {for (i = 2; i <= NF; i++) if ($i != e[i - 2] + 0) bad++;
            print (NF == 4097 && bad == 0 ? "ok" : "wrong")}' \
  $d/eta.txt $d/w_particle_position_z.dat > $R/surface-eta.txt
check "surface tracers: written z = surface height" $R/surface-eta.txt $R/shear-ok.txt
for np in 4 16; do
  check "surface tracers, writer: $np ranks = serial" $R/multilayer/surfacew-np{1,$np}/w_particle_position_z.dat
done
# Quadtree (tree grids use the original all-gather exchange): the same
# flow on the same uniform cells must give the multigrid's trajectories
for np in 1 4; do
  run $R/accuracy/quadtree-np$np $np after/tmigq snap=10
done
check "quadtree, serial = multigrid serial" $R/accuracy/{quadtree,periodic/after}-np1/traj.txt
check "quadtree, 4 ranks = serial" $R/accuracy/quadtree-np{1,4}/traj.txt
# the writer on a quadtree with MPI, at odd steps (stage positions): the
# rank holding a particle must find its cell, or its velocity is nan
for np in 1 4; do
  run $R/writer/quadtree-np$np $np after/tmigqw np=40 tend=0.25 write=10 woff=1
done
d=$R/writer/quadtree-np4
if grep -q nan $d/w_particle_velocity_x.dat $d/w_particle_velocity_y.dat; then
  echo wrong; else echo ok; fi > $R/writer-quadtree-nan.txt
check "writer, quadtree, odd steps, 4 ranks: no nan velocity" $R/writer-quadtree-nan.txt $R/shear-ok.txt
for q in position_x velocity_x; do
  check "writer, quadtree, odd steps: $q, 4 ranks = serial" $R/writer/quadtree-np{1,4}/w_particle_$q.dat
done
check "box, 2 ranks = 8 ranks" $R/accuracy/box/after-np{2,8}/traj.txt
check "box, serial = 2 ranks" $R/accuracy/box/after-np{1,2}/traj.txt

# Tags (M4): the same particle has the same tag at any rank count
for np in 4 16; do
  check "tags, square: $np ranks = serial" $R/accuracy/periodic/after-np{1,$np}/tags.txt
  check "tags, four lists: $np ranks = serial" $R/accuracy/lists/after-np{1,$np}/tags.txt
done
awk '$2 != $3' $R/accuracy/periodic/after-np1/tags.txt > $R/tags-mismatch.txt
check "tags, square: tag = lattice index" $R/tags-mismatch.txt /dev/null
for v in after; do
  for np in 1 4 16; do
    run $R/tags/cells-$v-np$np $np $v/tmig cells=1 tend=0.5
  done
done
for np in 1 4 16; do
  run $R/tags/mlcells-np$np $np after/tmlc
  run $R/tags/release-np$np $np after/tmig release=50 tend=1
done
for np in 4 16; do
  check "tags, cells: $np ranks = serial" $R/tags/cells-after-np{1,$np}/tags.txt
  check "tags, multilayer cells: $np ranks = serial" $R/tags/mlcells-np{1,$np}/out.txt
  check "tags, releases: $np ranks = serial" $R/tags/release-np{1,$np}/tags.txt
  check "tags, releases: same counter on $np ranks" $R/tags/release-np{1,$np}/ntags.txt
done
awk '{print $1}' $R/tags/mlcells-np1/out.txt | sort -u | wc -l > $R/mlcells-unique.txt
wc -l < $R/tags/mlcells-np1/out.txt > $R/mlcells-count.txt
check "tags, multilayer cells: one tag per particle and layer" $R/mlcells-{unique,count}.txt

# Output (M5): 1,000,000 particles, 3 rows, gathering against MPI-IO
: > $R/output.txt
for np in 1 4 16 64; do
  d=$R/output/np$np
  run $d $np after/tmig np=1000 n=128 tend=0.01 out=3
  cat $d/output.txt >> $R/output.txt
  check "output, $np ranks: MPI-IO = gathered" $d/rows-{mpiio,gather}.bin
  [ $np = 1 ] || check "output, $np ranks = serial" $d/rows-mpiio.bin $R/output/np1/rows-mpiio.bin
  rm -f $d/rows-gather.bin $d/out.txt
done
rm -f $R/output/np{4,16,64}/rows-mpiio.bin

# Dump on 16 ranks, restore on 1, 4 and 64 (M5)
d=$R/restore/dump-np16
run $d 16 after/tmig np=200 tend=0.1 release=500 dump=dump.bin
for np in 1 4 64; do
  run $R/restore/after-np$np $np after/tmig restore=$d/dump.bin tend=0
  check "restore on $np ranks = dumped state" $R/restore/after-np$np/restored.txt $d/out.txt
done
check "restore keeps the tag counter" $R/restore/after-np4/ntags.txt $d/ntags.txt
e=$R/restore/release-np4
run $e 4 after/tmig restore=$d/dump.bin release=500 tend=0.2 # releases at steps 10 and 30
awk '{print $2}' $e/tags.txt | sort -n | uniq -d > $R/restore-duplicates.txt
check "releases after a restore: no duplicate tags" $R/restore-duplicates.txt /dev/null
echo "ntags $(( $(awk '{print $2}' $d/ntags.txt) + 1000 ))" > $R/restore-expected.txt
check "releases after a restore: tags continue from the counter" $e/ntags.txt $R/restore-expected.txt

# Positions in degrees on the sphere (PARTICLE_GEOGRAPHY), North Pacific
# domain: exact trajectories for uniform zonal and meridional flows
echo ok > $R/geo-ok.txt
for f in zonal meridional; do
  for np in 1 2 8; do
    d=$R/geography/$f-np$np
    run $d $np after/tsph flow=$f
    awk '/^geography/ {print $3}' $d/stdout.txt > $d/ok.txt
    check "geography, $f flow on $np ranks: exact trajectory" $d/ok.txt $R/geo-ok.txt
  done
  for np in 2 8; do
    check "geography, $f flow: $np ranks = serial" $R/geography/$f-np{1,$np}/out.txt
  done
done
# the writer with WCP_LON_WRAPPED: longitudes in [-180:180), files only
for np in 1 2 8; do
  run $R/geography/writer-np$np $np after/tsphw
done
d=$R/geography/writer-np1
awk -v w=$(awk '/^geography/ {print $6}' $d/stdout.txt) 'NR == 1 {
  for (i = 2; i <= NF; i++) { if ($i < -180 || $i >= 180) bad++; if ($i < 0) neg++ }
  print (bad == 0 && neg == w ? "ok" : "wrong") }' $d/w_particle_position_x.dat > $R/geo-wrapped.txt
check "geography, writer: longitudes wrapped in the files only" $R/geo-wrapped.txt $R/geo-ok.txt
for np in 2 8; do
  check "geography, writer: $np ranks = serial" $R/geography/writer-np{1,$np}/w_particle_position_x.dat
done

# Text tables (write_class_particles.h): 40,000 particles + 2 releases
# of 500, a row every 10 steps
: > $R/writer.txt
for v in after; do
  for np in 1 4 16; do
    d=$R/writer/$v-np$np
    run $d $np $v/tmigw np=200 tend=0.25 release=500 write=10
    echo "$v $(grep '^writer' $d/stdout.txt)" >> $R/writer.txt
  done
done
for np in 4 16; do
  for f in position_x velocity_y; do # no slip files for tracers
    check "writer: $f on $np ranks = serial" $R/writer/after-np{1,$np}/w_particle_$f.dat
  done
done
$PYTHON check_writer.py $R/writer/after-np4 | tee -a $R/checks.txt

echo "$(grep -c PASS $R/checks.txt) passed, $(grep -c FAIL $R/checks.txt) failed (list: $R/checks.txt)"
! grep -q FAIL $R/checks.txt
