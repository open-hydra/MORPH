# MORPH

MORPH (Mesh mOtion, Remeshing and Positioning Handler) is the mesh-motion component of FUSS, as a
solver of its own. It works on structured multi-block meshes of hexahedra and does geometry only:

- moves the nodes of a block under a motion law (`static`, `prescribed`, `translation`, `nodeshift`);
  under `nodeshift` the nodes between the receding face and its opposite follow either their grid
  line or an elliptic smoothing of the face displacement (`interior`, below);
- rebuilds each cell's geometry: the signed volume (divergence theorem, four triangles per face), the
  metric tensor, the face normals and areas, and `dl`;
- returns the volume each face swept during the step, and the discrete geometric conservation law
  (GCL) residual `|(V_new - V_old) - sum_f dV_swept(f)|`, which is round-off (about 1e-16 relative);
- checks the mesh quality (an inverted cell, a Jacobian sign flip).

It has two parts:

| target | what | depends on |
|---|---|---|
| `MORPHL` (`libMORPHL.a`) | the library, `src/lib`: the kernels of FUSS's `src/lib/morph`, plus `morph_domain_t`, `Morph_Run` (the whole run in one call), the `given` law, the event list, and the wall as one surface across the blocks with its planes and fixed nodes | `iso_fortran_env` (and `omp_lib`) only |
| `MORPH` (`bin/MORPH`) | the driver, `src/app`: reads a mesh, moves it for a number of steps with one call to `Morph_Run`, writes it | MORPHL, ORION (mesh I/O), FiNeR (`input.ini`) |

The library never stops and does no I/O. Every entry point returns a `morph_status_t` and leaves the
decision to its caller. It uses no FUSS module. The guards check all three rules (see Test).

## Decisions

1. **A standalone copy, now FUSS's only one.** MORPH was copied from FUSS's `src/lib/morph`. FUSS
   now builds this repository instead, as its git submodule `lib/MORPH`, and keeps no copy of its own.
2. **A repository of its own.** MORPH is https://github.com/open-hydra/MORPH; FUSS pins a commit of it
   as its submodule `lib/MORPH`.
3. **Interior geometry only.** MORPH fills cells `1..dim` of each block. The ghost-cell metrics
   (`BC_Connect_Metrics`, `BC_Symmetry_Metrics` and `BC_Extrapolate_Metrics` in FUSS's `Lib_Metrics.f90`)
   need the block connectivity and stay in FUSS.

## Layout

```
CMakeLists.txt          MORPHL always; when top-level also ORION, FiNeR, the driver and the tests
install.sh              build / compile / update (FiNeR's dependencies)
cmake/                  FUSS's flag files, verbatim
src/lib/                the library (MORPHL)
  Morph_API.f90           the facade: the only module a solver needs to use
  base/ metrics/ motion/ quality/   FUSS's kernels; new: base/Morph_Events.f90, motion/Morph_Motion_Given.f90
  motion/Morph_Smooth.f90 the smoothed interior of the nodeshift law (new)
  core/Morph_Block_m.f90  Morph_Init and Morph_Update, moved out of the facade unchanged
  domain/Morph_Domain.f90 morph_domain_t: the multi-block mesh state, its boundaries and its wall (new)
  surface/                the wall as one surface across the blocks, its planes, the fixed-node masks (new)
  levelset/               the level-set kernels on a 2-D grid: the WENO-Z / Godunov / RK3 step and the
                          reinitialisation; the grid of blocks split among ranks, its ghosts and the
                          step on blocks (new; the regression that uses them comes later)
  parallel/               the exchange between ranks: the interface the library calls and its caller
                          implements (no MPI in the library), one rank's, and simulated ranks for the
                          tests (new)
  run/Morph_Run.f90       Morph_Run: the law installed, every step taken, checked and shown to an observer (new)
src/app/                the driver (MORPH)
  main.f90                read, call Morph_Run, write
  config/                 input.ini: the registry (FUSS's pattern) and Register_Morph
  io/Morph_IO.f90         the mesh through ORION: PLOT3D or Tecplot ASCII
  io/Morph_BC_IO.f90      the reader of the MORPH boundary-condition file, library MORPHBC (not read by the driver yet)
  io/Morph_Output.f90     the driver's observer: the GCL log, the step lines, the snapshots
test/
  CMakeLists.txt          the ctest suite
  unit/                   FUSS's three unit programs, verbatim, the new ones and their runner
  driver/                 run_case.sh and the multiblock case: generate.py, verify.py, input-*.ini
  run_morph_guards.sh     the architecture guards
  check_build_config.sh   the build tree against CMakePresets.json
```

## Build

Run `install.sh` from this directory:

```bash
git clone ../../../lib/third_party/FiNeR <somewhere>/FiNeR       # once: a FiNeR of MORPH's own (below)
./install.sh build --compilers=gnu --use-openmp --include-finer=<somewhere>/FiNeR
                       # wipes build/, fills FiNeR's dependencies, builds, writes CMakePresets.json
./install.sh compile   # incremental rebuild from CMakePresets.json
```

`CMakePresets.json` holds the absolute ORION and FiNeR paths of the last `install.sh build`, so it is
written by the build and ignored by git, as in FUSS. `test/check_build_config.sh` checks `build/`
against it.

`build` options are `--compilers=intel|gnu`, `--use-openmp`, `--use-tecio`, `--include-orion=<ORION tree>`
and `--include-finer=<FiNeR tree>`. Set `CMAKE_BUILD_PARALLEL_LEVEL=<n>` for a parallel build. There
is no MPI, because MORPH is rank-agnostic. By default ORION and FiNeR are hydra's
(`../../../lib/ORION` and `../../../lib/third_party/FiNeR`), the same layout ICE uses. The binary trees
of ORION and FiNeR are created inside MORPH's own build tree.

**Compilers.** MORPH is tested with gfortran 12.2.1. With ifx 2023.1 (`--compilers=intel`) it builds,
but 9 of its 31 tests fail and the 4 comparisons that depend on them do not run (the -Wall guards
are gfortran's and are not registered). The three reasons were found on 2026-10-24 and are not fixed:

- ORION writes PLOT3D values with list-directed output. gfortran writes 17 significant digits, which
  read back to the same bits; ifx writes 15 or 16. Every check of `verify.py` that compares the output
  bit for bit with the input or with the law then fails.
- `nodeshift` stopped at its first step with `allocatable array is already allocated` (`Set_Shift` in
  `main.f90` then; the allocation is now in `Morph_Run`'s law setup, `src/lib/run/Morph_Run.f90`, not
  tried again with ifx). This is an ifx bug: an `allocate( law, source = <structure constructor> )` that omits
  an allocatable component can leave that component with a non-null address, so `allocated()` says
  no and the next `allocate` says yes. The smoothed interior adds a second such component to the
  law, `dsp`, which its first step allocates after the same `allocated()` test; not checked with ifx.
- `prescribed` ends in a segmentation fault inside ifx's runtime, in the deallocation of the blocks
  by `domain%free`, under `-fiopenmp` only. Freeing each component of each block first avoids it; the
  cause is not known.

**FiNeR's dependencies, and why MORPH has its own FiNeR.** FiNeR v2.2.0 tracks none of its five
dependencies (PENF, StringiFor, FACE, FLAP, BeFoR64), although its CMake still adds them from
`src/third_party/`. Before configuring, `install.sh build` (or `install.sh update`) clones each
dependency into the FiNeR tree it is given, at the commit FUSS's `install.sh` pins. A dependency found
at another commit is **moved** to the pin, and the move is printed
(`finer_deps: <dep> moved <old> -> <new>`). hydra's own FiNeR holds other commits of all five, so
`./install.sh build` with the default FiNeR moves hydra's checkouts and changes what hydra's next
build compiles. Point `--include-finer` at a clone of your own, as above.

**hydra's checkout on 2026-10-24.** Its five FiNeR dependency checkouts were still the clones of
2026-05-25, at the commits FiNeR `cfc9194` tracked (PENF v1.3.13, StringiFor v1.1.7, FACE v1.1.12,
FLAP v1.2.13, BeFoR64 v1.1.8), not at the pins; hydra's own `figaro build` moves them to the same pins
(`scripts/utils/finer_deps.sh`). Its `main` was 9 commits behind `origin/main`, which pins ORION
`8cda645` where the checkout holds `5c32b73`: in the Fortran library `8cda645` changes the Tecplot ASCII
reader and adds the STRANDID it reads to `orion_data`. MORPH was accepted against `5c32b73`; a scratch
build against `8cda645` also passes all 33 tests and a Tecplot round trip.

ORION is only read. With `--use-tecio`, however, ORION builds TecIO inside its own source tree, so the
default is off: the driver reads and writes ASCII.

**MORPH as a library.** A solver can add this directory with `add_subdirectory` and link `MORPHL`.
In that case only `src/lib` is configured, with the parent's flags; ORION, FiNeR, the driver and the
tests are left out. The `.mod` files go to the parent's `CMAKE_Fortran_MODULE_DIRECTORY` when it sets
one, and to `<build>/…/src/lib/modules` when it does not; a target that links `MORPHL` finds them
either way. The `guard_embed` test builds such a parent both ways.

## Test

```bash
ctest --test-dir build -j 2              # 75 tests in an OpenMP build
ctest --test-dir build -L unit           # or -L guard, -L driver
./test/run_morph_guards.sh --build build # the guards without ctest, plus the unit programs built
                                         # standalone with gfortran -O2, as FUSS's script does
./test/check_build_config.sh             # build/ against CMakePresets.json
```

| label | tests | what they check |
|---|---|---|
| `unit` | 39 | `test_morph_block`, `test_morph_gcl` and `test_morph_nodeshift` (FUSS's), `test_morph_smooth` (the smoothed interior, below), `test_morph_events`, `test_morph_given`, `test_morph_run`, `test_morph_surface`, `test_morph_bc`, `test_morph_ls` (the level-set kernels, below), `test_morph_ls_blocks` (the level-set grid of blocks and its exchanges, below), `test_morph_ls_band` (its narrow band, below) and `test_morph_ls_reinit_blocks` (the reinitialisation on blocks, below) at 1 and 8 threads, plus the 1-thread and 8-thread outputs compared byte for byte. `test_morph_surface` and `test_morph_bc` write MORPH boundary files and read them with the driver's reader (MORPHBC); each run works in a directory of its own |
| `guard` | 8 | `no-fuss`: no source uses a `FUSS_` module. `no-stop`: no `stop` or `error stop` in `src/lib`. `no-io`: no `print`, no `write`/`read` on a unit, no `open`/`close`/`inquire`/`flush` in `src/lib`. `lib-isolated`: `src/lib` uses only `iso_fortran_env`, `omp_lib` and `Morph_*`. `wall-clean`: the library compiles standalone and `gfortran -Wall` prints nothing. `wall-clean-app`: the same for the driver against the tree's modules. `build-config`: the tree matches the preset. `embed`: MORPHL in another CMake project (`test/embed`). |
| `driver` | 28 | `bin/MORPH` on a three-block warped mesh (`generate.py`: 4x6x4 + 6x2x4 + 6x4x4 cells, three conforming interfaces), one run per law at 1 and 8 threads, three more for the nodeshift variants (the uniform recession smoothed, a bumped recession smoothed, the same bump under the line law), then the two outputs of each compared, plus seven runs that must be refused. A nodeshift run must also show its interior in `run.log`, and a smoothed one the smoothing's settings and its report under the step lines |

**On GitHub.** `.github/workflows/ctest.yml` runs the same build (`install.sh build --compilers=gnu
--use-openmp`, against ORION `5c32b73` and FiNeR `09adb98` cloned from GitHub) and all the tests above
on Ubuntu 24.04: for every pull request to `main`, every push to `main`, the merge queue, and by hand.
It fails when fewer than 75 tests are registered: without python3 >= 3.6, CMake drops the driver tests
with only a warning. Raise that floor when tests are added.

**The driver gate.** In every step of every run the discrete-GCL residual that MORPH logs must stay at
or below 1e-14 relative. FUSS's moving cases give ~1e-15; the measured maximum on this mesh is 4.1e-16
(gfortran 12.2.1, `-O3 -march=native -flto`, 1 and 8 threads), and up to 2.4e-15 with ifx 2023.1. A
real GCL defect does not hide at that level. In scratch builds, dropping one face's swept volume gave
2.0e-2, and the driver itself stopped. An error of 1e-9 relative in one face's swept volume gave
2.0e-11: the driver's own `gcl-tolerance` (1e-10) let it through, and this gate failed it. A
prescribed law evaluated at the start of the step instead of the end kept the GCL at 4.1e-16, and the
node check below caught it (1.0e-4 m).

`verify.py` then checks each law against the mesh MORPH read. Its bit-for-bit checks need the exact
PLOT3D round trip of a gfortran build (Compilers).

- `static`: the output is the input, bit for bit, and every residual is exactly 0.
- `translation`: `vel*dt` is added once per step, bit for bit.
- `prescribed`: the output matches `x_ref + amp*sin(k.x_ref + omega*t_N + phase)` within 1e-15 m. The
  measured difference is 5.6e-17 m, one ulp; the compiler contracts the argument into FMAs.
- `nodeshift`: the receding face moves by its displacement per step (`nodeshift-vel*dt`, bit for
  bit; with the bump of `nodeshift-bump-vel` within 1e-15 m, measured 0, because the compiler may
  contract the bump's arithmetic into FMAs). The opposite face and the other blocks do not move.
  Then, by `interior`:
  - `line` (`input-nodeshift.ini`, `input-nodeshift-bump-line.ini`): the lines in between are
    compressed uniformly, each on its own: they match a replica of the law within 1e-15 m (measured
    1.4e-17 m, one ulp, in both cases).
  - `smooth` (`input-nodeshift-smooth.ini`, `input-nodeshift-bump.ini`): the four edges match the
    replica of the line law (measured 0); the `y` and `z` coordinates, which no face node moves in,
    are untouched bit for bit; and, in every step, from the snapshots `out-diter = 1` writes, the
    solved nodes are in spring equilibrium: the residual `|sum_q w (d_q - d_p)| / sum_q w` with
    `w = 1/|e|` on the step's old nodes is at most 1e-11 of the largest face displacement (ten times
    `smooth-tolerance`; measured 2.7e-13 and 1.9e-13), and every node's displacement lies between the
    prescribed ones (the maximum principle; measured excess 0). On this warped mesh the smoothed
    nodes depart from the line law's `s/L` blend by up to 1.5e-7 m per step with the uniform
    displacement (the lines are not straight) and by up to 9.9e-4 m per step with the bump (the
    centre of the face moves 4.5e-3 m per step).
- for every law, the nodes shared by two blocks hold the same bits in both blocks.

The refused runs cover each exit status of the driver:

| exit status | refused run |
|---|---|
| 1 | an unknown key; a three-component key given two values; a `nodeshift-block` that names no block and a `smooth-omega` of 2 (both refused by `Morph_Run`) |
| 2 | a nodeshift that consumes a whole grid line; a smoothing allowed one sweep (`smooth-maxit = 1`, status 9, its message showing the `smooth-tolerance` of the input) -- MORPH statuses |
| 3 | a `gcl-tolerance` below round-off |

`guard_build-config` is skipped in a tree that is not the preset's `build/`. `wall-clean-app` drops one
gfortran 12 warning, `'.<name>' is used uninitialized`. It is a compiler false positive on the hidden
length of any deferred-length character array passed to an `intent(out)` dummy, as FiNeR's are. A
reproducer of 19 lines in two files, without FiNeR, shows it too. A Fortran name cannot start with
`.`, so no warning about the code itself can be dropped this way.

**The level-set kernels** (`src/lib/levelset`, `test_morph_ls`) are the first pieces of the surface
regression: one step of the level-set equation over a whole 2-D grid (fifth-order WENO-Z derivatives,
Godunov's upwind Hamiltonian, the three-stage strong-stability-preserving Runge-Kutta step), and the
safety reinitialisation with its trigger, the mean eikonal error near the front. Their order of
operations is fixed: built without fused multiply-adds (`-ffp-contract=off`), they give the same bits
as any implementation of the same formulas in the same order. The test uses grids of spacing
h = 1/32 with power-of-two origins and slopes, so that its exact checks are exact by construction.
Measured with gfortran 13.3 on aarch64, release flags:

- WENO-Z: on cubic data, the exact derivative within 0.58 eps of the largest one (bound 8 eps); next
  to a jump, the closed form of its weights within 1 eps (bound 32 eps); on flat data, 0; on data
  scaled by 2^20 or 2^-20, the value scaled, bit for bit (the epsilon scales with the data);
- Godunov's |grad phi|: on a linear field, exact at every node in both branches; at the bottom of a
  valley, 0 exactly for a front moving along +n and 1 - 2.08e-13 for one moving along -n (bound
  1e-12);
- one step: on a plane front, phi + dt; on a cubic profile, p(x + dt), which only a third-order step
  gives (a second-order one misses 5.96e-8); both exact in this build (bounds 4 and 8 eps of
  max |phi|). A zero velocity, +0 or -0, leaves phi untouched, bit for bit;
- a circle of radius 32 h, solid inside and fluid inside, over a stroke of h: the error at the nodes
  within 2h of the front is 1.50e-8 h and 1.75e-8 h, the same within 5.5e-6 (relative) for sub-steps
  of h/4, h/10 and h/100 (bounds 1e-6 h and 1e-3);
- the eikonal error: 0 on a linear field of unit slope, 1/32 on (x^2 - 1)/2, both exactly; the
  reinitialisation runs only when it is strictly above the threshold;
- the reinitialisation: the signed distance of a plane is a fixed point, bit for bit, and the -0 of
  an edge node keeps its sign; not run on the circle's distance (eikonal error 1.22e-4, threshold 0.03); on the circle's zero level with a
  gradient between 0.78 and 1.23, four iterations take the eikonal error from 0.0799 to 0.0226, the
  zero crossings along the grid edges move by 1.074e-3 h at most (bound 5e-3 h), and the eikonal
  error and the sum of phi recorded when the kernels were written come back within 1e-10 (relative;
  five gfortran builds spread them by 3.1e-16, a sum in another order could move them by 2.1e-12).

Each of the test's 40 checks fails under at least one of 67 deliberate changes of the kernels, and 66
of the 67 fail the test, before and after the reinitialisation's kernels were shared with the grid of
blocks (one of the changes had to be re-aimed: it took the correction's gradient from phi before the
step, which the shared kernel no longer receives, and now takes it from phi0). The other one cannot
change any result: the Godunov branch at the nodes
where V = 0, whose rate -V |grad phi| is a zero whose sign comes from V alone (|grad phi| >= 0).

**The level-set grid of blocks** (`levelset/Morph_LS_Grid.f90`, `parallel/`, `test_morph_ls_blocks`)
splits the same grid into blocks of nbn x nbn nodes (nbn >= 3), each block on the rank its caller
chooses. Before each stage of the step, the three ghost layers along each side of a block are filled
with the nodes of the block beyond it (copied when that block is on the same rank, swapped through
the exchange when it is on another) and, at the edge of the grid, with the extrapolation of the whole
grid; the corners are neither filled nor read. Every value a stencil reads is then the one it reads
on the whole grid, so the step on blocks gives the bits of the step on the whole grid, for any
blocks, any split among ranks and any number of threads. The library talks to the other ranks only
through `morph_exchange_t`: a swap of packed buffers with the contract of `MPI_Alltoallv`, an
allgather, and a sum of one partial per block taken in block order (the same bits for any number of
ranks). The caller implements it with its own MPI; `morph_exchange_serial_t` is one rank's.

The step is a collective: every rank checks its own arguments, then all share their verdicts and
whether V moves anything, in one allgather, so that either every rank leaves early or none does (a
rank that left alone would leave the others waiting in the next exchange). The test runs N ranks in
one process with `morph_exchange_sim_t`: each rank is a thread of a parallel region of N threads, and
a collective is: each rank writes its own mailbox, barrier, each reads all of them, barrier. Every
refusal of the simulated exchange is taken by all the ranks together: a team without one thread per
rank (`OMP_THREAD_LIMIT`, a call outside the region), a rank that is not its thread, counts or owners
that two ranks do not agree on. The library's own OpenMP loops, inside a rank's thread, are nested
regions, inactive by default. Built without OpenMP, the test skips the simulated ranks and says so.

Measured with gfortran 13.3 on aarch64, with the release flags (fused multiply-adds allowed) and with
`-ffp-contract=off`: in every case of the test the step on blocks gives the bits of the step on the
whole grid, at 1 and 8 threads (blocks of 3, 5 and 8 nodes; 1 to 20 blocks; 1 to 4 ranks, split as a
checkerboard, in halves, dealt in turn, with a rank that has no block, by rows, scattered). Under
`OMP_THREAD_LIMIT=2` the 9 checks that need 3 or 4 ranks fail, refused by the team check, and no rank
is left waiting. Each of the test's 50 checks fails under at least one of 120 deliberate changes of
the three new files, and 119 of the 120 fail the test (the 6 that remove a barrier of the simulated
exchange, in 20 runs out of 20). The other one cannot change any result: the zeroing of the step's
work array, whose corners no stencil reads.

**The narrow band** (`levelset/Morph_LS_Grid.f90`, `test_morph_ls_band`) keeps only some blocks active: they
hold phi and move, the others hold nothing. The active set is the same on every rank: all the blocks,
or the set given to `Morph_LS_Grid_Init`, then the set of each `Morph_LS_Grid_Update`, a collective in
which a block is active when at least one rank asks for it (the rank whose block the front approaches
asks for the blocks beyond it, whoever owns them). The ghosts toward an inactive neighbour come from
the block's own nodes and no message carries them: at every stage of a step `Morph_LS_Ghosts_Upwind`
takes them line by line from the speed, linear where the characteristics enter the block and quadratic
where they leave it (at the edge of the grid they stay linear, as on the whole grid). The update
returns, for each new slot, the slot its block had before (0 for a block just activated), and
`Morph_LS_Grid_Remap` moves any field of slots with it, leaving a NaN in the new slots for the caller to
fill. Every active node moves at its full speed. There is no window that slows the nodes far from the
front: a node whose speed is zero never changes, so it never comes back into the window, and the
levels ahead of a moving front pile up against it. The band then differs from the whole grid only by
the extrapolated ghosts at its edge and the values a block starts from; the caller chooses how far
from the front the inactive blocks stay, asks for the blocks and fills the new ones.

Why two kinds of ghosts: what the ghosts carry moves into the band where the characteristics enter
it (behind the front) and out of it where they leave (ahead of the front). In a Python emulation of
this scheme (240 x 240 nodes, blocks of 8 x 8, a block active when the front comes within 9.5 h of its
nodes, h the node spacing; times are the distance the front has moved, in h), on an expanding circle
after 39.9 h the band's front was 8.70e-5 h from the whole grid's with linear ghosts, 1.61e-6 h with
exact ghosts where the characteristics leave and linear ones where they enter, and still 8.70e-5 h
the other way round. Quadratic ghosts on every side brought it to 5.34e-7 h, but elsewhere the band
drifted 8100 h away from the whole grid (1330 h after 60 h on two merging circles). The switched
ghosts brought the front 30.4 to 109.5 times closer to the whole grid's than linear ones (an
expanding and a shrinking circle after 39.9 h and 75 h, an expanding circle with a speed from -1 to
1/2 after 39.9 h) and 2.8 times on the merging circles after 60 h, where the side behind the front
matters too; the largest distance elsewhere was the same as with linear ones, to three digits. With
the switched ghosts, blocks activated within 18.5 h of the front, filled with the exact solution and
dropped beyond 26.5 h, the band's front stayed within 4.3% of the whole grid's own error against the
exact solution every 6 h and at the end, on the expanding and the shrinking circle (to 39.9 h) and on
the merging circles (to 54 h); with linear ghosts no margin up to 24.5 h did (24.6% on the merging
circles at 24.5 h).

The test checks the band against an oracle written in the test on the whole grid (the ghosts of a
side from the whole grid when the block beyond is active, else extrapolated as above; then the
kernels' rate and stage), itself the step on the whole grid when every block is active, and
`Morph_LS_Ghosts_Upwind` alone on lines of integers, whose ghosts are exact in any build.

Measured with gfortran 13.3 on aarch64, with the release flags (fused multiply-adds allowed) and with
`-ffp-contract=off`: in every case of the test the band gives the oracle's bits, at 1 and 8 threads, and
the oracle with every block active gives the bits of the step on the whole grid. Each of the test's 62
checks but one fails under at least one of 72 deliberate changes of the band's code (25 in
`Morph_LS_Ghosts_Upwind`, 9 in the ghost fill it shares with the grid of blocks), and all 72 fail the
test: 69 with a failed check, 2 by a crash (an allocation not released before the next, a field read
before it is allocated) and 1 by the time limit (a rank that does not take another rank's refusal
leaves the others waiting). The other check is about the test's own fields: that the edges of two of
its bands hold lines where the characteristics enter and lines where they leave; it fails, alone, when
their speed is zero.

**The reinitialisation on blocks** (`levelset/Morph_LS_Reinit.f90`, `test_morph_ls_reinit_blocks`)
runs the reinitialisation and its trigger on the grid of blocks, every block active or a band. The
computations at the nodes are kernels on one block with its three ghost layers and a flag for each
side: a side that faces the edge of the grid or an inactive block is an edge, where the gradients are
one-sided and no node is cut, as on the edges of the whole grid. The routines on the whole grid call
the same kernels on one block with four edges. On blocks, `Morph_LS_Halo_Fill` fills the ghosts
without V (the nodes of the block beyond, linear toward an inactive block and at the edge of the
grid): those of phi0 once, then at every step those of phi and, with the correction, those of the
step; phi is assigned at the end only, so a refusal or a failed exchange leaves it untouched on every
rank. With every block active the reinitialisation on blocks gives the bits of the whole grid's, for
any blocks, any split among ranks and any number of threads. On a band, a node's value after k steps
depends only on the nodes within L1 distance 4k (WENO-Z reads 3 nodes along each axis, the
correction one more): the nodes farther than that from every inactive node get the bits of the whole
grid started from the same phi. When the band has been updated to the blocks within M of the front
just before, a node within 2 h of the front is farther than M - 2 h from every inactive node (L1 is at
least the distance), so it keeps those bits while 4 iters <= M/h - 2: with the band's margin of
18.5 h, up to 4 iterations, the default. These are distances, not values of phi: where |grad phi| is
g rather than 1, a band and nodes chosen by |phi| reach M/g and 2 h/g, and the bound becomes
4 iters <= (M/h - 2)/g.

The measure sums each block in node order, then the blocks in block order with the exchange's
ordered sum: the same bits on every rank, for any split and any number of ranks and threads, and on
one block covering the grid the bits of the whole grid. On several blocks the additions come in
another order than on the whole grid, which moves the last bits (within 2 (n+1) u E, n the nodes in
the band, u = 2^-53). Every rank checks its own arguments and all share their verdicts in one
allgather; iters, beta and the trigger's threshold are compared between the ranks bit for bit, since
a rank given other iters would make another number of exchanges and leave the others waiting.

In a Python emulation of these kernels on the cases of the narrow band (240 x 240 nodes, blocks of
8 x 8 active within 18.5 h of an expanding circle after 39.9 h, 304 of 900, of a shrinking one after
39.9 h, 180, and of two merging circles after 36 h, 364; phi the exact distance times a factor from 1
to 2, so that the reinitialisation has work to do), 4 and 6 iterations on the band and on the whole
grid from the same phi: no node farther than 4 iters from the inactive nodes differed in its bits, the
farthest that differed was at L1 distance 8 and 12, and none of the 676 to 1370 nodes with
|phi| <= 2 h differed; E on the band had the bits of E on every block. Built with `-ffp-contract=off`,
the Fortran code gives the emulation's bits in four cases (scattered blocks of 5 with and without the
correction, every block of 4, a checkerboard of blocks of 3), at 1 and 8 threads.

The test checks, bit for bit: every block active against the whole grid (blocks of 3, 5 and 8; 4
steps with the correction, 3 without, 0); bands (one block alone, a checkerboard, a frame, scattered
blocks, a ring) against an oracle written in the test on the whole grid, which gives each block the
ghosts the band reads and calls the public kernels; the reach after 2 and 4 steps; the measure against
the oracle's sum and the whole grid; 2 to 4 simulated ranks, one of them with no active block,
against one; the trigger below, at and above its threshold; the refusals, phi untouched (a grid whose
Init was refused among them, on one rank and on one of two); an exchange that fails at each of its
collectives, and the number of collectives of each kind. One check has a
single node in the band, at a corner of the grid, so that E is that node's term alone: across a
linearly extrapolated ghost the central difference is the one-sided one in exact arithmetic, and the
sums of the other checks absorb the rounding by which they differ. Measured with gfortran 13.3 on
aarch64, with the release flags and with `-ffp-contract=off`: its 50 checks pass at 1 and 8 threads,
with the same output (37 without OpenMP, the simulated ranks skipped). Each of the 50 fails under at
least one of 128 deliberate changes of `Morph_LS_Reinit`. The test catches 84 of them (77 with a
failed check, 4 by a crash, 3 by the time limit: ranks that do not compare their iters, or a rank
that ignores the others' refusals, make different numbers of exchanges and wait for each other), 67
that `test_morph_ls` does not; the two together catch 124. The other 4 cannot change any result: each drops the setting of values that
nothing reads (the ghosts of the measure and of phi0 on the whole grid, whose sides are all edges;
the ghost corners of phi0 on blocks; the distance of the subcell fix on edge nodes, which are never
cut), and each passes the whole suite.

## The driver

From a case directory holding `input.ini`, run `<MORPH>/bin/MORPH`. The computation is one library
call, `Morph_Run` (`src/lib/run/Morph_Run.f90`): the driver reads `input.ini` and the mesh, makes the
call, and writes what its observer is shown after each step and the mesh at the end. A solver that
holds a `morph_domain_t` can make the same call. Each run writes:

- `mesh-out`, the mesh at the end;
- every `out-diter` steps, `mesh-out` with `_<step>` before the extension;
- `gcl-log`, one row per step: the step, the time, the absolute and the relative GCL residual, the
  block where the relative one is largest, and the largest volume fraction a cell swept out.

The exit status is 0 on success, 1 for input errors, 2 when MORPH returned an error status and 3 when
a step's GCL residual exceeded `gcl-tolerance`.

```ini
[MORPH-IO]
mesh-in     = INPUT/mesh.p3d     ; required
mesh-out    = OUTPUT/mesh.p3d
mesh-format = p3d                ; p3d (17 digits and exact with gfortran) or tecplot (ASCII .tec/.dat, 15 digits)
gcl-log     = OUTPUT/gcl.dat
out-diter   = 0                  ; 0: the final mesh only
shell-diter = 1

[MORPH-Time]
t0    = 0.0
dt    = 1.0e-3                   ; required
steps = 20                       ; required

[MORPH-Motion]
law = prescribed                 ; static | prescribed | translation | nodeshift
amp = 2.0e-3 1.5e-3 1.0e-3       ; prescribed: x_c moves to xr_c + amp_c sin(kx_c xr + ky_c yr + kz_c zr + omega t + phase_c)
kx  = 10.0 5.0 8.0
ky  = 7.0 12.0 3.0
kz  = 4.0 6.0 9.0
omega = 50.0
taper-to-boundary = .false.
vel = 0.0 0.0 0.0                ; translation [m/s]
nodeshift-block = 1              ; nodeshift: this block's face moves, every other block is static
nodeshift-face  = 1              ; 1..6 = i-low, i-high, j-low, j-high, k-low, k-high
nodeshift-vel   = 2.5 0.0 0.0    ; [m/s], every node of the face
nodeshift-bump-vel   = 0.0 0.0 0.0  ; [m/s], added at the centre of the face times (1 - r^2/width^2)^2, r < width
nodeshift-bump-width = 0.25      ; r over the face mapped to the unit square
interior         = line          ; nodeshift: line | smooth, the nodes between the receding and the pinned face
smooth-tolerance = 1.0e-12       ; smooth: a sweep's largest update over the largest face displacement
smooth-maxit     = 100000        ; smooth: sweeps allowed per set of nodes, then status 9 (MORPH_ERR_SMOOTH)
smooth-omega     = 0.0           ; smooth: over-relaxation, below 2; 0 = 2/(1 + sin(pi/N))
gcl-tolerance   = 1.0e-10        ; asserted every step, as FUSS's Update_Mesh does
```

`[MORPH-Motion]` keeps the names of FUSS's `[FUSS-MeshMotion]` keys. The `nodeshift-*` keys are the
driver's own: in FUSS the receding wall gives the displacement. `[MORPH-IO]` and `[MORPH-Time]` belong to
the driver.

### The nodeshift interior

Under `nodeshift` one face of a block moves by a displacement given per node (FUSS: the ablating
wall's recession over the step), the opposite face is pinned, and `interior` says how the nodes
between them move. Both keep the receding face, the pinned face and the line-consumption check
(status `MORPH_ERR_NEGVOL` when a line is shorter than its displacement) the same.

- `line` (the default, `Morph_Motion_Nodeshift.f90`): a node at arc length `s` from the pinned end
  of its line, of length `L`, moves by `disp * s/L`. Each line is compressed uniformly on its own, a
  graded mesh keeps its grading, and a displacement that varies across the face is imprinted, scaled
  by `s/L`, on every node plane down to the pinned face. This is FUSS's law since plan 08: the
  default stays `line` so that FUSS's regression references hold bit for bit. The path's node
  arithmetic is unchanged by the smoothing (only the search for the consumed line moved into a
  procedure of its own): a build of the commit before it, `ae699fc`, writes the same bytes for the
  `nodeshift` driver case and the `test_morph_nodeshift` output (checked 2026-10-27).
- `smooth` (`Morph_Smooth.f90`): the displacement of every node between the faces is the equilibrium
  of lineal springs on the grid edges, stiffness `1/|e|` on the old nodes (Batina 1990): at each solved
  node `sum_q (d_q - d_p)/|e_pq| = 0`, the discrete Laplacian with inverse-distance weights. Along a
  straight line that makes `d` linear in the arc length, which is the line law, so on a block of
  straight lines with a uniform face displacement the two interiors agree to the solver's tolerance
  (`test_morph_smooth`: 0 and 1.4e-17 m on a Cartesian and on a graded block), and a graded line keeps
  its grading; index-space weights would not do that. With a displacement that varies across the face
  the deformation spreads into the block and dies away towards the pinned face instead of being
  imprinted on every plane. The solve runs lowest dimension first, so that an entity is fixed before
  the next is solved: the receding face (`disp`), the pinned face (0), the four edges of the block
  that join them (the line law's closed form, evaluated the same way in every block that holds the
  edge; against the line path's own arithmetic one rounding can differ, the compiler contracting that
  path into FMAs: 1.4e-17 m on 4 of the 20 edge nodes of the test mesh), the four side faces (the
  in-face springs only, Dirichlet on the edges), then the interior (all six springs). Each set is iterated with red-black
  SOR until a sweep's largest update is at most `smooth-tolerance` times the largest face
  displacement, `smooth-maxit` sweeps at most, with `smooth-omega` or, at 0, Young's
  `2/(1 + sin(pi/N))` for the longest extent `N` of the set. A node of one colour reads only the other
  colour and the stopping test is a max-norm, so the result does not depend on the number of threads
  (the `_threads_identical` tests). The start of each step's solve is the previous step's field (the
  line law on the first step): on the test mesh a steady recession converges in about 15-20 interior
  sweeps and 12 per side face. A set that does not converge, or two coincident neighbouring nodes,
  return status 9 (`MORPH_ERR_SMOOTH`) with a message; the nodes are then left where they were.

Two consequences of the smoothed interior to keep in mind:

- **2-D meshes and block seams.** In a 2-D mesh (two node planes in `k`) every node of a seam between
  two blocks is on an edge of each block, and the edges are the line law's closed form: the two
  blocks agree on the seam bit for bit, as FUSS's `Update_Mesh` asserts. The two node planes are
  solved as two side faces with the same data, weights and colouring and come out the same bit for
  bit (`test_morph_smooth`). In a 3-D mesh a seam face has nodes inside it, solved by each block
  with its own sweep order: the two blocks agree there only to the tolerance, not bitwise. FUSS's
  interface assertion would then refuse the step; a 3-D use needs the seam synchronised by the
  caller (or the seam nodes made Dirichlet) first.
- **Restarts.** The warm start makes the field reproducible across a restart only to the
  tolerance (a restart solves the first step from the line law). Nothing else is kept: the law
  still displaces incrementally from the old nodes and needs nothing beyond the mesh.

- Only `[MORPH-*]` sections are read. In those sections the driver refuses an unknown key, a repeated
  section, a value that cannot be read, a vector key without exactly 3 values and a string longer
  than its variable.
- As in every hydra INI file, no comment line contains `=`.
- A 2-D mesh must have two node planes in `k`, as FUSS's meshes do; a single plane is refused.
- `prescribed` takes the mesh as read for its reference `(xr, yr, zr)`. Given a mesh that has already
  been moved, it displaces it a second time; FUSS refuses a restart with this law for that reason,
  and the driver cannot tell the two meshes apart. Because the law is anchored at the reference mesh,
  the first step jumps from the reference to the displaced shape (a swept-outflow ratio of 0.12 on the
  test mesh); later steps move by about `amp*omega*dt`.

## Where MORPH comes from

The first commit (`40f4849`) copies 14 files from FUSS commit
`2943354329820c795d715323a42eace656c58924` (branch `worktree-porous-pyrolysis`, 2026-10-23, which is
on FUSS's origin), byte for byte: `src/lib/morph/` becomes `src/lib/`, and `src/lib/morph/test/`
becomes `test/unit/`. The copied files then had the same git blob hashes in both repositories. From
the same FUSS commit MORPH later took the three `cmake/` flag files, also verbatim, and adapted
these:

- the build: the LTO archiver block of FUSS's `CMakeLists.txt:44-71`, `install.sh`,
  `test/check_build_config.sh` and the guard script;
- the driver's input layer: FUSS's `Registry.f90`, `Backend_INI.f90` and `Read_Ini.f90`, renamed.

**FUSS builds this repository.** FUSS has no `src/lib/morph` any more: its `CMakeLists.txt` adds this
repository with `add_subdirectory` (the submodule `lib/MORPH`, or the tree `MORPH_PATH` names), which
configures `MORPHL` only, with FUSS's flags, and `FUSSL` links it. FUSS still calls MORPH through
`src/lib/numerics/space/Adapter_Morph.f90` alone, and its `test/run_morph_guards.sh` runs this
repository's guards on the submodule. An edit here reaches FUSS when FUSS moves its submodule pin; the
drift guard that compared the two copies (`test/check_sync.sh`) is gone with FUSS's copy.

## Next steps

1. **The duplicated-kernel cleanup.** FUSS's `Lib_Metrics.f90` still holds its own copies of the
   kernels: `Check_Mesh_Type`, `Compute_Norm_Area` and `Compute_Metric_Tensor`. The last one is live:
   the three `BC_*_Metrics` call it from 18 sites. The two copies already
   disagree on the ghost volume (`abs()` five-tetrahedron against MORPH's signed volume), and they
   agree on `M`, the only ghost metric that is read. The ghost metrics stay in FUSS (decision 3), so
   the cleanup is no longer a port of `BC_*_Metrics` into MORPH. Instead:
   - make FUSS's `BC_*_Metrics` call MORPH's per-cell kernel. `Morph_Metric_Tensor` is public in
     `Morph_Metrics` but not exported by `Morph_API`, so the API must export it first;
   - then delete the three legacy kernels.

   The gate: FUSS's `run_regression.sh -p 1 -k` unchanged, where only `volg`, which
   nothing reads, may move.
2. **Who owns `blk%node_old`.** At FUSS `2943354`, FUSS allocates `blk%node_old` on a moving mesh
   (`Mod_Allocate_Data.f90:118`) and frees it on remote blocks (`:336`). The adapter writes it
   (`Adapter_Morph.f90:233`), but it builds its own previous positions from the current nodes
   (`:179`), and nothing else in FUSS reads `blk%node_old`. In MORPH's `morph_domain_t` the previous
   positions belong to MORPH: `Morph_Update` copies `node` into `node_old` before applying the law.
   The options are:
   - FUSS holds a `morph_domain_t`, and `blk%node_old` goes;
   - FUSS keeps passing its own arrays to `Morph_Update` one block at a time, as it does today.

   FUSS's switch to this repository left it open. The first option waits on the ifx crash in
   `domain%free` (Compilers): hydra builds FUSS with ifx too (`figaro build hydra-MF --compilers=intel`).

Smaller items noticed during the copy:

- The driver's registry is a copy of FUSS's `Registry.f90`, and ICE carries a registry of its own: a
  fix to one reaches none of the others.
- `generate_markdown` in the registry has no caller yet: FUSS documents its keys through
  `src/app/docgen.f90`, MORPH through this README.
- `amp`'s `>= 0` rule is FUSS's, and, as in FUSS, the registry does not check rules on vector keys.

## License

MORPH is free and open-source software released under the [GNU General Public License v3.0](LICENSE), as FUSS is.
