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
| `MORPHL` (`libMORPHL.a`) | the library, `src/lib`: the kernels of FUSS's `src/lib/morph`, unchanged, plus `morph_domain_t` | `iso_fortran_env` (and `omp_lib`) only |
| `MORPH` (`bin/MORPH`) | the driver, `src/app`: reads a mesh, moves it for a number of steps, writes it | MORPHL, ORION (mesh I/O), FiNeR (`input.ini`) |

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
  base/ metrics/ motion/ quality/   FUSS's kernels, verbatim
  motion/Morph_Smooth.f90 the smoothed interior of the nodeshift law (new)
  domain/Morph_Domain.f90 morph_domain_t: the multi-block mesh state (new)
src/app/                the driver (MORPH)
  main.f90                read, move, write
  config/                 input.ini: the registry (FUSS's pattern) and Register_Morph
  io/Morph_IO.f90         the mesh through ORION: PLOT3D or Tecplot ASCII
test/
  CMakeLists.txt          the ctest suite
  unit/                   FUSS's three unit programs, verbatim, test_morph_smooth, and their runner
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
- `nodeshift` stops at its first step with `allocatable array is already allocated` (`Set_Shift` in
  `main.f90`). This is an ifx bug: an `allocate( law, source = <structure constructor> )` that omits
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
ctest --test-dir build -j 2              # 46 tests in an OpenMP build
ctest --test-dir build -L unit           # or -L guard, -L driver
./test/run_morph_guards.sh --build build # the guards without ctest, plus the unit programs built
                                         # standalone with gfortran -O2, as FUSS's script does
./test/check_build_config.sh             # build/ against CMakePresets.json
```

| label | tests | what they check |
|---|---|---|
| `unit` | 12 | `test_morph_block`, `test_morph_gcl` and `test_morph_nodeshift` (FUSS's) and `test_morph_smooth` (the smoothed interior, below) at 1 and 8 threads, plus the 1-thread and 8-thread outputs compared byte for byte |
| `guard` | 8 | `no-fuss`: no source uses a `FUSS_` module. `no-stop`: no `stop` or `error stop` in `src/lib`. `no-io`: no `print`, no `write`/`read` on a unit, no `open`/`close`/`inquire`/`flush` in `src/lib`. `lib-isolated`: `src/lib` uses only `iso_fortran_env`, `omp_lib` and `Morph_*`. `wall-clean`: the library compiles standalone and `gfortran -Wall` prints nothing. `wall-clean-app`: the same for the driver against the tree's modules. `build-config`: the tree matches the preset. `embed`: MORPHL in another CMake project (`test/embed`). |
| `driver` | 26 | `bin/MORPH` on a three-block warped mesh (`generate.py`: 4x6x4 + 6x2x4 + 6x4x4 cells, three conforming interfaces), one run per law at 1 and 8 threads, three more for the nodeshift variants (the uniform recession smoothed, a bumped recession smoothed, the same bump under the line law), then the two outputs of each compared, plus five runs that must be refused |

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
| 1 | an unknown key; a three-component key given two values |
| 2 | a nodeshift that consumes a whole grid line; a smoothing allowed one sweep (`smooth-maxit = 1`, status 6) -- MORPH statuses |
| 3 | a `gcl-tolerance` below round-off |

`guard_build-config` is skipped in a tree that is not the preset's `build/`. `wall-clean-app` drops one
gfortran 12 warning, `'.<name>' is used uninitialized`. It is a compiler false positive on the hidden
length of any deferred-length character array passed to an `intent(out)` dummy, as FiNeR's are. A
reproducer of 19 lines in two files, without FiNeR, shows it too. A Fortran name cannot start with
`.`, so no warning about the code itself can be dropped this way.

## The driver

From a case directory holding `input.ini`, run `<MORPH>/bin/MORPH`. Each run writes:

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
smooth-maxit     = 100000        ; smooth: sweeps allowed per set of nodes, then status 6 (MORPH_ERR_SMOOTH)
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
  return status 6 (`MORPH_ERR_SMOOTH`) with a message; the nodes are then left where they were.

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
