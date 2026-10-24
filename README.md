# MORPH

MORPH (Mesh mOtion, Remeshing and Positioning Handler) is the mesh-motion component of FUSS, as a
solver of its own. It works on structured multi-block meshes of hexahedra and does geometry only:

- moves the nodes of a block under a motion law (`static`, `prescribed`, `translation`, `nodeshift`);
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

## Decisions (the user, 2026-10-24)

The record is FUSS's `plan-bucket/12-coupled-wall-follow-up.md` (DECISIONS 2026-10-24) and the note in
plan 09 §8.

1. **A standalone copy.** MORPH was copied from FUSS's `src/lib/morph`. FUSS keeps its own copy,
   unchanged, and still builds from it. FUSS switches to this repository later (Next steps).
2. **A local repository only.** There is no remote and nothing is pushed. MORPH sits in hydra's
   checkout at `src/lib/MORPH` but is not part of hydra: hydra shows it as untracked (`?? src/lib/MORPH/`),
   and that is intended.
3. **Interior geometry only.** MORPH fills cells `1..dim` of each block. The ghost-cell metrics
   (`BC_Connect_Metrics`, `BC_Symmetry_Metrics` and `BC_Extrapolate_Metrics` in FUSS's `Lib_Metrics.f90`)
   need the block connectivity and stay in FUSS (plan 09 §7).

## Layout

```
CMakeLists.txt          MORPHL always; when top-level also ORION, FiNeR, the driver and the tests
install.sh              build / compile / update (FiNeR's dependencies)
cmake/                  FUSS's flag files, verbatim
src/lib/                the library (MORPHL)
  Morph_API.f90           the facade: the only module a solver needs to use
  base/ metrics/ motion/ quality/   FUSS's kernels, verbatim
  domain/Morph_Domain.f90 morph_domain_t: the multi-block mesh state (new)
src/app/                the driver (MORPH)
  main.f90                read, move, write
  config/                 input.ini: the registry (FUSS's pattern) and Register_Morph
  io/Morph_IO.f90         the mesh through ORION: PLOT3D or Tecplot ASCII
test/
  CMakeLists.txt          the ctest suite
  unit/                   FUSS's three unit programs, verbatim, and their runner
  driver/                 run_case.sh and the multiblock case: generate.py, verify.py, input-*.ini
  run_morph_guards.sh     the architecture guards
  check_build_config.sh   the build tree against CMakePresets.json
  check_sync.sh           the drift guard against FUSS (not a ctest test)
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
is no MPI, because MORPH is rank-agnostic (plan 09 §9). By default ORION and FiNeR are hydra's
(`../../../lib/ORION` and `../../../lib/third_party/FiNeR`), the same layout ICE uses. The binary trees
of ORION and FiNeR are created inside MORPH's own build tree.

**FiNeR's dependencies, and why MORPH has its own FiNeR.** FiNeR v2.2.0 tracks none of its five
dependencies (PENF, StringiFor, FACE, FLAP, BeFoR64), although its CMake still adds them from
`src/third_party/`. Before configuring, `install.sh build` (or `install.sh update`) clones each
dependency into the FiNeR tree it is given, at the commit FUSS's `install.sh` pins. A dependency found
at another commit is **moved** to the pin, and the move is printed
(`finer_deps: <dep> moved <old> -> <new>`). hydra's own FiNeR holds other commits of all five, so
`./install.sh build` with the default FiNeR moves hydra's checkouts and changes what hydra's next
build compiles. Point `--include-finer` at a clone of your own, as above.

ORION is only read. With `--use-tecio`, however, ORION builds TecIO inside its own source tree, so the
default is off: the driver reads and writes ASCII.

**MORPH as a library.** A solver can add this directory with `add_subdirectory` and link `MORPHL`.
In that case only `src/lib` is configured, with the parent's flags; ORION, FiNeR, the driver and the
tests are left out. The `.mod` files go to the parent's `CMAKE_Fortran_MODULE_DIRECTORY` when it sets
one, and to `<build>/…/src/lib/modules` when it does not; a target that links `MORPHL` finds them
either way. The `guard_embed` test builds such a parent both ways.

## Test

```bash
ctest --test-dir build -j 2              # 33 tests in an OpenMP build
ctest --test-dir build -L unit           # or -L guard, -L driver
./test/run_morph_guards.sh --build build # the guards without ctest, plus the unit programs built
                                         # standalone with gfortran -O2, as FUSS's script does
./test/check_build_config.sh             # build/ against CMakePresets.json
```

| label | tests | what they check |
|---|---|---|
| `unit` | 9 | `test_morph_block`, `test_morph_gcl` and `test_morph_nodeshift` (FUSS's) at 1 and 8 threads, plus the 1-thread and 8-thread outputs compared byte for byte |
| `guard` | 8 | `no-fuss`: no source uses a `FUSS_` module. `no-stop`: no `stop` or `error stop` in `src/lib`. `no-io`: no `print`, no `write`/`read` on a unit, no `open`/`close`/`inquire`/`flush` in `src/lib`. `lib-isolated`: `src/lib` uses only `iso_fortran_env`, `omp_lib` and `Morph_*`. `wall-clean`: the library compiles standalone and `gfortran -Wall` prints nothing. `wall-clean-app`: the same for the driver against the tree's modules. `build-config`: the tree matches the preset. `embed`: MORPHL in another CMake project (`test/embed`). |
| `driver` | 16 | `bin/MORPH` on a three-block warped mesh (`generate.py`: 4x6x4 + 6x2x4 + 6x4x4 cells, three conforming interfaces), one run per law at 1 and 8 threads, then the two outputs compared, plus four runs that must be refused |

**The driver gate.** In every step of every run the discrete-GCL residual that MORPH logs must stay at
or below 1e-14 relative. Plan 09 records ~1e-15 for FUSS's moving cases; the measured maximum on this
mesh is 4.1e-16 (gfortran 12.2.1, `-O3 -march=native -flto`, 1 and 8 threads). A real GCL defect does
not hide at that level. In scratch builds (FUSS `plan-bucket/records/2026-10-24-morph-standalone/red/`),
dropping one face's swept volume gave 2.0e-2, and the driver itself stopped. An error of 1e-9
relative in one face's swept volume gave 2.0e-11: the driver's own `gcl-tolerance` (1e-10) let it
through, and this gate failed it. A prescribed law evaluated at the start of the step instead of the
end kept the GCL at 4.1e-16, and the node check below caught it (1.0e-4 m).

`verify.py` then checks each law against the mesh MORPH read:

- `static`: the output is the input, bit for bit, and every residual is exactly 0.
- `translation`: `vel*dt` is added once per step, bit for bit.
- `prescribed`: the output matches `x_ref + amp*sin(k.x_ref + omega*t_N + phase)` within 1e-15 m. The
  measured difference is 5.6e-17 m, one ulp; the compiler contracts the argument into FMAs.
- `nodeshift`: the receding face moves by `nodeshift-vel*dt` per step, bit for bit. The opposite face
  and the other blocks do not move. The lines in between are compressed uniformly: they match a
  replica of the law within 1e-15 m (measured 1.4e-17 m).
- for every law, the nodes shared by two blocks hold the same bits in both blocks.

The refused runs cover each exit status of the driver:

| exit status | refused run |
|---|---|
| 1 | an unknown key; a three-component key given two values |
| 2 | a nodeshift that consumes a whole grid line (a MORPH status) |
| 3 | a `gcl-tolerance` below round-off |

`guard_build-config` is skipped in a tree that is not the preset's `build/`. `wall-clean-app` drops one
gfortran 12 warning, `'.<name>' is used uninitialized`. It is a compiler false positive on the hidden
length of any deferred-length character array passed to an `intent(out)` dummy, as FiNeR's are; a
15-line reproducer without FiNeR shows it too. A Fortran name cannot start with `.`, so no warning
about the code itself can be dropped this way.

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
mesh-format = p3d                ; p3d (exact, 17 digits) or tecplot (ASCII .tec/.dat, 15 digits)
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
gcl-tolerance   = 1.0e-10        ; asserted every step, as FUSS's Update_Mesh does
```

`[MORPH-Motion]` keeps the names of FUSS's `[FUSS-MeshMotion]` keys. The `nodeshift-*` keys are the
driver's own: in FUSS the receding wall gives the displacement. `[MORPH-IO]` and `[MORPH-Time]` belong to
the driver.

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

## Where MORPH comes from, and the sync guard

The first commit (`40f4849`) copies 14 files from FUSS commit
`2943354329820c795d715323a42eace656c58924` (branch `worktree-porous-pyrolysis`, 2026-10-23, which is
on FUSS's origin), byte for byte: `src/lib/morph/` becomes `src/lib/`, and `src/lib/morph/test/`
becomes `test/unit/`. The copied files then had the same git blob hashes in both repositories. From
the same FUSS commit MORPH later took the three `cmake/` flag files, also verbatim, and adapted
these:

- the build: the LTO archiver block of FUSS's `CMakeLists.txt:44-71`, `install.sh`,
  `test/check_build_config.sh` and the guard script;
- the driver's input layer: FUSS's `Registry.f90`, `Backend_INI.f90` and `Read_Ini.f90`, renamed.

**FUSS still uses its own copy.** It compiles `src/lib/morph` into `FUSSL`, calls it through
`src/lib/numerics/space/Adapter_Morph.f90`, and links nothing from this repository. So an edit here
does not reach FUSS, and an edit in FUSS does not reach MORPH. To see how far the two have drifted,
fetch FUSS, then run:

```bash
git -C ../FUSS fetch
./test/check_sync.sh [--fuss <FUSS repository>] [--diff] [REF]   # REF default: worktree-porous-pyrolysis
```

The guard reads FUSS through git only and changes neither repository. It reports three things:

- which verbatim files MORPH has changed since `2943354`;
- which of FUSS's originals, verbatim or adapted, FUSS has changed between `2943354` and REF;
- which files exist under FUSS's `src/lib/morph` at REF but not in MORPH.

It exits 1 when anything differs. It is not part of ctest, because FUSS moving is not a MORPH defect.

## Next steps (plan 09 §8)

1. **FUSS's switch to MORPH.** Replace FUSS's in-tree `src/lib/morph` with an `add_subdirectory` of
   this repository, discovered the way ORION is (`MORPH_PATH`), and link `MORPHL`. `Adapter_Morph.f90`
   is the only FUSS file that changes (plan 09 §7). After the switch, FUSS's regression baseline and
   `run_morph_guards.sh` must be unchanged. Run `check_sync.sh` first: a FUSS-side change still to
   port must reach MORPH before FUSS starts building from it.
2. **The duplicated-kernel cleanup.** FUSS's `Lib_Metrics.f90` still holds its own copies of the
   kernels: `Check_Mesh_Type`, `Compute_Norm_Area` and `Compute_Metric_Tensor`. The last one is live:
   the three `BC_*_Metrics` call it from 18 sites. According to plan 09 §11, the two copies already
   disagree on the ghost volume (`abs()` five-tetrahedron against MORPH's signed volume), and they
   agree on `M`, the only ghost metric that is read. The ghost metrics stay in FUSS (decision 3), so
   the cleanup is no longer the port of `BC_*_Metrics` into MORPH that plan 09 §8 planned. Instead:
   - make FUSS's `BC_*_Metrics` call MORPH's per-cell kernel. `Morph_Metric_Tensor` is public in
     `Morph_Metrics` but not exported by `Morph_API`, so the API must export it first;
   - then delete the three legacy kernels.

   The gate is plan 09 §11's: FUSS's `run_regression.sh -p 1 -k` unchanged, where only `volg`, which
   nothing reads, may move.
3. **Who owns `blk%node_old`.** At FUSS `2943354`, FUSS allocates `blk%node_old` on a moving mesh
   (`Mod_Allocate_Data.f90:118`) and frees it on remote blocks (`:336`). The adapter writes it
   (`Adapter_Morph.f90:233`), but it builds its own previous positions from the current nodes
   (`:179`), and nothing else in FUSS reads `blk%node_old`. In MORPH's `morph_domain_t` the previous
   positions belong to MORPH: `Morph_Update` copies `node` into `node_old` before applying the law.
   The options are:
   - FUSS holds a `morph_domain_t`, and `blk%node_old` goes;
   - FUSS keeps passing its own arrays to `Morph_Update` one block at a time, as it does today.

   This is to be decided with the switch (plan 09 §11).

Smaller items noticed during the copy:

- The driver's registry is a copy of FUSS's `Registry.f90`, and ICE carries a registry of its own: a
  fix to one reaches none of the others. `check_sync.sh` reports FUSS-side changes.
- `generate_markdown` in the registry has no caller yet: FUSS documents its keys through
  `src/app/docgen.f90`, MORPH through this README.
- `amp`'s `>= 0` rule is FUSS's, and, as in FUSS, the registry does not check rules on vector keys.
