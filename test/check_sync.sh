#!/bin/bash
#===============================================================================
#         FILE: test/check_sync.sh
#        USAGE: ./test/check_sync.sh [--fuss DIR] [--diff] [REF]
#  DESCRIPTION: Drift guard between MORPH and the FUSS sources it was copied
#               from (README.md, "Where MORPH comes from").
#
#               MORPH is a copy of FUSS's src/lib/morph at FUSS 2943354 (user
#               decision 2026-10-24), and FUSS keeps building its own copy until
#               it switches to MORPH. The two can drift apart; this reports how:
#
#                 1. the files copied verbatim that MORPH has changed since
#                    (MORPH's working tree vs FUSS 2943354);
#                 2. the FUSS files MORPH copied or adapted that FUSS has
#                    changed since (FUSS 2943354 vs REF);
#                 3. the files under FUSS's src/lib/morph at REF that MORPH
#                    does not have (a new law, a new kernel).
#
#               An adapted file (the build, the guard script, the registry
#               copy) differs from its origin by design: only (2) is reported
#               for it, so that a FUSS change worth porting is seen.
#
#               FUSS is read through git only (cat-file, show, diff of two
#               commits, ls-tree): neither repository is changed, and FUSS's
#               working tree is not looked at, so commit there first. To
#               compare with what is pushed, fetch FUSS and pass origin/<branch>.
#
#               --fuss DIR  the FUSS repository (default: ../FUSS beside MORPH,
#                           hydra/src/lib/FUSS)
#               --diff      print the differences, not only their list
#               REF         a commit-ish of FUSS (default: the branch the copy
#                           was taken from, worktree-porous-pyrolysis)
#
#               Exit status 0: nothing differs; 1: something does; 2: FUSS, the
#               recorded commit or REF not found. Not a ctest test: FUSS moving
#               is not a MORPH defect.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RECORDED=2943354329820c795d715323a42eace656c58924
FUSS="$ROOT/../FUSS"
REF=worktree-porous-pyrolysis
SHOWDIFF=0
while [ $# -gt 0 ]; do
  case $1 in
    --fuss) FUSS=$2; shift 2 ;;
    --diff) SHOWDIFF=1; shift ;;
    -*)     sed -n '3,4p' "$0" >&2; exit 2 ;;
    *)      REF=$1; shift ;;
  esac
done

fgit() { git -C "$FUSS" "$@"; }

if ! fgit rev-parse --git-dir > /dev/null 2>&1; then
  echo "check_sync: no git repository at $FUSS (pass --fuss <FUSS repository>)" >&2; exit 2
fi
FUSS=$(cd "$FUSS" && pwd)
if ! fgit cat-file -e "$RECORDED^{commit}" 2>/dev/null; then
  echo "check_sync: FUSS at $FUSS does not have the recorded commit $RECORDED (fetch it)" >&2; exit 2
fi
if ! REFHASH=$(fgit rev-parse --verify --quiet "$REF^{commit}"); then
  echo "check_sync: FUSS has no commit '$REF'" >&2; exit 2
fi

# MORPH path | FUSS path at the recorded commit
VERBATIM="
src/lib/Morph_API.f90|src/lib/morph/Morph_API.f90
src/lib/base/Morph_Types_m.f90|src/lib/morph/base/Morph_Types_m.f90
src/lib/metrics/Morph_GCL.f90|src/lib/morph/metrics/Morph_GCL.f90
src/lib/metrics/Morph_Metrics.f90|src/lib/morph/metrics/Morph_Metrics.f90
src/lib/motion/Morph_Motion_m.f90|src/lib/morph/motion/Morph_Motion_m.f90
src/lib/motion/Morph_Motion_Static.f90|src/lib/morph/motion/Morph_Motion_Static.f90
src/lib/motion/Morph_Motion_Prescribed.f90|src/lib/morph/motion/Morph_Motion_Prescribed.f90
src/lib/motion/Morph_Motion_Translation.f90|src/lib/morph/motion/Morph_Motion_Translation.f90
src/lib/motion/Morph_Motion_Nodeshift.f90|src/lib/morph/motion/Morph_Motion_Nodeshift.f90
src/lib/quality/Morph_Quality.f90|src/lib/morph/quality/Morph_Quality.f90
test/unit/test_morph_block.f90|src/lib/morph/test/test_morph_block.f90
test/unit/test_morph_gcl.f90|src/lib/morph/test/test_morph_gcl.f90
test/unit/test_morph_nodeshift.f90|src/lib/morph/test/test_morph_nodeshift.f90
cmake/SetFortranFlags.cmake|cmake/SetFortranFlags.cmake
cmake/SetCompileFlag.cmake|cmake/SetCompileFlag.cmake
cmake/SetParallelizationLibrary.cmake|cmake/SetParallelizationLibrary.cmake
"
ADAPTED="
test/run_morph_guards.sh|test/run_morph_guards.sh
test/check_build_config.sh|test/check_build_config.sh
install.sh|install.sh
CMakeLists.txt|CMakeLists.txt
src/app/config/Morph_Registry.f90|src/lib/config/Registry.f90
src/app/config/Morph_Backend_INI.f90|src/lib/config/Backend_INI.f90
src/app/config/Morph_Read_Ini.f90|src/lib/config/Read_Ini.f90
"

ndiff=0
fuss_side() {   # "same" when FUSS's file is the same at the recorded commit and at REF
  local fp=$1
  if ! fgit cat-file -e "$REFHASH:$fp" 2>/dev/null; then echo "GONE at REF"; return 1; fi
  if fgit diff --quiet "$RECORDED" "$REFHASH" -- "$fp"; then echo same; return 0; fi
  echo CHANGED; return 1
}

echo "MORPH sync check"
echo "  MORPH  $ROOT"
echo "  FUSS   $FUSS"
echo "  copied from FUSS $RECORDED"
echo "  compared with REF $REF = $REFHASH"
if ! fgit merge-base --is-ancestor "$RECORDED" "$REFHASH"; then
  echo "  NOTE: REF does not descend from the recorded commit: the FUSS-side column also shows what REF lacks"
fi
echo

printf '  %-48s %-22s %s\n' "copied verbatim" "MORPH vs FUSS@${RECORDED:0:7}" "FUSS ${RECORDED:0:7} vs REF"
for pair in $VERBATIM; do
  mp=${pair%%|*}; fp=${pair#*|}
  if [ ! -f "$ROOT/$mp" ]; then
    m="MISSING"; ndiff=$((ndiff+1))
  elif cmp -s "$ROOT/$mp" <(fgit show "$RECORDED:$fp"); then
    m="same"
  else
    m="CHANGED"; ndiff=$((ndiff+1))
  fi
  f=$(fuss_side "$fp") || ndiff=$((ndiff+1))
  printf '  %-48s %-22s %s\n' "$mp" "$m" "$f"
  if [ $SHOWDIFF -eq 1 ] && [ "$m" = CHANGED ]; then
    diff -u <(fgit show "$RECORDED:$fp") "$ROOT/$mp" --label "FUSS@${RECORDED:0:7}:$fp" --label "MORPH:$mp" | sed 's/^/      /'
  fi
  if [ $SHOWDIFF -eq 1 ] && [ "$f" = CHANGED ]; then
    fgit diff "$RECORDED" "$REFHASH" -- "$fp" | sed 's/^/      /'
  fi
done
echo

printf '  %-38s    %-30s %s\n' "adapted" "from FUSS" "FUSS ${RECORDED:0:7} vs REF"
for pair in $ADAPTED; do
  mp=${pair%%|*}; fp=${pair#*|}
  f=$(fuss_side "$fp") || ndiff=$((ndiff+1))
  printf '  %-38s <- %-30s %s\n' "$mp" "$fp" "$f"
  if [ $SHOWDIFF -eq 1 ] && [ "$f" = CHANGED ]; then
    fgit diff "$RECORDED" "$REFHASH" -- "$fp" | sed 's/^/      /'
  fi
done
echo

known=$(for pair in $VERBATIM; do echo "${pair#*|}"; done | sort)
new=$(comm -13 <(echo "$known") <(fgit ls-tree -r --name-only "$REFHASH" -- src/lib/morph | sort))
if [ -z "$new" ]; then
  echo "  in FUSS's src/lib/morph at REF, not in MORPH: none"
else
  echo "  in FUSS's src/lib/morph at REF, not in MORPH:"
  echo "$new" | sed 's/^/    /'
  ndiff=$((ndiff + $(echo "$new" | grep -c .)))
fi
echo

if [ $ndiff -eq 0 ]; then
  echo "check_sync: in sync with FUSS $REF"
  exit 0
fi
echo "check_sync: $ndiff difference(s); --diff prints them"
exit 1
