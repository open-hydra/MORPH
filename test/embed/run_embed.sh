#!/bin/bash
#===============================================================================
#         FILE: test/embed/run_embed.sh
#        USAGE: run_embed.sh <MORPH tree> <work dir> <Fortran compiler>
#  DESCRIPTION: MORPH embedded in another CMake project, as FUSS will embed it
#               (README, next steps): the parent of test/embed adds the MORPH
#               tree with add_subdirectory, links MORPHL and runs probe.f90.
#               Twice: once with a module directory of the parent's (FUSS sets
#               one) and once without. Only MORPHL may be configured: no ORION,
#               no FiNeR, no driver, no tests.
#===============================================================================
set -u
[ $# -eq 3 ] || { sed -n '3,3p' "$0" >&2; exit 2; }
SRC=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MORPH=$1; WORK=$2; FC=$3
mkdir -p "$WORK" || exit 2
fail=0
for mod in ON OFF; do
  b="$WORK/module-dir-$mod"
  rm -rf "$b"
  if ! cmake -S "$SRC" -B "$b" -DMORPH_DIR="$MORPH" -DWITH_MODULE_DIR=$mod -DCMAKE_Fortran_COMPILER="$FC" \
       > "$b.configure.log" 2>&1; then
    echo "  FAIL  configure, parent module dir $mod:"; tail -20 "$b.configure.log"; fail=1; continue
  fi
  if grep -qE 'ORION|FiNeR|Python3' "$b.configure.log" || [ -d "$b/morph/src/app" ] || [ -d "$b/morph/test" ]; then
    echo "  FAIL  parent module dir $mod: more than MORPHL was configured"; fail=1
  fi
  if ! cmake --build "$b" -j 4 > "$b.build.log" 2>&1; then
    echo "  FAIL  build, parent module dir $mod:"; grep -E 'Error' -A3 "$b.build.log" | head -20; fail=1; continue
  fi
  if out=$("$b/probe" 2>&1); then
    echo "  PASS  parent module dir $mod: $out"
  else
    echo "  FAIL  probe, parent module dir $mod:"; echo "$out"; fail=1
  fi
done
exit $fail
