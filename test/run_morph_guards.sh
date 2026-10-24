#!/bin/bash
#===============================================================================
#         FILE: test/run_morph_guards.sh
#        USAGE: ./test/run_morph_guards.sh [--build DIR] [GUARD ...]
#  DESCRIPTION: Architecture guards for the MORPH mesh component (plan 05 5.0,
#               plan 09 section 2).
#
#               The entire value of MORPH being a separate component rests on
#               one rule -- it must never depend on FUSS -- and on it staying
#               true as later phases add code. A design rule nobody checks
#               decays within two phases, so it is checked here, mechanically,
#               in milliseconds.
#
#               FUSS's test/run_morph_guards.sh at 2943354, adapted to the
#               standalone layout (library src/lib, driver src/app, unit
#               programs test/unit) and split into named guards, each one a
#               ctest test (label guard, test/CMakeLists.txt):
#
#                 no-fuss         no MORPH source uses a FUSS module   (FUSS a)
#                 no-stop         the library never stops: no stop or error
#                                 stop statement under src/lib         (FUSS c)
#                 no-io           the library does no I/O: no print, no
#                                 write to a unit, no open/close/inquire  (new)
#                 lib-isolated    the library uses no module but
#                                 iso_fortran_env, omp_lib and its own: no
#                                 ORION, no FiNeR (the driver's)          (new)
#                 wall-clean      the library compiles standalone, with no
#                                 other source on the path, and -Wall prints
#                                 nothing at all                       (FUSS d)
#                 wall-clean-app  the driver compiles against the modules of
#                                 the build tree DIR, -Wall prints nothing
#                                 but gfortran's hidden-length false
#                                 positive (below)                        (new)
#                 units           the unit programs built here with the
#                                 library sources, at 1 and 8 threads, same
#                                 output                               (FUSS e)
#
#               No GUARD runs them all; wall-clean-app is then skipped unless
#               --build names a configured and built tree. ctest builds the unit
#               programs itself, so it does not register `units`. FUSS's guard
#               (b), exactly one FUSS file naming MORPH types, is about FUSS's
#               tree and stays in FUSS.
#
#               FC selects the compiler (default gfortran; the flags are GNU's).
#               Exit status 1 when a guard fails.
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LIB="$ROOT/src/lib"
APP="$ROOT/src/app"
FC=${FC:-gfortran}
BUILD=""
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

fail=0
pass() { printf '  \033[0;32mPASS\033[0m  %s\n' "$1"; }
bad()  { printf '  \033[0;31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }
skip() { printf '  SKIP  %s\n' "$1"; }

GUARDS=()
while [ $# -gt 0 ]; do
  case $1 in
    --build) BUILD=$(cd "$2" 2>/dev/null && pwd) || { echo "no build tree $2" >&2; exit 2; }; shift 2 ;;
    no-fuss|no-stop|no-io|lib-isolated|wall-clean|wall-clean-app|units) GUARDS+=("$1"); shift ;;
    *) sed -n '3,4p' "$0" >&2; echo "unknown guard $1" >&2; exit 2 ;;
  esac
done
ALL=0
[ ${#GUARDS[@]} -eq 0 ] && ALL=1 && GUARDS=(no-fuss no-stop no-io lib-isolated wall-clean wall-clean-app units)

# The library in dependency order: Morph_Metrics uses Morph_GCL for the cell
# volume, so GCL must be compiled first or its .mod does not exist yet; the
# domain type uses the API. A new library file must be added here: the guards
# that compile refuse a list that does not cover src/lib.
LIB_SRC="
base/Morph_Types_m.f90
metrics/Morph_GCL.f90
metrics/Morph_Metrics.f90
motion/Morph_Motion_m.f90
motion/Morph_Motion_Static.f90
motion/Morph_Motion_Prescribed.f90
motion/Morph_Motion_Translation.f90
motion/Morph_Motion_Nodeshift.f90
quality/Morph_Quality.f90
Morph_API.f90
domain/Morph_Domain.f90
"
APP_SRC="
config/Morph_Config_Types.f90
config/Morph_Registry.f90
config/Register_Morph.f90
config/Morph_Backend_INI.f90
config/Morph_Read_Ini.f90
io/Morph_IO.f90
main.f90
"

# The Fortran sources of a tree, relative to it, sorted
sources() { (cd "$1" && find . -name '*.f90' -o -name '*.F90' | sed 's|^\./||' | sort); }

# 0 when LIST (relative to DIR) is exactly the Fortran sources of DIR
covers() {
  local dir=$1 list=$2 what=$3 diffs
  diffs=$(diff <(printf '%s\n' $list | sort) <(sources "$dir"))
  if [ -n "$diffs" ]; then
    bad "$what: the source list of run_morph_guards.sh is not the content of ${dir#$ROOT/}:"
    echo "$diffs" | grep '^[<>]' | sed -e 's/^</          listed, not found:/' -e 's/^>/          found, not listed:/'
    return 1
  fi
  return 0
}

# Module names used by Fortran files: "file:line: module", comments stripped
uses_of() {
  awk '{ l = tolower($0); sub(/!.*/, "", l) }
       l ~ /^[ \t]*use([ \t]|,|::)/ {
         m = l
         sub(/^[ \t]*use[ \t]*(,[ \t]*(non_)?intrinsic[ \t]*)?(::)?[ \t]*/, "", m)
         sub(/[^a-z0-9_].*$/, "", m)
         print FILENAME ":" FNR ": " m
       }' "$@"
}

echo "MORPH architecture guards"
echo "  library: $LIB"
echo

for guard in "${GUARDS[@]}"; do
  case $guard in

  # ---------------------------------------------------------------------------
  # (a) MORPH must not use any FUSS module. This is THE rule.
  # ---------------------------------------------------------------------------
  no-fuss)
    hits=$(uses_of $(find "$ROOT/src" "$ROOT/test" -name '*.f90' -o -name '*.F90' | sort) | grep -E ': fuss_' || true)
    if [ -z "$hits" ]; then
      pass "MORPH uses no FUSS module (src/, test/)"
    else
      bad "MORPH depends on FUSS -- it is no longer a standalone component:"
      echo "$hits" | sed "s|$ROOT/||; s/^/          /"
    fi
    ;;

  # ---------------------------------------------------------------------------
  # (c) MORPH must never abort. A library that stops cannot be embedded.
  #     The driver (src/app) and the test programs are executables: they may.
  #     Comments are stripped first; a `stop` inside a string still counts.
  # ---------------------------------------------------------------------------
  no-stop)
    stops=$(awk '{ l = tolower($0); sub(/!.*/, "", l) }
                 l ~ /(^|[^a-z0-9_])stop([^a-z0-9_]|$)/ { print FILENAME ":" FNR ": " $0 }' \
                 $(find "$LIB" -name '*.f90' -o -name '*.F90' | sort))
    if [ -z "$stops" ]; then
      pass "the library contains no 'stop' or 'error stop'"
    else
      bad "the library aborts the program instead of returning a status:"
      echo "$stops" | sed "s|$ROOT/||; s/^/          /"
    fi
    ;;

  # ---------------------------------------------------------------------------
  # Nor may it print or touch a file: what a solver embedding it writes, and
  # where, is the solver's (plan 09 section 9: "No I/O"). An internal write
  # into a character variable is not I/O and is allowed.
  # ---------------------------------------------------------------------------
  no-io)
    io=$(awk '{ l = tolower($0); sub(/!.*/, "", l) }
              l ~ /(^|[^a-z0-9_])print([^a-z0-9_]|$)/ ||
              l ~ /(^|[^a-z0-9_])(write|read)[ \t]*\([ \t]*(unit[ \t]*=[ \t]*)?(\*|[0-9]+|output_unit|error_unit|input_unit)/ ||
              l ~ /(^|[^a-z0-9_])(open|close|inquire|flush|rewind|backspace|endfile)[ \t]*\(/ {
                print FILENAME ":" FNR ": " $0 }' \
              $(find "$LIB" -name '*.f90' -o -name '*.F90' | sort))
    if [ -z "$io" ]; then
      pass "the library does no I/O (print, write/read on a unit, open, close, inquire, flush)"
    else
      bad "the library does I/O, which belongs to its caller:"
      echo "$io" | sed "s|$ROOT/||; s/^/          /"
    fi
    ;;

  # ---------------------------------------------------------------------------
  # The library's only dependencies are the compiler's: what a solver must
  # supply to embed it is nothing. ORION and FiNeR belong to the driver.
  # ---------------------------------------------------------------------------
  lib-isolated)
    foreign=$(uses_of $(find "$LIB" -name '*.f90' -o -name '*.F90' | sort) \
              | grep -vE ': (iso_fortran_env|omp_lib|morph_[a-z0-9_]*)$' || true)
    if [ -z "$foreign" ]; then
      pass "the library uses only iso_fortran_env, omp_lib and Morph_* modules"
    else
      bad "the library uses a module it does not own:"
      echo "$foreign" | sed "s|$ROOT/||; s/^/          /"
    fi
    ;;

  # ---------------------------------------------------------------------------
  # (d) MORPH compiles standalone, with no other source on the path, and is
  #     -Wall clean: any compiler output is a failure.
  # ---------------------------------------------------------------------------
  wall-clean)
    covers "$LIB" "$LIB_SRC" wall-clean || continue
    mkdir -p "$WORK/lib"
    if ( cd "$WORK/lib" && $FC -c -O2 -fopenmp -fimplicit-none -Wall \
           -Wno-unused-dummy-argument $(printf "$LIB/%s " $LIB_SRC) ) > "$WORK/lib.log" 2>&1 \
       && [ ! -s "$WORK/lib.log" ]; then
      pass "the library compiles standalone, -Wall clean ($FC -O2 -fopenmp)"
    else
      bad "the library does not compile standalone without a word:"
      sed 's/^/          /' "$WORK/lib.log" | head -30
    fi
    ;;

  # ---------------------------------------------------------------------------
  # The driver against the modules of a build tree (ORION and FiNeR are built
  # there). One warning is gfortran's, not the code's: any deferred-length
  # character array passed to an intent(out) allocatable dummy (FiNeR's
  # get_sections_list, loop) gives "'.<name>' is used uninitialized" for the
  # hidden length the compiler made (gfortran 12.2.1; a reproducer of 19 lines
  # in two files, with no FiNeR, does it). A Fortran name cannot start with
  # '.', so only those lines are dropped; every other line fails the guard.
  # ---------------------------------------------------------------------------
  wall-clean-app)
    if [ -z "$BUILD" ]; then
      [ $ALL -eq 1 ] && { skip "wall-clean-app: no --build tree given"; continue; }
      bad "wall-clean-app needs --build <configured and built tree>"; continue
    fi
    covers "$APP" "$APP_SRC" wall-clean-app || continue
    incs=$(find "$BUILD" -name '*.mod' -printf '%h\n' | sort -u | sed 's/^/-I/')
    if [ -z "$incs" ]; then
      bad "wall-clean-app: no .mod file under $BUILD (not built?)"; continue
    fi
    mkdir -p "$WORK/app"
    ( cd "$WORK/app" && $FC -c -O2 -fopenmp -fimplicit-none -cpp -ffree-line-length-none -Wall \
        -Wno-unused-dummy-argument -J "$WORK/app" $incs $(printf "$APP/%s " $APP_SRC) ) > "$WORK/app.log" 2>&1
    rc=$?
    # gfortran prints a warning as: location line, blank, source line, caret
    # line, "Warning: ..." [, then a note block]. Keep the "Warning:" and
    # "Error:" lines that are not the hidden-length case.
    left=$(grep -E '^(Warning|Error|Fatal Error):' "$WORK/app.log" \
           | grep -vE "^Warning: .\.[A-Za-z0-9_]+. is used uninitialized \[-Wuninitialized\]$" || true)
    if [ $rc -eq 0 ] && [ -z "$left" ]; then
      n=$(grep -cE "^Warning: .\.[A-Za-z0-9_]+. is used uninitialized" "$WORK/app.log")
      pass "the driver compiles -Wall clean against $BUILD ($n hidden-length false positive(s) dropped)"
    else
      bad "the driver does not compile -Wall clean (exit status $rc):"
      sed 's/^/          /' "$WORK/app.log" | head -40
    fi
    ;;

  # ---------------------------------------------------------------------------
  # (e) Unit tests: the GCL identity, single cell and whole block, and the
  #     nodeshift law; 1 and 8 threads must print the same.
  # ---------------------------------------------------------------------------
  units)
    covers "$LIB" "$LIB_SRC" units || continue
    mkdir -p "$WORK/units"
    for prog in test_morph_gcl test_morph_block test_morph_nodeshift; do
      if ( cd "$WORK/units" && $FC -O2 -fopenmp -fimplicit-none -o "$prog" \
             $(printf "$LIB/%s " $LIB_SRC) "$ROOT/test/unit/$prog.f90" ) >> "$WORK/units.log" 2>&1; then
        for nt in 1 8; do
          if ( cd "$WORK/units" && OMP_NUM_THREADS=$nt "./$prog" ) > "$WORK/units/$prog.$nt.out" 2>&1; then
            pass "$prog (OMP_NUM_THREADS=$nt)"
          else
            bad "$prog failed with OMP_NUM_THREADS=$nt:"
            sed 's/^/          /' "$WORK/units/$prog.$nt.out" | grep -vE '^\s+#|Error termination' | head -15
          fi
        done
        if ! diff -q "$WORK/units/$prog.1.out" "$WORK/units/$prog.8.out" >/dev/null 2>&1; then
          bad "$prog output differs between 1 and 8 threads (race in the metric refresh?)"
        fi
      else
        bad "could not build $prog:"
        sed 's/^/          /' "$WORK/units.log" | tail -15
      fi
    done
    ;;
  esac
done

echo
if [ "$fail" -eq 0 ]; then
  echo "all MORPH guards passed"
else
  echo "$fail MORPH guard(s) FAILED"
  exit 1
fi
