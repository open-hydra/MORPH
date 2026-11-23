#!/bin/bash
#===============================================================================
#         FILE: test/driver/run_case.sh
#        USAGE: run_case.sh <case-dir> <ini> <workdir> <threads> <MORPH> [<exit> <pattern>...]
#  DESCRIPTION: One ctest test of the MORPH driver (test/CMakeLists.txt).
#
#               <workdir> is wiped and filled: generate.py of <case-dir> writes
#               INPUT/mesh.p3d, <case-dir>/<ini> becomes input.ini, MORPH runs
#               there with OMP_NUM_THREADS=<threads> (its output in run.log),
#               and verify.py of <case-dir> checks what it wrote.
#
#               With <exit> and one or more <pattern>s, MORPH must end with
#               exit status <exit> and every <pattern> (an extended regex) must
#               match a line of run.log. A non-zero <exit> is a run that must
#               FAIL (MORPH ends with `error stop <exit>`), and verify.py is
#               then not run; with <exit> 0 the run must succeed, and verify.py
#               checks it as above.
#
#               PYTHON selects the interpreter (default python3).
#===============================================================================
set -u
if [ $# -ne 5 ] && [ $# -lt 7 ]; then
  sed -n '3,4p' "$0" >&2; exit 2
fi
CASE=$(cd "$1" && pwd) || exit 2
INI=$2; WORK=$3; NT=$4; EXE=$5
EXPECT=${6:-0}
shift 5
[ $# -gt 0 ] && shift           # what is left: the patterns
PY=${PYTHON:-python3}

rm -rf "$WORK" && mkdir -p "$WORK" && cd "$WORK" || exit 2
"$PY" "$CASE/generate.py" INPUT/mesh.p3d || { echo "run_case: FAIL generate.py"; exit 1; }
cp "$CASE/$INI" input.ini || exit 2

OMP_NUM_THREADS=$NT "$EXE" > run.log 2>&1
rc=$?
cat run.log
echo "run_case: MORPH exit status $rc ($INI, OMP_NUM_THREADS=$NT, in $WORK)"

# An OpenMP build reports its team size: the run must have had the threads asked for
if grep -q 'parallel execution' run.log; then
  got=$(sed -n 's/.*Number of threads --> *//p' run.log)
  if [ "$got" != "$NT" ]; then echo "run_case: FAIL ran on '$got' threads, not $NT"; exit 1; fi
fi

if [ $# -gt 0 ]; then
  if [ "$rc" -ne "$EXPECT" ]; then echo "run_case: FAIL exit status $rc, expected $EXPECT"; exit 1; fi
  for pattern in "$@"; do
    if ! grep -qE -- "$pattern" run.log; then echo "run_case: FAIL no line of run.log matches: $pattern"; exit 1; fi
  done
  if [ "$EXPECT" -ne 0 ]; then echo "run_case: refused as expected (exit status $EXPECT)"; exit 0; fi
  echo "run_case: run.log shows the $# line(s) expected"
fi

[ "$rc" -eq 0 ] || { echo "run_case: FAIL MORPH exit status $rc"; exit 1; }
"$PY" "$CASE/verify.py" "$WORK"
