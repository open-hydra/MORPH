#!/bin/bash
#===============================================================================
#         FILE: test/unit/run_unit.sh
#        USAGE: run_unit.sh <program> <threads> <out>
#  DESCRIPTION: One ctest test of a unit program: run it with
#               OMP_NUM_THREADS=<threads>, keep what it prints in <out> (the
#               1- and 8-thread outputs are then compared, same_files.sh) and
#               show it. The exit status is the program's: the unit programs
#               end with `error stop 1` on any failed check.
#===============================================================================
set -u
[ $# -eq 3 ] && [ -n "$3" ] || { sed -n '3,3p' "$0" >&2; exit 2; }
mkdir -p "$(dirname "$3")" || exit 2
# Each run in a directory of its own, <out>.work: some programs write files
# where they run, and ctest may run the 1- and 8-thread runs of one program at
# the same time. <program> and <out> are absolute paths (test/CMakeLists.txt).
work="$3.work"
rm -rf "$work" && mkdir -p "$work" || exit 2
( cd "$work" && OMP_NUM_THREADS=$2 "$1" ) > "$3" 2>&1
rc=$?
cat "$3"
exit $rc
