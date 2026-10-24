#!/bin/bash
#===============================================================================
#         FILE: test/same_files.sh
#        USAGE: same_files.sh <a1> <b1> [<a2> <b2> ...]
#  DESCRIPTION: Each pair of files must be byte-identical: the ctest tests
#               that compare a 1-thread run with an 8-thread run. A missing
#               file is a difference.
#===============================================================================
set -u
if [ $# -lt 2 ] || [ $(( $# % 2 )) -ne 0 ]; then sed -n '3,3p' "$0" >&2; exit 2; fi
fail=0
while [ $# -gt 0 ]; do
  if cmp -s "$1" "$2"; then
    echo "  PASS  same bytes: $1 $2"
  else
    echo "  FAIL  $1 and $2 differ:"
    diff "$1" "$2" 2>&1 | head -10 | sed 's/^/          /'
    fail=1
  fi
  shift 2
done
exit $fail
