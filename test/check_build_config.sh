#!/bin/bash
#===============================================================================
#         FILE: test/check_build_config.sh
#        USAGE: ./test/check_build_config.sh [--quiet] [--preset NAME] [--tree DIR]
#  DESCRIPTION: Assert that a build tree was configured the way CMakePresets.json
#               says (install.sh build writes the preset from the tree it built).
#
#               FUSS's test/check_build_config.sh at 2943354, with MORPH's
#               options. Its reason, from FUSS: build/ was once rebuilt with a
#               hand-written cmake line instead of the preset, CMake's defaults
#               were the opposite of the preset on two options that change
#               results, and the binary configured, compiled and ran without a
#               single complaint. MORPH's CMake defaults are USE_OPENMP OFF and
#               USE_TECIO OFF; install.sh's preset also names the ORION and the
#               FiNeR trees the build was made against.
#
#               --tree DIR (ctest's guard_build-config): check DIR only when it
#               is the preset's binaryDir; otherwise, or with no preset, exit 77
#               (ctest: skipped). Without --tree the preset's own tree is
#               checked and a missing preset is not an error (exit 0).
#===============================================================================
set -u

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
QUIET=0
NAME=default
TREE=""
while test $# -gt 0; do
  case $1 in
    --quiet)  QUIET=1; shift ;;
    --preset) NAME=$2; shift 2 ;;
    --tree)   TREE=$(cd "$2" 2>/dev/null && pwd -P) || { echo "  build-config: no directory $2" >&2; exit 1; }; shift 2 ;;
    *) echo "usage: $0 [--quiet] [--preset NAME] [--tree DIR]" >&2; exit 1 ;;
  esac
done

PRESET="$ROOT/CMakePresets.json"

if [ ! -f "$PRESET" ]; then
  echo "  build-config: no CMakePresets.json, nothing to check against" >&2
  [ -n "$TREE" ] && exit 77
  exit 0
fi

python3 - "$PRESET" "$ROOT" "$NAME" "$QUIET" "$TREE" <<'PY'
import json, os, re, sys

preset_path, root, name, quiet, tree = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == "1", sys.argv[5]

with open(preset_path) as fh:
    doc = json.load(fh)
presets = {p.get("name"): p for p in doc.get("configurePresets", [])}

def resolved(n, seen=()):
    p = presets.get(n)
    if p is None or n in seen:
        return {}, None
    cv, bdir = {}, None
    inh = p.get("inherits", [])
    for parent in ([inh] if isinstance(inh, str) else inh):
        pcv, pdir = resolved(parent, seen + (n,))
        cv.update(pcv)
        bdir = pdir or bdir
    cv.update(p.get("cacheVariables", {}))
    return cv, p.get("binaryDir", bdir)

if name not in presets:
    print("  build-config: no '%s' preset, nothing to check" % name)
    sys.exit(77 if tree else 0)
want, bdir = resolved(name)
bdir = (bdir or "${sourceDir}/build").replace("${sourceDir}", root)
if tree and os.path.realpath(tree) != os.path.realpath(bdir):
    print("  build-config: %s is not the tree of preset %s (%s): skipped" % (tree, name, bdir))
    sys.exit(77)
cache_path = os.path.join(bdir, "CMakeCache.txt")
if not os.path.isfile(cache_path):
    print("  build-config: %s is not configured. Run: cmake --preset %s" % (os.path.relpath(bdir, root), name))
    sys.exit(1)

# Read the cache as NAME:TYPE=VALUE
have = {}
with open(cache_path) as fh:
    for line in fh:
        m = re.match(r"^([A-Za-z_0-9]+):[A-Z]+=(.*)$", line.strip())
        if m:
            have[m.group(1)] = m.group(2)

def norm(v):
    v = str(v).strip()
    low = v.lower()
    if low in ("on", "true", "1", "yes"):   return "true"
    if low in ("off", "false", "0", "no"):  return "false"
    return v.rstrip("/")          # paths in the preset may carry a trailing /

bad = []
for k, v in want.items():
    if isinstance(v, dict):       # {"type": ..., "value": ...} form
        v = v.get("value", "")
    if k not in have:
        bad.append((k, norm(v), "<absent from cache>"))
    elif norm(have[k]) != norm(v):
        bad.append((k, norm(v), norm(have[k])))

tree = os.path.relpath(bdir, root)
if bad:
    print("  build-config: %s/ DOES NOT match CMakePresets.json (preset %s)" % (tree, name))
    for k, w, h in bad:
        print("      %-28s preset=%-12s cache=%s" % (k, w, h))
    print("  Results from this build are not those of the preset.")
    print("  Fix with:  rm -rf %s && cmake --preset %s && cmake --build %s -j 8" % (tree, name, tree))
    sys.exit(1)

if not quiet:
    print("  build-config: %s/ matches preset %s of CMakePresets.json (%d options checked)" % (tree, name, len(want)))
sys.exit(0)
PY
