#!/usr/bin/env python3
"""Check one MORPH driver run on the mesh of generate.py.

    verify.py <workdir> [--gcl-tol X] [--pos-tol Y]

<workdir> is where run_case.sh ran MORPH: input.ini, run.log, INPUT/mesh.p3d,
OUTPUT/mesh.p3d and OUTPUT/gcl.dat. The law is read from input.ini. One PASS
or FAIL line per check; exit status 1 when any check fails. The two options
exist for the RED proof of the gate (records); a test never passes them.

THE GATE, every law: the per-step discrete-GCL residual MORPH logged, the max
over all cells of all blocks of |(V_new - V_old) - sum_f dV_swept(f)| / V, is at
most GCL_TOL in every step. Plan 09 records ~1e-15 for FUSS's moving cases; this
mesh gives 4.1e-16 at most (2026-10-24, gfortran 12 -O3 -march=native -flto, 1
and 8 threads). A GCL defect is not a round-off effect: a face swept with the
wrong sign or not at all leaves a residual of the order of the swept-outflow
ratio, 1e-3 and more here.

Then what each law must have done to the nodes, from the mesh MORPH read:
  static       nothing, bit for bit, and every logged residual exactly 0
  translation  x + vel*dt added once per step, bit for bit
  prescribed   x_ref + amp*sin(k.x_ref + omega*t_N + phase), within POS_TOL:
               the compiler contracts the argument into FMAs, Python does not
  nodeshift    the receding face moved by nodeshift-vel*dt per step and the
               opposite face pinned, both bit for bit; the lines between them
               compressed uniformly (a replica of the law's arithmetic, within
               POS_TOL for the same reason); every other block untouched
and, for every law, the nodes the blocks share still the same bits in each.

Python 3.6, standard library only.
"""
import configparser
import math
import os
import sys

sys.dont_write_bytecode = True   # no __pycache__ in the source tree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate  # noqa: E402  (the mesh: BLOCKS, INTERFACES, read_p3d)

GCL_TOL = 1.0e-14   # relative; measured max 4.1e-16 (above)
POS_TOL = 1.0e-15   # m, on a 0.4 m domain; measured below in the PASS lines

#: the prescribed law's phases, Morph_Motion_Prescribed.f90
PHASE = (0.0, 2.0943951023931953, 4.1887902047863905)

failures = []


def check(ok, what):
    print("  %s  %s" % ("PASS" if ok else "FAIL", what))
    if not ok:
        failures.append(what)
    return ok


def ulp(x):
    """Spacing of the doubles at |x| (math.ulp is Python 3.9+)."""
    if x == 0.0:
        return 5e-324
    return 2.0 ** (math.frexp(abs(x))[1] - 53)


def read_ini(path):
    cp = configparser.ConfigParser(inline_comment_prefixes=(";",), comment_prefixes=(";", "#"))
    cp.optionxform = str
    with open(path) as f:
        cp.read_file(f)
    return cp


def vec(cp, sec, key, default):
    raw = cp.get(sec, key, fallback=None)
    if raw is None:
        return list(default)
    return [float(t) for t in raw.split()]


def read_gcl_log(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            tok = line.split()
            rows.append((int(tok[0]), float(tok[1]), float(tok[2]), float(tok[3]), int(tok[4]), float(tok[5])))
    return rows


def max_diff(got, want):
    """Largest |got - want| over every node and coordinate, and its size in ulps
    of the coordinate where it is found."""
    worst, worst_ulp = 0.0, 0.0
    for b, (_, nodes) in enumerate(got):
        for key, x in nodes.items():
            for c in range(3):
                d = abs(x[c] - want[b][key][c])
                if d > worst:
                    worst, worst_ulp = d, d / ulp(want[b][key][c])
    return worst, worst_ulp


def same_bits(got, want, blocks=None, where=None):
    """Nodes of `got` equal to `want` bit for bit (in `blocks`, at the keys `where(b, key)` accepts)."""
    n, bad = 0, 0
    for b, (_, nodes) in enumerate(got):
        if blocks is not None and b not in blocks:
            continue
        for key, x in nodes.items():
            if where is not None and not where(b, key):
                continue
            n += 1
            if tuple(x) != tuple(want[b][key]):
                bad += 1
    return n, bad


def main():
    args = sys.argv[1:]
    global GCL_TOL, POS_TOL
    if "--gcl-tol" in args:
        i = args.index("--gcl-tol"); GCL_TOL = float(args[i + 1]); del args[i:i + 2]
    if "--pos-tol" in args:
        i = args.index("--pos-tol"); POS_TOL = float(args[i + 1]); del args[i:i + 2]
    if len(args) != 1:
        print(__doc__)
        return 2
    wd = args[0]

    cp = read_ini(os.path.join(wd, "input.ini"))
    law = cp.get("MORPH-Motion", "law", fallback="static").strip()
    steps = int(cp.get("MORPH-Time", "steps"))
    dt = float(cp.get("MORPH-Time", "dt"))
    t0 = float(cp.get("MORPH-Time", "t0", fallback="0.0"))
    print("verify.py: law %s, %d steps of %g s, in %s" % (law, steps, dt, wd))

    # ---- the run itself
    with open(os.path.join(wd, "run.log")) as f:
        log = f.read()
    check(" MORPH: done" in log and "[ERROR" not in log, "run.log ends with ' MORPH: done' and holds no [ERROR")

    # ---- the GCL log: one row per step, the time accumulated as MORPH does
    rows = read_gcl_log(os.path.join(wd, "OUTPUT", "gcl.dat"))
    check([r[0] for r in rows] == list(range(1, steps + 1)), "gcl.dat has one row per step, 1..%d" % steps)
    t, times_ok = t0, True
    for r in rows:
        t = t + dt
        times_ok = times_ok and (r[1] == t)
    check(times_ok, "gcl.dat times are t0 + dt summed step by step, to the last bit")

    worst = max([r[3] for r in rows] or [0.0])
    worst_step = max(rows, key=lambda r: r[3])[0] if rows else 0
    check(worst <= GCL_TOL,
          "GATE per-step discrete-GCL residual <= %.1e: max %.3e (step %d)" % (GCL_TOL, worst, worst_step))
    print("        max absolute residual %.3e m^3" % max([r[2] for r in rows] or [0.0]))

    # ---- the meshes
    mesh_in = generate.read_p3d(os.path.join(wd, "INPUT", "mesh.p3d"))
    mesh_out = generate.read_p3d(os.path.join(wd, "OUTPUT", "mesh.p3d"))
    if not check([m[0] for m in mesh_out] == [m[0] for m in mesh_in] == [generate.dims(b) for b in range(len(generate.BLOCKS))],
                 "mesh-out has the blocks and the dimensions of mesh-in (%d blocks)" % len(mesh_in)):
        return finish()
    x_in = [m[1] for m in mesh_in]

    if law == "static":
        n, bad = same_bits(mesh_out, x_in)
        check(bad == 0, "static: all %d nodes are the bits read (%d differ)" % (n, bad))
        check(all(r[2] == 0.0 and r[3] == 0.0 and r[5] == 0.0 for r in rows),
              "static: every logged residual and swept-outflow ratio is exactly 0")

    elif law == "translation":
        vel = vec(cp, "MORPH-Motion", "vel", [0.0, 0.0, 0.0])
        d = [v * dt for v in vel]
        want = []
        for nodes in x_in:
            moved = {}
            for key, x in nodes.items():
                y = list(x)
                for _ in range(steps):
                    y = [y[c] + d[c] for c in range(3)]
                moved[key] = tuple(y)
            want.append(moved)
        n, bad = same_bits(mesh_out, want)
        check(bad == 0, "translation: all %d nodes moved by vel*dt per step, bit for bit (%d differ)" % (n, bad))

    elif law == "prescribed":
        amp = vec(cp, "MORPH-Motion", "amp", [0.0] * 3)
        kx = vec(cp, "MORPH-Motion", "kx", [0.0] * 3)
        ky = vec(cp, "MORPH-Motion", "ky", [0.0] * 3)
        kz = vec(cp, "MORPH-Motion", "kz", [0.0] * 3)
        omega = float(cp.get("MORPH-Motion", "omega", fallback="0.0"))
        taper = cp.get("MORPH-Motion", "taper-to-boundary", fallback=".false.").strip().lower()
        if taper not in (".false.", "false", "f"):
            check(False, "prescribed: verify.py does not model taper-to-boundary")
            return finish()
        tN = t0
        for _ in range(steps):
            tN = tN + dt
        want = []
        for nodes in x_in:
            moved = {}
            for key, xr in nodes.items():
                moved[key] = tuple(
                    xr[c] + amp[c] * math.sin(kx[c] * xr[0] + ky[c] * xr[1] + kz[c] * xr[2] + omega * tN + PHASE[c])
                    for c in range(3))
            want.append(moved)
        dmax, umax = max_diff(mesh_out, want)
        check(dmax <= POS_TOL, "prescribed: nodes within %.1e m of x_ref + amp*sin(k.x_ref + omega*t_N + phase): "
                               "max %.2e m (%.1f ulp)" % (POS_TOL, dmax, umax))

    elif law == "nodeshift":
        b = int(cp.get("MORPH-Motion", "nodeshift-block")) - 1
        f = int(cp.get("MORPH-Motion", "nodeshift-face"))
        dirn = (f + 1) // 2 - 1                  # 0-based grid direction of the receding lines
        side = 2 - f % 2                         # 1: the low-index face recedes, 2: the high one
        d = [v * dt for v in vec(cp, "MORPH-Motion", "nodeshift-vel", [0.0] * 3)]
        dims = generate.dims(b)
        n = dims[dirn]
        moving = 0 if side == 1 else n
        pinned = n if side == 1 else 0

        def line_key(m, key):
            """the node m steps from the pinned end on the line through `key`"""
            k = list(key)
            k[dirn] = (n - m) if side == 1 else m
            return tuple(k)

        # replica of Morph_Motion_Nodeshift, step by step
        x = dict(x_in[b])
        line_starts = [key for key in x if key[dirn] == pinned]
        for _ in range(steps):
            old = dict(x)
            dd = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
            for key0 in line_starts:
                arc = [0.0]
                for m in range(1, n + 1):
                    p, q = old[line_key(m, key0)], old[line_key(m - 1, key0)]
                    seg = [p[c] - q[c] for c in range(3)]
                    arc.append(arc[-1] + math.sqrt(seg[0] * seg[0] + seg[1] * seg[1] + seg[2] * seg[2]))
                L = arc[n]
                if L - dd <= 0.0:
                    check(False, "nodeshift: verify.py's replica consumed a line; the case is wrong")
                    return finish()
                for m in range(1, n + 1):
                    key = line_key(m, key0)
                    x[key] = tuple(old[key][c] + d[c] * (arc[m] / L) for c in range(3))
        want = [dict(nodes) for nodes in x_in]
        want[b] = x
        # the surface: x_in + d once per step, exactly (s/L = 1 bitwise in the law)
        surf = [dict(nodes) for nodes in x_in]
        for key in surf[b]:
            if key[dirn] == moving:
                y = list(surf[b][key])
                for _ in range(steps):
                    y = [y[c] + d[c] for c in range(3)]
                surf[b][key] = tuple(y)
        nm, bad = same_bits(mesh_out, surf, blocks=[b], where=lambda bb, key: key[dirn] == moving)
        check(bad == 0, "nodeshift: the %d nodes of the receding face moved by nodeshift-vel*dt per step, "
                        "bit for bit (%d differ)" % (nm, bad))
        npin, bad = same_bits(mesh_out, x_in, blocks=[b], where=lambda bb, key: key[dirn] == pinned)
        check(bad == 0, "nodeshift: the %d nodes of the opposite face did not move (%d differ)" % (npin, bad))
        others = [bb for bb in range(len(x_in)) if bb != b]
        no, bad = same_bits(mesh_out, x_in, blocks=others)
        check(bad == 0, "nodeshift: the %d nodes of the other blocks did not move (%d differ)" % (no, bad))
        dmax, umax = max_diff([mesh_out[b]], [want[b]])
        check(dmax <= POS_TOL, "nodeshift: lines compressed uniformly, within %.1e m of the replica: "
                               "max %.2e m (%.1f ulp)" % (POS_TOL, dmax, umax))

    else:
        check(False, "verify.py knows no law '%s'" % law)
        return finish()

    # ---- every law: the blocks still share their interface nodes, bit for bit
    out = [m[1] for m in mesh_out]
    ok, nshared = True, 0
    for a, c, (I0, I1), (J0, J1), (K0, K1) in generate.INTERFACES:
        oa = [generate.BLOCKS[a][1][0], generate.BLOCKS[a][2][0], generate.BLOCKS[a][3][0]]
        oc = [generate.BLOCKS[c][1][0], generate.BLOCKS[c][2][0], generate.BLOCKS[c][3][0]]
        for K in range(K0, K1 + 1):
            for J in range(J0, J1 + 1):
                for I in range(I0, I1 + 1):
                    nshared += 1
                    pa = out[a][(I - oa[0], J - oa[1], K - oa[2])]
                    pc = out[c][(I - oc[0], J - oc[1], K - oc[2])]
                    ok = ok and (pa == pc)
    check(ok, "the %d interface nodes (3 interfaces) are the same bits in the blocks that share them" % nshared)

    return finish()


def finish():
    if failures:
        print("verify.py: %d check(s) FAILED" % len(failures))
        return 1
    print("verify.py: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
