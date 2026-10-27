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
  nodeshift    the receding face moved by its displacement per step
               (nodeshift-vel*dt, plus the bump of nodeshift-bump-vel when set:
               bit for bit without a bump, within POS_TOL with one) and the
               opposite face pinned, bit for bit; every other block untouched;
               then, by interior:
    line       the lines between them compressed uniformly, each on its own
               (a replica of the law's arithmetic, within POS_TOL)
    smooth     the four edges the replica of the line law; the solved nodes
               (the four side faces with their in-face springs, the interior
               with all six) in spring equilibrium in EVERY step, from the
               snapshots out-diter = 1 writes: the residual
               |sum_q w (d_q - d_p)| / sum_q w, w = 1/|e| on the step's old
               nodes, at most SMOOTH_RES_TOL times the largest face
               displacement; every node's displacement between the prescribed
               ones (the maximum principle); a coordinate no face node moves
               in untouched, bit for bit
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
SMOOTH_RES_TOL = 1.0e-11   # the equilibrium residual over the face displacement: 10 x smooth-tolerance

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


def max_diff(got, want, where=None):
    """Largest |got - want| over every node (at the keys `where(b, key)` accepts)
    and coordinate, and its size in ulps of the coordinate where it is found."""
    worst, worst_ulp = 0.0, 0.0
    for b, (_, nodes) in enumerate(got):
        for key, x in nodes.items():
            if where is not None and not where(b, key):
                continue
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
        vel = vec(cp, "MORPH-Motion", "nodeshift-vel", [0.0] * 3)
        bump_vel = vec(cp, "MORPH-Motion", "nodeshift-bump-vel", [0.0] * 3)
        width = float(cp.get("MORPH-Motion", "nodeshift-bump-width", fallback="0.25"))
        interior = cp.get("MORPH-Motion", "interior", fallback="line").strip()
        out_diter = int(cp.get("MORPH-IO", "out-diter", fallback="0"))
        bump = any(v != 0.0 for v in bump_vel)
        dims = generate.dims(b)
        n = dims[dirn]
        t1, t2 = [c for c in range(3) if c != dirn]   # the transverse grid directions
        moving = 0 if side == 1 else n
        pinned = n if side == 1 else 0
        print("        interior %s, %s" % (interior, "a bump of %g m/s, width %g" % (bump_vel[0], width) if bump
                                           else "a uniform displacement"))

        def line_key(m, key):
            """the node m steps from the pinned end on the line through `key`"""
            k = list(key)
            k[dirn] = (n - m) if side == 1 else m
            return tuple(k)

        def disp(key):
            """this step's displacement of the surface node of the line through
            `key`: Set_Shift's arithmetic (main.f90), the old expression without a bump"""
            if not bump:
                return [v * dt for v in vel]
            r2 = 0.0
            if dims[t1] > 1:
                xi = key[t1] / dims[t1] - 0.5
                r2 = r2 + xi * xi
            if dims[t2] > 1:
                eta = key[t2] / dims[t2] - 0.5
                r2 = r2 + eta * eta
            g = max(0.0, 1.0 - r2 / (width * width))
            g = g * g
            return [(vel[c] + bump_vel[c] * g) * dt for c in range(3)]

        def is_edge(key):
            return key[t1] in (0, dims[t1]) and key[t2] in (0, dims[t2])

        def line_replica(line_starts):
            """replica of the line law (Morph_Motion_Nodeshift), step by step, on
            the lines through `line_starts` (nodes of the pinned face); None when
            a line is consumed"""
            x = dict(x_in[b])
            for _ in range(steps):
                old = dict(x)
                for key0 in line_starts:
                    d = disp(key0)
                    dd = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                    if dd == 0.0:
                        continue
                    arc = [0.0]
                    for m in range(1, n + 1):
                        p, q = old[line_key(m, key0)], old[line_key(m - 1, key0)]
                        seg = [p[c] - q[c] for c in range(3)]
                        arc.append(arc[-1] + math.sqrt(seg[0] * seg[0] + seg[1] * seg[1] + seg[2] * seg[2]))
                    L = arc[n]
                    if L - dd <= 0.0:
                        return None
                    for m in range(1, n + 1):
                        key = line_key(m, key0)
                        x[key] = tuple(old[key][c] + d[c] * (arc[m] / L) for c in range(3))
            return x

        # the surface: x_in + d once per step (s/L = 1 bitwise in the law)
        surf = [dict(nodes) for nodes in x_in]
        for key in surf[b]:
            if key[dirn] == moving:
                d = disp(key)
                y = list(surf[b][key])
                for _ in range(steps):
                    y = [y[c] + d[c] for c in range(3)]
                surf[b][key] = tuple(y)
        on_face = lambda bb, key: key[dirn] == moving
        if bump:
            dmax, umax = max_diff([mesh_out[b]], [surf[b]], where=on_face)
            check(dmax <= POS_TOL, "nodeshift: the receding face moved by its displacement per step, within %.1e m: "
                                   "max %.2e m (%.1f ulp; the bump's arithmetic is contracted into FMAs)"
                                   % (POS_TOL, dmax, umax))
        else:
            nm, bad = same_bits(mesh_out, surf, blocks=[b], where=on_face)
            check(bad == 0, "nodeshift: the %d nodes of the receding face moved by nodeshift-vel*dt per step, "
                            "bit for bit (%d differ)" % (nm, bad))
        npin, bad = same_bits(mesh_out, x_in, blocks=[b], where=lambda bb, key: key[dirn] == pinned)
        check(bad == 0, "nodeshift: the %d nodes of the opposite face did not move (%d differ)" % (npin, bad))
        others = [bb for bb in range(len(x_in)) if bb != b]
        no, bad = same_bits(mesh_out, x_in, blocks=others)
        check(bad == 0, "nodeshift: the %d nodes of the other blocks did not move (%d differ)" % (no, bad))

        line_starts = [key for key in x_in[b] if key[dirn] == pinned]
        if interior == "line":
            x = line_replica(line_starts)
            if x is None:
                check(False, "nodeshift: verify.py's replica consumed a line; the case is wrong")
                return finish()
            dmax, umax = max_diff([mesh_out[b]], [x])
            check(dmax <= POS_TOL, "nodeshift line: every line compressed uniformly on its own, within %.1e m of "
                                   "the replica: max %.2e m (%.1f ulp)" % (POS_TOL, dmax, umax))

        elif interior == "smooth":
            smooth_tol = float(cp.get("MORPH-Motion", "smooth-tolerance", fallback="1.0e-12"))
            # the four edges: the line law
            x = line_replica([key for key in line_starts if is_edge(key)])
            if x is None:
                check(False, "nodeshift: verify.py's replica consumed a line; the case is wrong")
                return finish()
            dmax, umax = max_diff([mesh_out[b]], [x], where=lambda bb, key: is_edge(key))
            check(dmax <= POS_TOL, "nodeshift smooth: the four edges follow the line law, within %.1e m of the "
                                   "replica: max %.2e m (%.1f ulp)" % (POS_TOL, dmax, umax))
            # a coordinate no face node moves in is untouched
            still = [c for c in range(3) if vel[c] == 0.0 and bump_vel[c] == 0.0]
            nb, bad = 0, 0
            for key, xo in mesh_out[b][1].items():
                nb += 1
                if any(xo[c] != x_in[b][key][c] for c in still):
                    bad += 1
            check(bad == 0, "nodeshift smooth: the coordinates %s of the %d nodes did not move, bit for bit (%d differ)"
                            % ("".join("xyz"[c] for c in still), nb, bad))
            # every step from the snapshots: the equilibrium and the maximum principle
            if out_diter != 1:
                check(False, "nodeshift smooth: the per-step checks need out-diter = 1 in input.ini")
                return finish()
            meshes = [x_in[b]]
            for s in range(1, steps):
                meshes.append(generate.read_p3d(os.path.join(wd, "OUTPUT", "mesh_%d.p3d" % s))[b][1])
            meshes.append(mesh_out[b][1])
            worst_res, worst_step, worst_mp, spread = 0.0, 0, 0.0, 0.0
            for s in range(1, steps + 1):
                old, new = meshes[s - 1], meshes[s]
                dsp = {key: tuple(new[key][c] - old[key][c] for c in range(3)) for key in new}
                face = [key for key in new if key[dirn] == moving]
                scale = max(max(abs(v) for v in disp(key)) for key in face)
                lo = [min(0.0, min(disp(key)[c] for key in face)) for c in range(3)]
                hi = [max(0.0, max(disp(key)[c] for key in face)) for c in range(3)]
                for key in new:
                    if key[dirn] in (moving, pinned) or is_edge(key):
                        continue                       # prescribed
                    on1 = key[t1] in (0, dims[t1])
                    on2 = key[t2] in (0, dims[t2])
                    use = {dirn: True, t1: not on1, t2: not on2}
                    acc, sw = [0.0, 0.0, 0.0], 0.0
                    for dd in range(3):
                        if not use[dd]:
                            continue
                        for sgn in (-1, 1):
                            q = list(key); q[dd] += sgn; q = tuple(q)
                            e = [old[q][c] - old[key][c] for c in range(3)]
                            w = 1.0 / math.sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2])
                            for c in range(3):
                                acc[c] += w * (dsp[q][c] - dsp[key][c])
                            sw += w
                    res = max(abs(a) for a in acc) / sw / scale
                    if res > worst_res:
                        worst_res, worst_step = res, s
                    for c in range(3):
                        worst_mp = max(worst_mp, lo[c] - dsp[key][c], dsp[key][c] - hi[c])
                    # how far the step's field is from the line law's s/L blend at this node
                    arc = [0.0]
                    for m in range(1, n + 1):
                        p, q = old[line_key(m, key)], old[line_key(m - 1, key)]
                        seg = [p[c] - q[c] for c in range(3)]
                        arc.append(arc[-1] + math.sqrt(seg[0] * seg[0] + seg[1] * seg[1] + seg[2] * seg[2]))
                    m = key[dirn] if side == 2 else n - key[dirn]
                    dl = disp(key)
                    spread = max(spread, max(abs(dsp[key][c] - dl[c] * (arc[m] / arc[n])) for c in range(3)))
            check(worst_res <= SMOOTH_RES_TOL,
                  "nodeshift smooth: GATE every step's field is the spring equilibrium at every solved node, "
                  "residual <= %.1e of the face displacement: max %.2e (step %d; smooth-tolerance %.1e)"
                  % (SMOOTH_RES_TOL, worst_res, worst_step, smooth_tol))
            check(worst_mp <= 10.0 * smooth_tol * scale + POS_TOL,
                  "nodeshift smooth: every node's displacement lies between the prescribed ones (maximum "
                  "principle): worst excess %.2e m" % worst_mp)
            print("        the solved nodes depart from the line law's s/L blend by up to %.2e m per step" % spread)

        else:
            check(False, "verify.py knows no interior '%s'" % interior)
            return finish()

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
