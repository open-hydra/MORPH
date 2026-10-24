#!/usr/bin/env python3
"""The mesh of the driver test: three conforming blocks of warped hexahedra.

One global structured lattice of 10 x 6 x 4 cells, graded in i and j and
distorted by a smooth map so that no face is planar, is cut into three blocks
of different sizes:

    block 1   i 0..4    j 0..6   k 0..4     4 x 6 x 4 cells
    block 2   i 4..10   j 0..2   k 0..4     6 x 2 x 4
    block 3   i 4..10   j 2..6   k 0..4     6 x 4 x 4

Block 1's i-high face meets block 2's and block 3's i-low faces; block 2's
j-high face meets block 3's j-low face. A node on an interface is computed once
from its lattice indices, so the blocks that share it hold the same bits.

Writes PLOT3D ASCII multiblock with ONE value per line: ORION reads each value
with its own list-directed READ, which starts a new record. repr() gives the
shortest decimal that reads back to the same double.

Python 3.6, standard library only.
"""
import math
import os
import sys

LX, LY, LZ = 0.4, 0.3, 0.2
NI, NJ, NK = 10, 6, 4

#: (name, i range, j range, k range) in lattice node indices, inclusive
BLOCKS = [
    ("block1", (0, 4), (0, 6), (0, 4)),
    ("block2", (4, 10), (0, 2), (0, 4)),
    ("block3", (4, 10), (2, 6), (0, 4)),
]

#: pairs of block faces that share their nodes: (block a, i,j,k ranges on a), same for b
#: given as lattice node ranges; every lattice node in the range belongs to both blocks
INTERFACES = [
    (0, 1, (4, 4), (0, 2), (0, 4)),
    (0, 2, (4, 4), (2, 6), (0, 4)),
    (1, 2, (4, 10), (2, 2), (0, 4)),
]


def lattice_node(I, J, K):
    """Position of lattice node (I, J, K): graded, then smoothly distorted."""
    xi, eta, zeta = I / NI, J / NJ, K / NK
    x = LX * (xi + 0.15 * math.sin(math.pi * xi) / math.pi)
    y = LY * eta ** 1.3
    z = LZ * zeta
    X = x + 0.004 * math.sin(2.0 * math.pi * y / LY) * math.sin(math.pi * z / LZ)
    Y = y + 0.003 * math.sin(2.0 * math.pi * x / LX) * math.cos(math.pi * z / LZ)
    Z = z + 0.002 * math.sin(math.pi * x / LX) * math.sin(2.0 * math.pi * y / LY)
    return (X, Y, Z)


def block_nodes(b):
    """{(i, j, k) local: (x, y, z)} of block b (0-based)."""
    _, (i0, i1), (j0, j1), (k0, k1) = BLOCKS[b]
    return {(I - i0, J - j0, K - k0): lattice_node(I, J, K)
            for K in range(k0, k1 + 1)
            for J in range(j0, j1 + 1)
            for I in range(i0, i1 + 1)}


def dims(b):
    """Cells per direction of block b."""
    _, (i0, i1), (j0, j1), (k0, k1) = BLOCKS[b]
    return (i1 - i0, j1 - j0, k1 - k0)


def write_p3d(path, blocks):
    """blocks: list of (dims, {(i,j,k): (x,y,z)}). PLOT3D ASCII, one value per line."""
    with open(path, "w") as f:
        f.write("%d\n" % len(blocks))
        for (ni, nj, nk), _ in blocks:
            f.write("%d %d %d\n" % (ni + 1, nj + 1, nk + 1))
        for (ni, nj, nk), nodes in blocks:
            for d in range(3):
                for k in range(nk + 1):
                    for j in range(nj + 1):
                        for i in range(ni + 1):
                            f.write(repr(nodes[(i, j, k)][d]) + "\n")


def read_p3d(path):
    """[(dims, {(i,j,k): (x,y,z)})] from a PLOT3D ASCII multiblock file (3-D)."""
    with open(path) as f:
        tok = f.read().split()
    pos = 0
    nb = int(tok[pos]); pos += 1
    shape = []
    for _ in range(nb):
        shape.append(tuple(int(t) - 1 for t in tok[pos:pos + 3])); pos += 3
    out = []
    for ni, nj, nk in shape:
        n = (ni + 1) * (nj + 1) * (nk + 1)
        comp = []
        for d in range(3):
            comp.append([float(t) for t in tok[pos:pos + n]]); pos += n
        nodes = {}
        m = 0
        for k in range(nk + 1):
            for j in range(nj + 1):
                for i in range(ni + 1):
                    nodes[(i, j, k)] = (comp[0][m], comp[1][m], comp[2][m])
                    m += 1
        out.append(((ni, nj, nk), nodes))
    if pos != len(tok):
        raise ValueError("%s: %d values left over" % (path, len(tok) - pos))
    return out


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join("INPUT", "mesh.p3d")
    d = os.path.dirname(out)
    if d:
        os.makedirs(d, exist_ok=True)
    write_p3d(out, [(dims(b), block_nodes(b)) for b in range(len(BLOCKS))])
    print("generate.py: %d blocks, %s cells -> %s"
          % (len(BLOCKS), " + ".join("x".join(str(n) for n in dims(b)) for b in range(len(BLOCKS))), out))


if __name__ == "__main__":
    main()
