"""Numba-accelerated Barnes–Hut traversal; same octree and softened law as barnes_hut.py."""
import numpy as np
from numba import njit
from .barnes_hut import build_tree


def packed_tree(x, m):
    """Flatten our existing Python octree into arrays for compiled traversal."""
    nodes = []
    children = []
    leaves = []

    def visit(node):
        i = len(nodes)
        nodes.append(node)
        children.append(np.full(8, -1, dtype=np.int32))
        leaves.append(node.ids.copy() if not node.children else np.empty(0, dtype=np.int32))
        for j, child in enumerate(node.children):
            children[i][j] = visit(child)
        return i

    visit(build_tree(x, m))
    starts, counts, ids = [], [], []
    for leaf in leaves:
        starts.append(len(ids))
        counts.append(len(leaf))
        ids.extend(leaf)
    return (
        np.array([n.center for n in nodes]),
        np.array([n.size for n in nodes]),
        np.array([n.com for n in nodes]),
        np.array([n.mass for n in nodes]),
        np.array(children),
        np.array(starts, dtype=np.int32),
        np.array(counts, dtype=np.int32),
        np.array(ids, dtype=np.int32),
    )


@njit(cache=True)
def _traverse(x, m, centers, sizes, coms, masses, children,
              starts, counts, ids, eps, theta, G):
    out = np.zeros_like(x)
    n_nodes = len(sizes)
    for i in range(len(m)):
        stack = np.empty(n_nodes, dtype=np.int32)
        stack[0] = 0
        top = 1
        while top:
            top -= 1
            k = stack[top]
            dx, dy, dz = coms[k, 0]-x[i, 0], coms[k, 1]-x[i, 1], coms[k, 2]-x[i, 2]
            r2 = dx*dx+dy*dy+dz*dz
            contains = (abs(x[i, 0]-centers[k, 0]) <= sizes[k]/2+1e-14 and
                        abs(x[i, 1]-centers[k, 1]) <= sizes[k]/2+1e-14 and
                        abs(x[i, 2]-centers[k, 2]) <= sizes[k]/2+1e-14)
            if counts[k]:
                for q in range(starts[k], starts[k]+counts[k]):
                    j = ids[q]
                    if j != i:
                        a, b, c = x[j, 0]-x[i, 0], x[j, 1]-x[i, 1], x[j, 2]-x[i, 2]
                        fac = G*m[j]/((a*a+b*b+c*c+eps*eps)**1.5)
                        out[i, 0] += fac*a
                        out[i, 1] += fac*b
                        out[i, 2] += fac*c
            elif not contains and r2 > 0 and sizes[k]/np.sqrt(r2) < theta:
                fac = G*masses[k]/((r2+eps*eps)**1.5)
                out[i, 0] += fac*dx
                out[i, 1] += fac*dy
                out[i, 2] += fac*dz
            else:
                for j in range(8):
                    child = children[k, j]
                    if child != -1:
                        stack[top] = child
                        top += 1
    return out


def acceleration(x, m, eps=0.01, theta=0.5, G=1., cancel=None):
    """Numba Barnes–Hut; cancellation is checked between whole force evaluations."""
    if cancel is not None and cancel.is_set():
        raise InterruptedError('Simulation cancelled')
    packed = packed_tree(x, m)
    result = _traverse(x, m, *packed, eps, theta, G)
    if cancel is not None and cancel.is_set():
        raise InterruptedError('Simulation cancelled')
    return result
