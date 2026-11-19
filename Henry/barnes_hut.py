"""3D Barnes-Hut tree, adapted from our earlier OctreeNode implementation."""
import numpy as np


class Node:
    def __init__(self, center, size, ids, x, m):
        self.center, self.size, self.ids = center, size, ids
        self.mass = m[ids].sum()
        self.com = np.average(x[ids], weights=m[ids], axis=0)
        self.children = []
        if len(ids) <= 1 or size < 1e-10:
            return

        bits = (x[ids] >= center).astype(int)
        octants = bits @ np.array([1, 2, 4])
        for octant in range(8):
            subset = ids[octants == octant]
            if len(subset):
                signs = 2 * np.array([(octant >> bit) & 1 for bit in range(3)]) - 1
                child_center = center + signs * size / 4
                self.children.append(Node(child_center, size / 2, subset, x, m))


def build_tree(x, m):
    low, high = x.min(axis=0), x.max(axis=0)
    size = max(np.max(high - low) * 1.000001, 1e-10)
    return Node((low + high) / 2, size, np.arange(len(x)), x, m)


def from_node(i, node, x, m, eps, theta, G):
    d = node.com - x[i]
    dist = np.linalg.norm(d)
    # Only approximate cells that cannot contain the target particle.
    contains = np.all(np.abs(x[i] - node.center) <= node.size/2 + 1e-14)
    if not node.children:
        ids = node.ids[node.ids != i]
        delta = x[ids] - x[i]
        r2 = np.sum(delta*delta, axis=1) + eps**2
        return G * np.sum(m[ids, None] * delta / r2[:, None]**1.5, axis=0)
    if not contains and dist > 0 and node.size / dist < theta:
        return G * node.mass * d / (dist**2 + eps**2)**1.5
    return sum((from_node(i, c, x, m, eps, theta, G) for c in node.children),
               start=np.zeros(3))


def acceleration(x, m, eps=0.01, theta=0.5, G=1.0, cancel=None):
    root = build_tree(x, m)
    acc = []
    for i in range(len(m)):
        if cancel is not None and cancel.is_set():
            raise InterruptedError('Simulation cancelled')
        acc.append(from_node(i, root, x, m, eps, theta, G))
    return np.array(acc)
