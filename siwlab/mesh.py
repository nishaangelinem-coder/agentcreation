"""Structured triangular mesher for planar metal plates with rectangular apertures.
Grid lines are placed on every aperture edge so apertures are represented exactly."""
import numpy as np


def _fill(lines, hmax):
    lines = np.array(sorted(set(np.round(lines, 9))))
    out = [lines[0]]
    for a, b in zip(lines[:-1], lines[1:]):
        n = max(1, int(np.ceil((b - a) / hmax - 1e-9)))
        out.extend(list(a + (b - a) * np.arange(1, n + 1) / n))
    return np.array(out)


def plate_mesh(x0, x1, y0, y1, apertures, hmax):
    """Return (nodes[N,3], tris[M,3]) for the plate [x0,x1]x[y0,y1] with
    rectangular apertures [(ax0, ax1, ay0, ay1), ...] removed. Lengths in metres."""
    xl = [x0, x1]; yl = [y0, y1]
    for (ax0, ax1, ay0, ay1) in apertures:
        xl += [ax0, ax1]; yl += [ay0, ay1]
    xs = _fill(xl, hmax); ys = _fill(yl, hmax)
    nx, ny = len(xs), len(ys)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    nid = np.arange(nx * ny).reshape(nx, ny)
    xc = 0.5 * (xs[1:] + xs[:-1]); yc = 0.5 * (ys[1:] + ys[:-1])
    XC, YC = np.meshgrid(xc, yc, indexing='ij')
    keep = np.ones(XC.shape, bool)
    for (ax0, ax1, ay0, ay1) in apertures:
        keep &= ~((XC > ax0) & (XC < ax1) & (YC > ay0) & (YC < ay1))
    tris = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            if not keep[i, j]:
                continue
            n00, n10, n01, n11 = nid[i, j], nid[i + 1, j], nid[i, j + 1], nid[i + 1, j + 1]
            if (i + j) % 2 == 0:
                tris.append((n00, n10, n11)); tris.append((n00, n11, n01))
            else:
                tris.append((n00, n10, n01)); tris.append((n10, n11, n01))
    tris = np.array(tris)
    used = np.unique(tris)
    remap = -np.ones(nx * ny, int); remap[used] = np.arange(len(used))
    nodes = np.c_[X.ravel()[used], Y.ravel()[used], np.zeros(len(used))]
    return nodes, remap[tris]


def rwg_edges(nodes, tris):
    """Interior edges -> RWG basis functions.
    Returns dict with arrays: edge nodes (E,2), tp, tm (triangle +/-), vp, vm (free vertex
    index), length, triangle areas, centroids."""
    edges = {}
    for t, tri in enumerate(tris):
        for a in range(3):
            n1, n2 = tri[a], tri[(a + 1) % 3]
            key = (min(n1, n2), max(n1, n2))
            edges.setdefault(key, []).append((t, tri[(a + 2) % 3]))
    en, tp, tm, vp, vm = [], [], [], [], []
    for key, lst in edges.items():
        if len(lst) == 2:
            en.append(key); tp.append(lst[0][0]); vp.append(lst[0][1]); tm.append(lst[1][0]); vm.append(lst[1][1])
    en = np.array(en); tp = np.array(tp); tm = np.array(tm); vp = np.array(vp); vm = np.array(vm)
    length = np.linalg.norm(nodes[en[:, 0]] - nodes[en[:, 1]], axis=1)
    v = nodes[tris]
    area = 0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1)
    cen = v.mean(axis=1)
    return dict(en=en, tp=tp, tm=tm, vp=vp, vm=vm, length=length, area=area, cen=cen, nE=len(en))
