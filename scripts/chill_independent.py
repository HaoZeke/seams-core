"""Independent CHILL and CHILL+ reference.

The spherical harmonics are SymPy's Y_3^m, lambdified and evaluated with
NumPy. Neighbour searches, correlations and the eight label rules are
written here and are not taken from the d-SEAMS sources.

A is the shell gate:
  shell    ice-like labels require exactly four neighbours inside the cutoff
  nearest  the four nearest neighbours inside the cutoff are scored even
           when a fifth neighbour also sits inside the cutoff

B is the set averaged into q_3 of every molecule:
  shell    every neighbour inside the cutoff
  nearest  the four nearest neighbours inside the cutoff

C is the interfacial-neighbour test:
  True     Nguyen Table 1 for CHILL+, Moore's wording for CHILL
  False    two staggered bonds, or three staggered and no eclipsed, with no
           look at the neighbour's own staggered count

Bonds are the neighbours selected by A. q_3 is the average selected by B.
A molecule with fewer than four scored bonds is liquid under every variant.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
CUTOFF = 3.5
LABELS = (
    "cubic",
    "hexagonal",
    "interfacial",
    "clathrate",
    "interclathrate",
    "liquid",
)

def _harmonic_expr(degree, order, theta, phi):
    """Condon-Shortley Y_l^m from SymPy's associated Legendre function."""
    absm = abs(order)
    cos_theta = sp.cos(theta)
    norm = sp.sqrt(
        (2 * degree + 1)
        / (4 * sp.pi)
        * sp.factorial(degree - absm)
        / sp.factorial(degree + absm)
    )
    radial = sp.assoc_legendre(degree, absm, cos_theta)
    phase = sp.exp(sp.I * order * phi)
    sign = 1 if order >= 0 else (-1) ** absm
    return sp.simplify(sign * norm * radial * phase)


_THETA, _PHI = sp.symbols("theta phi", real=True)
_Y_EXPRS = [_harmonic_expr(3, m, _THETA, _PHI) for m in range(-3, 4)]
_Y_FUNC = sp.lambdify((_THETA, _PHI), _Y_EXPRS, modules="numpy")


def y3m(theta, phi):
    """Evaluate Y_3^m for m = -3..3. Returns a complex array of shape (7, ...)."""
    values = _Y_FUNC(np.asarray(theta, dtype=float), np.asarray(phi, dtype=float))
    return np.stack([np.asarray(v, dtype=np.complex128) for v in values], axis=0)


def minimum_image(pos, box):
    """Orthorhombic minimum-image deltas, shape (N, N, 3), and distances."""
    pos = np.asarray(pos, dtype=float)
    box = np.asarray(box, dtype=float).reshape(3)
    delta = pos[:, None, :] - pos[None, :, :]
    delta -= box * np.round(delta / box)
    dist = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
    np.fill_diagonal(dist, np.inf)
    return delta, dist


def shells(dist, cutoff):
    """Neighbours with distance strictly below cutoff, nearest first.

    Equal distances keep the lower atom index. The cutoff is open at the
    top because the OVITO manual describes an O-O distance less than 3.5 A,
    and Nguyen's footnote says neighbours within 3.5 A. A neighbour that
    lands on the cutoff exactly is not in this shell.
    """
    n = dist.shape[0]
    rows = []
    for i in range(n):
        js = np.flatnonzero(dist[i] < cutoff)
        order = np.argsort(dist[i, js], kind="mergesort")
        rows.append(js[order])
    return rows


def _select(row, mode):
    if mode == "nearest":
        return row[:4]
    if mode == "shell":
        return row
    raise ValueError(mode)


def q3(delta, dist, rows, mode):
    """q_3m averaged over the neighbour set named by mode. Shape (N, 7)."""
    n = len(rows)
    q = np.zeros((n, 7), dtype=np.complex128)
    for i, row in enumerate(rows):
        chosen = _select(row, mode)
        if chosen.size == 0:
            continue
        dr = delta[i, chosen]
        r = dist[i, chosen]
        theta = np.arccos(np.clip(dr[:, 2] / r, -1.0, 1.0))
        phi = np.arctan2(dr[:, 1], dr[:, 0])
        harmonics = y3m(theta, phi)
        q[i] = harmonics.mean(axis=1)
    return q


def correlation(qi, qj):
    """Real part of the normalised inner product in Nguyen equation 2."""
    ni = float(np.vdot(qi, qi).real)
    nj = float(np.vdot(qj, qj).real)
    if ni <= 0.0 or nj <= 0.0:
        return np.nan
    return float(np.vdot(qj, qi).real / np.sqrt(ni * nj))


def bond_class(value, scheme):
    if not np.isfinite(value):
        return "out"
    if value <= -0.8:
        return "S"
    if scheme == "plus":
        if -0.35 <= value <= 0.25:
            return "E"
    elif scheme == "chill":
        if -0.2 <= value <= -0.05:
            return "E"
    else:
        raise ValueError(scheme)
    return "out"


def classify(pos, box, a_mode, b_mode, interfacial, scheme="plus", cutoff=CUTOFF):
    """Return labels, per-atom c_ij values, shell sizes and staggered counts.

    c_ij is reported for the bonds selected by A, in nearest-first order.
    """
    pos = np.asarray(pos, dtype=float)
    delta, dist = minimum_image(pos, box)
    rows = shells(dist, cutoff)
    q = q3(delta, dist, rows, b_mode)
    n = len(rows)
    cij = [np.empty(0, dtype=float) for _ in range(n)]
    kinds = [[] for _ in range(n)]
    shell_n = np.array([len(row) for row in rows], dtype=int)
    for i, row in enumerate(rows):
        chosen = _select(row, a_mode)
        if a_mode == "nearest":
            chosen = chosen[:4]
        values = np.empty(chosen.size, dtype=float)
        tagged = []
        for k, j in enumerate(chosen):
            values[k] = correlation(q[i], q[j])
            tagged.append(bond_class(values[k], scheme))
        cij[i] = values
        kinds[i] = tagged
    staggered = np.array([sum(k == "S" for k in tagged) for tagged in kinds], dtype=int)
    eclipsed = np.array([sum(k == "E" for k in tagged) for tagged in kinds], dtype=int)
    labels = []
    for i, row in enumerate(rows):
        chosen = _select(row, a_mode)
        n_bonds = int(chosen.size)
        eligible = n_bonds == 4
        if a_mode == "shell" and shell_n[i] != 4:
            eligible = False
        if not eligible:
            labels.append("liquid")
            continue
        ns = int(staggered[i])
        ne = int(eclipsed[i])
        label = "liquid"
        if scheme == "plus" and ne == 0 and ns == 4:
            label = "cubic"
        elif scheme == "plus" and ne == 1 and ns == 3:
            label = "hexagonal"
        elif scheme == "chill" and ns == 4:
            label = "cubic"
        elif scheme == "chill" and ne == 1 and ns == 3:
            label = "hexagonal"
        else:
            interfacial_hit = False
            if interfacial:
                neigh = [_select(rows[int(j)], a_mode) for j in chosen]
                neigh_s = [int(staggered[int(j)]) for j in chosen]
                if ns == 2 and any(s > 2 for s in neigh_s):
                    interfacial_hit = True
                if scheme == "plus" and ns == 3 and ne == 0 and any(s > 1 for s in neigh_s):
                    interfacial_hit = True
                if scheme == "chill" and ns == 3 and ne == 0 and any(s == 2 for s in neigh_s):
                    interfacial_hit = True
                del neigh
            else:
                if ns == 2 or (ns == 3 and ne == 0):
                    interfacial_hit = True
            if interfacial_hit:
                label = "interfacial"
            elif scheme == "plus" and ne == 4 and ns == 0:
                label = "clathrate"
            elif scheme == "plus" and ne == 3:
                label = "interclathrate"
        labels.append(label)
    return {
        "label": np.array(labels),
        "cij": cij,
        "shell_n": shell_n,
        "staggered": staggered,
        "eclipsed": eclipsed,
        "q": q,
    }


def variants(scheme="plus"):
    out = []
    for a_mode in ("shell", "nearest"):
        for b_mode in ("shell", "nearest"):
            for interfacial in (True, False):
                name = "A-{0}_B-{1}_C-{2}".format(
                    a_mode, b_mode, "on" if interfacial else "off"
                )
                out.append((name, a_mode, b_mode, interfacial, scheme))
    return out


def label_all(pos, box, scheme="plus", cutoff=CUTOFF):
    found = {}
    for name, a_mode, b_mode, interfacial, scheme_name in variants(scheme):
        found[name] = classify(pos, box, a_mode, b_mode, interfacial, scheme_name, cutoff)
    return found
