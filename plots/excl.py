"""Excluded mass interval from r(m) = sigma95 / sigma on a mass grid (log-linear interpolation)."""
import math
import numpy as np


def _x(m0, m1, r0, r1, thr):
    return m0 + (m1 - m0) * math.log(thr / r0) / math.log(r1 / r0)


def interval(ms, rs, thr):
    """(lo, hi) of the highest-mass contiguous run with r < thr; lo = lowest grid mass if the run starts there; (0, 0) if none."""
    order = np.argsort(ms)
    ms, rs = np.array(ms, float)[order], np.array(rs, float)[order]
    ok = np.nonzero(rs < thr)[0]
    if len(ok) == 0 or len(ms) < 2:
        return 0., 0.
    ihi = ok[-1]
    hi = _x(ms[ihi], ms[ihi + 1], rs[ihi], rs[ihi + 1], thr) if ihi + 1 < len(ms) else ms[-1]
    ilo = ihi
    while ilo > 0 and rs[ilo - 1] < thr:
        ilo -= 1
    lo = _x(ms[ilo - 1], ms[ilo], rs[ilo - 1], rs[ilo], thr) if ilo > 0 else ms[0]
    return lo, hi
