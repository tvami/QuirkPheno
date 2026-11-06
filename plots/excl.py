"""Excluded mass interval from r(m) = sigma95 / sigma on a mass grid (log-linear interpolation)."""
import math
import numpy as np


def _x(m0, m1, r0, r1, thr):
    return m0 + (m1 - m0) * math.log(thr / r0) / math.log(r1 / r0)


def interval(ms, rs, thr):
    """(lo, hi) of the highest-mass contiguous run with r < thr; lo = lowest grid mass if the run starts there; (0, 0) if none."""
    o = np.argsort(ms)
    ms, rs = np.array(ms, float)[o], np.array(rs, float)[o]
    ok = np.nonzero(rs < thr)[0]
    if len(ok) == 0 or len(ms) < 2:
        return 0., 0.
    i = ok[-1]
    hi = _x(ms[i], ms[i + 1], rs[i], rs[i + 1], thr) if i + 1 < len(ms) else ms[-1]
    j = i
    while j > 0 and rs[j - 1] < thr:
        j -= 1
    lo = _x(ms[j - 1], ms[j], rs[j - 1], rs[j], thr) if j > 0 else ms[0]
    return lo, hi
