"""Smooth exclusion in the (m_Q, Lambda) plane: log(r) interpolated bilinearly in (m, log10 Lambda) on a fine mesh."""
import array
import math
import numpy as np
import ROOT
import pstyle as P

RCAP = 1e3  # r assigned to grid points without any selected signal
SMOOTH = 3  # passes of the 5x5 smoothing kernel on the log r mesh


def rows(res, key):
    """Per Lambda: the simulated masses and log r (RCAP where no signal survives)."""
    out = {}
    for r in res:
        v = r["best"].get(key) if r["best"] else None
        out.setdefault(r["L"], []).append((r["m"], math.log(min(max(v, 1e-6), RCAP)) if v else math.log(RCAP)))
    return {L: (np.array([p[0] for p in sorted(v)], float), np.array([p[1] for p in sorted(v)])) for L, v in out.items()}


def fill_holes(h, thr):
    """Mark as excluded every non-excluded pocket fully enclosed by the excluded region (not reachable from the mesh edge)."""
    lt = math.log(thr)
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    out = np.array([[h.GetBinContent(i + 1, j + 1) >= lt for j in range(ny)] for i in range(nx)])
    seen = np.zeros_like(out)
    st = [(i, j) for i in range(nx) for j in (0, ny - 1)] + [(i, j) for i in (0, nx - 1) for j in range(ny)]
    while st:  # flood the non-excluded cells from the edge
        i, j = st.pop()
        if 0 <= i < nx and 0 <= j < ny and out[i, j] and not seen[i, j]:
            seen[i, j] = True
            st += [(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)]
    for i, j in zip(*np.nonzero(out & ~seen)):
        h.SetBinContent(int(i) + 1, int(j) + 1, lt - 0.01)


def mesh(res, key, name, nm=300, nl=300, thr=None):
    """TH2D of log r on a fine (m, log10 Lambda) mesh: each Lambda row interpolated over its own masses, then in log10 Lambda."""
    rw = rows(res, key)
    Ls = np.array(sorted(rw), float)
    lL = np.log10(Ls)
    m0 = min(v[0][0] for v in rw.values())
    m1 = max(v[0][-1] for v in rw.values())
    h = ROOT.TH2D(name, "", nm, m0, m1, nl, lL[0], lL[-1])
    h.SetDirectory(0)
    for i in range(1, nm + 1):
        for j in range(1, nl + 1):
            h.SetBinContent(i, j, math.log(RCAP))
    big = math.log(RCAP)
    for i in range(1, nm + 1):
        m = h.GetXaxis().GetBinCenter(i)
        cov = [k for k, L in enumerate(Ls) if rw[L][0][0] <= m <= rw[L][0][-1]]  # rows that simulated this mass
        if not cov:
            continue
        col = np.array([np.interp(m, rw[Ls[k]][0], rw[Ls[k]][1]) for k in cov])
        for j in range(1, nl + 1):
            y = h.GetYaxis().GetBinCenter(j)
            # interpolate in log10 Lambda across the covering rows; outside them, take the row edge only if the gap is one grid step
            h.SetBinContent(i, j, float(np.interp(y, lL[cov], col, left=big, right=big)) if lL[cov][0] <= y <= lL[cov][-1] else big)
    h.Smooth(SMOOTH, "k5a")  # mild 5x5 kernel against point-to-point MC scatter
    if thr:
        fill_holes(h, thr)
    P.keep.append(h)
    return h


def to_lambda(h, name):
    """Copy of the mesh with a log-scale Lambda axis (variable bins), for drawing on a log-y canvas."""
    ny = h.GetNbinsY()
    edges = array.array("d", [10 ** h.GetYaxis().GetBinLowEdge(j) for j in range(1, ny + 2)])
    xe = array.array("d", [h.GetXaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)])
    o = ROOT.TH2D(name, "", h.GetNbinsX(), xe, ny, edges)
    o.SetDirectory(0)
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, ny + 1):
            o.SetBinContent(i, j, h.GetBinContent(i, j))
    P.keep.append(o)
    return o


def contour_graphs(h, thr):
    """TGraphs of the log r = log thr contour (via CONT LIST on a scratch canvas)."""
    c = ROOT.TCanvas("scratch_" + h.GetName(), "", 10, 10)
    hc = h.Clone(h.GetName() + "_c")
    hc.SetContour(1, array.array("d", [math.log(thr)]))
    hc.Draw("CONT LIST")
    c.Update()
    out = []
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    if conts and conts.GetSize() > 0:
        for g in conts.At(0):
            n = g.GetN()
            xs = [g.GetPointX(k) for k in range(n)]
            ys = [10 ** g.GetPointY(k) for k in range(n)]  # back to Lambda
            out.append(P.graph(xs, ys))
    c.Close()
    return out


def draw(res, key, thr, name, fill=None, color=ROOT.kBlack, style=1, width=2):
    """Contour of log r = log thr in the (m, Lambda) plane; optional transparent fill of the excluded side."""
    cname = ROOT.gPad.GetCanvas().GetName()  # gPad is a moving proxy: remember the canvas by name
    h = mesh(res, key, name, thr=thr)
    gs = contour_graphs(h, thr)
    ROOT.gROOT.GetListOfCanvases().FindObject(cname).cd()
    if fill is not None:  # shade the excluded mesh cells, row by row
        lt = math.log(thr)
        ax, ay = h.GetXaxis(), h.GetYaxis()
        for j in range(1, h.GetNbinsY() + 1):
            y0, y1 = 10 ** ay.GetBinLowEdge(j), 10 ** ay.GetBinUpEdge(j)
            i = 1
            while i <= h.GetNbinsX():
                if h.GetBinContent(i, j) < lt:
                    k = i
                    while k + 1 <= h.GetNbinsX() and h.GetBinContent(k + 1, j) < lt:
                        k += 1
                    b = ROOT.TBox(ax.GetBinLowEdge(i), y0, ax.GetBinUpEdge(k), y1)
                    b.SetFillColorAlpha(fill, 0.25)
                    b.SetLineWidth(0)
                    b.Draw()
                    P.keep.append(b)
                    i = k + 1
                else:
                    i += 1
    for g in gs:
        g.SetLineColor(color)
        g.SetLineStyle(style)
        g.SetLineWidth(width)
        g.Draw("L same")
    return gs


def draw_band(res, thr, color, name, hatch=3354):
    """Hatched 68% band of the expected exclusion: cells excluded with the 16% quantile but not with the 84% one."""
    cname = ROOT.gPad.GetCanvas().GetName()
    hlo, hhi = mesh(res, "m1", name + "_m1", thr=thr), mesh(res, "p1", name + "_p1", thr=thr)
    ROOT.gROOT.GetListOfCanvases().FindObject(cname).cd()
    lt = math.log(thr)
    ax, ay = hlo.GetXaxis(), hlo.GetYaxis()
    first = None
    for j in range(1, hlo.GetNbinsY() + 1):
        y0, y1 = 10 ** ay.GetBinLowEdge(j), 10 ** ay.GetBinUpEdge(j)
        inb = [hlo.GetBinContent(i, j) < lt <= hhi.GetBinContent(i, j) for i in range(1, hlo.GetNbinsX() + 1)]
        i = 0
        while i < len(inb):
            if inb[i]:
                k = i
                while k + 1 < len(inb) and inb[k + 1]:
                    k += 1
                b = ROOT.TBox(ax.GetBinLowEdge(i + 1), y0, ax.GetBinUpEdge(k + 1), y1)
                b.SetFillColor(color)
                b.SetFillStyle(hatch)  # screen-aligned hatching, so adjacent boxes join seamlessly
                b.SetLineWidth(0)
                b.Draw()
                P.keep.append(b)
                first = first or b
                i = k + 1
            else:
                i += 1
    return first
