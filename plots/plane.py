"""Smooth exclusion in the (m_Q, Lambda) plane: log(r) interpolated bilinearly in (m, log10 Lambda) on a fine mesh."""
import array
import math
import numpy as np
import ROOT
import pstyle

RCAP = 1e3  # r assigned to grid points without any selected signal
SMOOTH = 3  # passes of the 5x5 smoothing kernel on the log r mesh


def rows(res, key):
    """Per Lambda: the simulated masses and log r (RCAP where no signal survives)."""
    out = {}
    for row in res:
        rval = row["best"].get(key) if row["best"] else None
        out.setdefault(row["L"], []).append((row["m"], math.log(min(max(rval, 1e-6), RCAP)) if rval else math.log(RCAP)))
    return {lam: (np.array([pair[0] for pair in sorted(pairs)], float), np.array([pair[1] for pair in sorted(pairs)]))
            for lam, pairs in out.items()}


def fill_holes(hmesh, thr):
    """Mark as excluded every non-excluded pocket fully enclosed by the excluded region (not reachable from the mesh edge)."""
    lt = math.log(thr)
    nx, ny = hmesh.GetNbinsX(), hmesh.GetNbinsY()
    out = np.array([[hmesh.GetBinContent(ix + 1, iy + 1) >= lt for iy in range(ny)] for ix in range(nx)])
    seen = np.zeros_like(out)
    st = [(ix, iy) for ix in range(nx) for iy in (0, ny - 1)] + [(ix, iy) for ix in (0, nx - 1) for iy in range(ny)]
    while st:  # flood the non-excluded cells from the edge
        ix, iy = st.pop()
        if 0 <= ix < nx and 0 <= iy < ny and out[ix, iy] and not seen[ix, iy]:
            seen[ix, iy] = True
            st += [(ix + 1, iy), (ix - 1, iy), (ix, iy + 1), (ix, iy - 1)]
    for ix, iy in zip(*np.nonzero(out & ~seen)):
        hmesh.SetBinContent(int(ix) + 1, int(iy) + 1, lt - 0.01)


def mesh(res, key, name, nm=300, nl=300, thr=None):
    """TH2D of log r on a fine (m, log10 Lambda) mesh: each Lambda row interpolated over its own masses, then in log10 Lambda."""
    rw = rows(res, key)
    Ls = np.array(sorted(rw), float)
    lL = np.log10(Ls)
    m0 = min(row[0][0] for row in rw.values())
    m1 = max(row[0][-1] for row in rw.values())
    hmesh = ROOT.TH2D(name, "", nm, m0, m1, nl, lL[0], lL[-1])
    hmesh.SetDirectory(0)
    for ix in range(1, nm + 1):
        for iy in range(1, nl + 1):
            hmesh.SetBinContent(ix, iy, math.log(RCAP))
    big = math.log(RCAP)
    for ix in range(1, nm + 1):
        mass = hmesh.GetXaxis().GetBinCenter(ix)
        cov = [idx for idx, lam in enumerate(Ls) if rw[lam][0][0] <= mass <= rw[lam][0][-1]]  # rows that simulated this mass
        if not cov:
            continue
        col = np.array([np.interp(mass, rw[Ls[idx]][0], rw[Ls[idx]][1]) for idx in cov])
        for iy in range(1, nl + 1):
            ly = hmesh.GetYaxis().GetBinCenter(iy)
            # interpolate in log10 Lambda across the covering rows; outside them, take the row edge only if the gap is one grid step
            hmesh.SetBinContent(ix, iy, float(np.interp(ly, lL[cov], col, left=big, right=big))
                                if lL[cov][0] <= ly <= lL[cov][-1] else big)
    hmesh.Smooth(SMOOTH, "k5a")  # mild 5x5 kernel against point-to-point MC scatter
    if thr:
        fill_holes(hmesh, thr)
    pstyle.keep.append(hmesh)
    return hmesh


def to_lambda(hmesh, name):
    """Copy of the mesh with a log-scale Lambda axis (variable bins), for drawing on a log-y canvas."""
    ny = hmesh.GetNbinsY()
    edges = array.array("d", [10 ** hmesh.GetYaxis().GetBinLowEdge(iy) for iy in range(1, ny + 2)])
    xe = array.array("d", [hmesh.GetXaxis().GetBinLowEdge(ix) for ix in range(1, hmesh.GetNbinsX() + 2)])
    out = ROOT.TH2D(name, "", hmesh.GetNbinsX(), xe, ny, edges)
    out.SetDirectory(0)
    for ix in range(1, hmesh.GetNbinsX() + 1):
        for iy in range(1, ny + 1):
            out.SetBinContent(ix, iy, hmesh.GetBinContent(ix, iy))
    pstyle.keep.append(out)
    return out


def contour_graphs(hmesh, thr):
    """TGraphs of the log r = log thr contour (via CONT LIST on a scratch canvas)."""
    canv = ROOT.TCanvas("scratch_" + hmesh.GetName(), "", 10, 10)
    hc = hmesh.Clone(hmesh.GetName() + "_c")
    hc.SetContour(1, array.array("d", [math.log(thr)]))
    hc.Draw("CONT LIST")
    canv.Update()
    out = []
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    if conts and conts.GetSize() > 0:
        for cont in conts.At(0):
            npts = cont.GetN()
            xs = [cont.GetPointX(ipt) for ipt in range(npts)]
            ys = [10 ** cont.GetPointY(ipt) for ipt in range(npts)]  # back to Lambda
            out.append(pstyle.graph(xs, ys))
    canv.Close()
    return out


def draw(res, key, thr, name, fill=None, color=ROOT.kBlack, style=1, width=2):
    """Contour of log r = log thr in the (m, Lambda) plane; optional transparent fill of the excluded side."""
    cname = ROOT.gPad.GetCanvas().GetName()  # gPad is a moving proxy: remember the canvas by name
    hmesh = mesh(res, key, name, thr=thr)
    gs = contour_graphs(hmesh, thr)
    ROOT.gROOT.GetListOfCanvases().FindObject(cname).cd()
    if fill is not None:  # shade the excluded mesh cells, row by row
        lt = math.log(thr)
        ax, ay = hmesh.GetXaxis(), hmesh.GetYaxis()
        for iy in range(1, hmesh.GetNbinsY() + 1):
            y0, y1 = 10 ** ay.GetBinLowEdge(iy), 10 ** ay.GetBinUpEdge(iy)
            ix = 1
            while ix <= hmesh.GetNbinsX():
                if hmesh.GetBinContent(ix, iy) < lt:
                    iend = ix
                    while iend + 1 <= hmesh.GetNbinsX() and hmesh.GetBinContent(iend + 1, iy) < lt:
                        iend += 1
                    box = ROOT.TBox(ax.GetBinLowEdge(ix), y0, ax.GetBinUpEdge(iend), y1)
                    box.SetFillColorAlpha(fill, 0.25)
                    box.SetLineWidth(0)
                    box.Draw()
                    pstyle.keep.append(box)
                    ix = iend + 1
                else:
                    ix += 1
    for gr in gs:
        gr.SetLineColor(color)
        gr.SetLineStyle(style)
        gr.SetLineWidth(width)
        gr.Draw("L same")
    return gs


def draw_band(res, thr, color, name, hatch=3354):
    """Hatched 68% band of the expected exclusion: cells excluded with the 16% quantile but not with the 84% one."""
    cname = ROOT.gPad.GetCanvas().GetName()
    hlo, hhi = mesh(res, "m1", name + "_m1", thr=thr), mesh(res, "p1", name + "_p1", thr=thr)
    ROOT.gROOT.GetListOfCanvases().FindObject(cname).cd()
    lt = math.log(thr)
    ax, ay = hlo.GetXaxis(), hlo.GetYaxis()
    first = None
    for iy in range(1, hlo.GetNbinsY() + 1):
        y0, y1 = 10 ** ay.GetBinLowEdge(iy), 10 ** ay.GetBinUpEdge(iy)
        inb = [hlo.GetBinContent(ix, iy) < lt <= hhi.GetBinContent(ix, iy) for ix in range(1, hlo.GetNbinsX() + 1)]
        ix = 0
        while ix < len(inb):
            if inb[ix]:
                iend = ix
                while iend + 1 < len(inb) and inb[iend + 1]:
                    iend += 1
                box = ROOT.TBox(ax.GetBinLowEdge(ix + 1), y0, ax.GetBinUpEdge(iend + 1), y1)
                box.SetFillColor(color)
                box.SetFillStyle(hatch)  # screen-aligned hatching, so adjacent boxes join seamlessly
                box.SetLineWidth(0)
                box.Draw()
                pstyle.keep.append(box)
                first = first or box
                ix = iend + 1
            else:
                ix += 1
    return first
