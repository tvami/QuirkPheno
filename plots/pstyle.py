# Paper style: ROOT + cmsstyle, square canvases, no experiment label
import array
import ROOT
import cmsstyle as CMS

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
CMS.SetExtraText("")
CMS.SetCmsText("")  # pheno study: no experiment label
COLORS = [ROOT.kAzure + 1, ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kRed + 1, ROOT.kViolet + 1, ROOT.kGray + 2]
keep = []  # protect ROOT objects from GC
OUT = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "figs")


def label(text="Geant4 11.4.2", energy=None):
    CMS.SetLumi(None, run=text)
    if energy is None:
        CMS.SetEnergy(0, unit="")
    else:
        CMS.SetEnergy(energy)


def canvas(name, xlo, xhi, ylo, yhi, xt, yt, logy=False, logx=False, iPos=11, yoff=1.5):
    c = CMS.cmsCanvas(name, xlo, xhi, ylo, yhi, xt, yt, square=True, iPos=iPos, extraSpace=0.01, yTitOffset=yoff)
    c.SetLogy(logy)
    c.SetLogx(logx)
    h = CMS.GetCmsCanvasHist(c)
    h.GetXaxis().SetNdivisions(505)
    if logx:
        h.GetXaxis().SetMoreLogLabels(False)
    keep.append(c)
    return c


def dicanvas(name, xlo, xhi, ylo, yhi, rlo, rhi, xt, yt, rt, logy=False):
    c = CMS.cmsDiCanvas(name, xlo, xhi, ylo, yhi, rlo, rhi, xt, yt, rt, square=True, iPos=11, extraSpace=0.02)
    if logy:
        c.cd(1).SetLogy(True)
    for pad in (1, 2):
        CMS.GetCmsCanvasHist(c.cd(pad)).GetXaxis().SetNdivisions(505)
    c.cd(1)
    keep.append(c)
    return c


def legend(x1=0.55, y1=0.72, x2=0.93, y2=0.88, size=0.034):
    l = CMS.cmsLeg(x1, y1, x2, y2, textSize=size)
    keep.append(l)
    return l


def note(text, x=0.18, y=0.70, size=0.032, align=11):
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(42)
    t.SetTextSize(size)
    t.SetTextAlign(align)
    t.DrawLatex(x, y, text)
    keep.append(t)


def graph(xs, ys, eys=None):
    a = lambda v: array.array('d', [float(t) for t in v])
    n = len(xs)
    g = ROOT.TGraphErrors(n, a(xs), a(ys), a([0.] * n), a(eys if eys is not None else [0.] * n))
    keep.append(g)
    return g


def hist(name, vals, nb, lo, hi, weights=None):
    h = ROOT.TH1D(name, "", nb, lo, hi)
    h.SetDirectory(0)
    for i, v in enumerate(vals):
        h.Fill(min(max(v, lo + 1e-9 * (hi - lo)), hi - 1e-9 * (hi - lo)), 1. if weights is None else weights[i])
    keep.append(h)
    return h


def save(c, stem):
    import os
    c.SaveAs(os.path.join(OUT, stem + ".pdf"))
    c.SaveAs(os.path.join(OUT, stem + ".png"))
