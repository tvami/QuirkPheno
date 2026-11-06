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
__import__("os").makedirs(OUT, exist_ok=True)


def label(text="Geant4 11.4.2", energy=None):
    CMS.SetLumi(None, run=text)
    if energy is None:
        CMS.SetEnergy(0, unit="")
    else:
        CMS.SetEnergy(energy)


def canvas(name, xlo, xhi, ylo, yhi, xt, yt, logy=False, logx=False, iPos=11, yoff=1.5):
    canv = CMS.cmsCanvas(name, xlo, xhi, ylo, yhi, xt, yt, square=True, iPos=iPos, extraSpace=0.01, yTitOffset=yoff)
    canv.SetLogy(logy)
    canv.SetLogx(logx)
    frame = CMS.GetCmsCanvasHist(canv)
    frame.GetXaxis().SetNdivisions(505)
    if logx:
        frame.GetXaxis().SetMoreLogLabels(False)
    keep.append(canv)
    return canv


def dicanvas(name, xlo, xhi, ylo, yhi, rlo, rhi, xt, yt, rt, logy=False):
    canv = CMS.cmsDiCanvas(name, xlo, xhi, ylo, yhi, rlo, rhi, xt, yt, rt, square=True, iPos=11, extraSpace=0.02)
    if logy:
        canv.cd(1).SetLogy(True)
    for pad in (1, 2):
        CMS.GetCmsCanvasHist(canv.cd(pad)).GetXaxis().SetNdivisions(505)
    canv.cd(1)
    keep.append(canv)
    return canv


def legend(x1=0.55, y1=0.72, x2=0.93, y2=0.88, size=0.034):
    leg = CMS.cmsLeg(x1, y1, x2, y2, textSize=size)
    keep.append(leg)
    return leg


def note(text, x=0.18, y=0.70, size=0.032, align=11):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(size)
    latex.SetTextAlign(align)
    latex.DrawLatex(x, y, text)
    keep.append(latex)


def graph(xs, ys, eys=None):
    to_arr = lambda vals: array.array('d', [float(val) for val in vals])
    npts = len(xs)
    gr = ROOT.TGraphErrors(npts, to_arr(xs), to_arr(ys), to_arr([0.] * npts), to_arr(eys if eys is not None else [0.] * npts))
    keep.append(gr)
    return gr


def hist(name, vals, nb, lo, hi, weights=None):
    histo = ROOT.TH1D(name, "", nb, lo, hi)
    histo.SetDirectory(0)
    for idx, val in enumerate(vals):
        histo.Fill(min(max(val, lo + 1e-9 * (hi - lo)), hi - 1e-9 * (hi - lo)), 1. if weights is None else weights[idx])
    keep.append(histo)
    return histo


def save(canv, stem):
    import os
    canv.SaveAs(os.path.join(OUT, stem + ".pdf"))
    canv.SaveAs(os.path.join(OUT, stem + ".png"))
