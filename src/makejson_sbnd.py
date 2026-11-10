#!/usr/bin/env python3
"""
Convert SBND YZ nonuniformity maps (CzyHist_<plane>_<tpc> TH2s) into the JSON
format read by WireCell's YZMap service, for use by Gen::Scaler in simulation.

Based on Joseph Zennamo's makejson_file_correct.py (ICARUS), adapted for SBND:
  - SBND has one anode per TPC: TPC 0 (East, -x) -> "AnodePlane:apa0",
                                TPC 1 (West, +x) -> "AnodePlane:apa1".
    These keys MUST match wc.tn(anode) in the jsonnet, otherwise
    YZMap::value() returns 0 and Scaler removes all charge.
  - Plane keys are "0", "1", "2" (U, V, Y).
  - YZMap reads grid[binz][biny], so the outer list must run over Z.
    This script assumes Z is on the histogram X axis and Y on the Y axis,
    and prints the axis titles/ranges so you can confirm that.

Usage (on a GPVM with ROOT):
  python3 makejson_sbnd.py INPUT.root -o yzmap_gain_sbnd_run2.json
  python3 makejson_sbnd.py INPUT.root -o out.json --invert   # write 1/C instead of C

Open questions to settle before production (see notes printed at the end):
  * Correction factor C or its inverse 1/C for simulation (--invert)?
  * What value empty / zero bins should take (--fill, default 1.0)?
"""

import argparse
import json
import math

ANODE_KEY = {"0": "AnodePlane:apa0",   # TPC 0, East
             "1": "AnodePlane:apa1"}   # TPC 1, West
PLANES = ["0", "1", "2"]
HIST_TEMPLATE = "CzyHist_{plane}_{tpc}"


def axis_info(axis):
    """Return (nbins, low, high, width, title) and check binning is uniform."""
    n = axis.GetNbins()
    lo, hi = axis.GetXmin(), axis.GetXmax()
    width = (hi - lo) / n
    for b in range(1, n + 1):
        if not math.isclose(axis.GetBinWidth(b), width, rel_tol=1e-6):
            raise ValueError("Non-uniform binning on axis '%s'; YZMap assumes fixed bins"
                             % axis.GetTitle())
    return n, lo, hi, width, axis.GetTitle()


def hist_to_grid(h, fill, invert, ndigits=3):
    """Convert a TH2 (X = Z, Y = Y) into grid[binz][biny], counting replaced bins."""
    nz, ny = h.GetNbinsX(), h.GetNbinsY()
    grid, n_nan, n_zero = [], 0, 0
    for i in range(1, nz + 1):
        row = []
        for j in range(1, ny + 1):
            c = h.GetBinContent(i, j)
            if math.isnan(c) or math.isinf(c):
                c, n_nan = fill, n_nan + 1
            elif c == 0.0:
                c, n_zero = fill, n_zero + 1
            elif invert:
                c = 1.0 / c
            row.append(round(c, ndigits))
        grid.append(row)
    return grid, n_nan, n_zero


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="ROOT file containing CzyHist_<plane>_<tpc>")
    ap.add_argument("-o", "--output", default="yzmap_gain_sbnd.json")
    ap.add_argument("--invert", action="store_true",
                    help="store 1/C instead of C (see open question in docstring)")
    ap.add_argument("--fill", type=float, default=1.0,
                    help="value for NaN/inf/zero bins (default 1.0 = no scaling)")
    args = ap.parse_args()

    import ROOT  # imported here so --help works without ROOT

    f = ROOT.TFile.Open(args.input, "READ")
    if not f or f.IsZombie():
        raise SystemExit("Cannot open %s" % args.input)

    out = {key: {} for key in ANODE_KEY.values()}
    binning = None

    for tpc, anode in ANODE_KEY.items():
        for plane in PLANES:
            name = HIST_TEMPLATE.format(plane=plane, tpc=tpc)
            h = f.Get(name)
            if not h:
                raise SystemExit("Missing histogram %s in %s" % (name, args.input))

            this = (axis_info(h.GetXaxis()), axis_info(h.GetYaxis()))
            if binning is None:
                binning = this
            elif [a[:4] for a in this] != [a[:4] for a in binning]:
                raise SystemExit("%s has different binning from the other maps" % name)

            grid, n_nan, n_zero = hist_to_grid(h, args.fill, args.invert)
            out[anode][plane] = grid
            vals = [v for r in grid for v in r]
            print("%-16s plane %s  <- %-13s  min %.3f  max %.3f  NaN->fill %4d  zero->fill %4d"
                  % (anode, plane, name, min(vals), max(vals), n_nan, n_zero))

    f.Close()

    with open(args.output, "w") as fh:
        json.dump(out, fh)
    print("\nWrote %s (%s)" % (args.output, "1/C" if args.invert else "C as stored"))

    (nz, zlo, zhi, zw, zt), (ny, ylo, yhi, yw, yt) = binning
    print("\nHistogram X axis: '%s'  %d bins  [%g, %g]  -> treated as Z" % (zt, nz, zlo, zhi))
    print("Histogram Y axis: '%s'  %d bins  [%g, %g]  -> treated as Y" % (yt, ny, ylo, yhi))
    print("Check these are Z and Y respectively, in cm.\n")

    print("Matching YZMap config (assuming axes in cm):")
    print("""local yzmap_gain = {
    type: "YZMap", name: "yzmap_gain",
    data: {
        filename: '%s',
        bin_width:  %g*wc.cm,   // Z bin size
        bin_height: %g*wc.cm,   // Y bin size
        zoffset: %g*wc.cm,
        yoffset: %g*wc.cm,
        nbinsz: %d,
        nbinsy: %d,
    }
};""" % (args.output, zw, yw, -zlo, -ylo, nz, ny))


if __name__ == "__main__":
    main()
