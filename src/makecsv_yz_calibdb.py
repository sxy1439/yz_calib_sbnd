#!/usr/bin/env python3
"""
Export SBND YZ nonuniformity correction factors (CzyHist_<plane>_<tpc> TH2s
written by yz_median_planes_sce.C) to a CSV for the SBND Calibration Database
table "tpc_yz_correction", as documented in sbnd_calibration_database.md:

  variable  type     description
  tpc       text     W = west TPC, E = east TPC
  ybin      integer  0 to 79
  zbin      integer  0 to 99
  ylow      real     lower y edge [cm]
  yhigh     real     upper y edge [cm]
  zlow      real     lower z edge [cm]
  zhigh     real     upper z edge [cm]
  scale     real     YZ non-uniformity scale (dQ/dx_corr = dQ/dx * scale)

The table holds a single map per TPC (2 x 80 x 100 = 16000 rows), so one plane
is exported per CSV (--plane, default 2 = collection).

CSV layout (no header row, like the DB examples), first column is the channel:

  channel,tpc,ybin,zbin,ylow,yhigh,zlow,zhigh,scale

  channel = itpc*(ny*nz) + ybin*nz + zbin, with itpc 0 = E (x<0), 1 = W (x>=0)
  (the channel has no physical meaning, it only has to be unique)

Bins without enough statistics are 0 in CzyHist; they are written as --fill
(default 0, which NormalizeYZ treats as "no correction").

By default the CSV is named after the input file and written to a csv_files/
folder next to it (created if needed), so
  output_files/run_2_sample/NAME.root -> output_files/run_2_sample/csv_files/NAME.csv
Planes other than the default get a suffix: NAME_plane0.csv. Use -o to override.

Usage:
  python3 makecsv_yz_calibdb.py INPUT.root                                     # collection plane
  python3 makecsv_yz_calibdb.py INPUT.root --plane 0
  python3 makecsv_yz_calibdb.py INPUT.root -o check.csv --header               # inspection only

Works with PyROOT if available, otherwise with uproot (pip install uproot).
"""

import argparse
import math
import os

TPC_LABEL = {0: "E", 1: "W"}   # index as in CzyHist_<plane>_<tpc> and NormalizeYZ
HIST_TEMPLATE = "CzyHist_{plane}_{tpc}"
COLUMNS = ["channel", "tpc", "ybin", "zbin", "ylow", "yhigh", "zlow", "zhigh", "scale"]


def read_hist_pyroot(path, name):
    """Return (values[iz][iy], zedges, yedges) for TH2 `name`, X axis = z, Y axis = y."""
    import ROOT
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise SystemExit("Cannot open %s" % path)
    h = f.Get(name)
    if not h:
        raise SystemExit("Missing histogram %s in %s" % (name, path))
    xa, ya = h.GetXaxis(), h.GetYaxis()
    zedges = [xa.GetBinLowEdge(i) for i in range(1, xa.GetNbins() + 2)]
    yedges = [ya.GetBinLowEdge(j) for j in range(1, ya.GetNbins() + 2)]
    vals = [[h.GetBinContent(i, j) for j in range(1, ya.GetNbins() + 1)]
            for i in range(1, xa.GetNbins() + 1)]
    f.Close()
    return vals, zedges, yedges


def read_hist_uproot(path, name):
    """Same as read_hist_pyroot, using uproot."""
    import uproot
    f = uproot.open(path)
    if name not in f:
        raise SystemExit("Missing histogram %s in %s" % (name, path))
    vals, zedges, yedges = f[name].to_numpy()
    return vals.tolist(), list(zedges), list(yedges)


def read_hist(path, name):
    try:
        import ROOT  # noqa: F401
        return read_hist_pyroot(path, name)
    except ImportError:
        return read_hist_uproot(path, name)


def default_output(input_path, plane):
    """csv_files/<input name>.csv next to the input file, with a suffix for non-collection planes."""
    indir, base = os.path.split(input_path)
    stem = os.path.splitext(base)[0]
    suffix = "" if plane == 2 else "_plane%d" % plane
    return os.path.join(indir, "csv_files", stem + suffix + ".csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="ROOT file from yz_median_planes_sce.C (contains CzyHist_<plane>_<tpc>)")
    ap.add_argument("-o", "--output", default=None,
                    help="output CSV (default: csv_files/<input name>.csv next to the input file)")
    ap.add_argument("--plane", type=int, choices=[0, 1, 2], default=2,
                    help="plane to export: 0 = U, 1 = V, 2 = Y/collection (default 2)")
    ap.add_argument("--fill", type=float, default=0.0,
                    help="value written for empty/NaN bins (default 0 = no correction)")
    ap.add_argument("--precision", type=int, default=6, help="decimals for scale (default 6)")
    ap.add_argument("--header", action="store_true",
                    help="write a header row (not for the DB upload, only for inspection)")
    args = ap.parse_args()

    hists = {itpc: read_hist(args.input, HIST_TEMPLATE.format(plane=args.plane, tpc=itpc))
             for itpc in TPC_LABEL}

    # both TPC maps must share the same binning, otherwise one channel formula can't describe them
    _, zedges, yedges = hists[0]
    _, ze, ye = hists[1]
    if len(ze) != len(zedges) or len(ye) != len(yedges) or \
       not all(math.isclose(a, b, abs_tol=1e-6) for a, b in zip(ze + ye, zedges + yedges)):
        raise SystemExit("East and West maps have different binning")
    nz, ny = len(zedges) - 1, len(yedges) - 1

    if args.output is None:
        args.output = default_output(args.input, args.plane)
    outdir = os.path.dirname(args.output)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    nrows = 0
    with open(args.output, "w") as fh:
        if args.header:
            fh.write(",".join(COLUMNS) + "\n")
        for itpc, label in TPC_LABEL.items():
            vals = hists[itpc][0]
            n_empty, filled = 0, []
            for iy in range(ny):
                for iz in range(nz):
                    c = vals[iz][iy]
                    if math.isnan(c) or math.isinf(c) or c == 0.0:
                        c, n_empty = args.fill, n_empty + 1
                    else:
                        filled.append(c)
                    channel = itpc * (ny * nz) + iy * nz + iz
                    fh.write("%d,%s,%d,%d,%g,%g,%g,%g,%.*f\n"
                             % (channel, label, iy, iz,
                                yedges[iy], yedges[iy + 1], zedges[iz], zedges[iz + 1],
                                args.precision, c))
                    nrows += 1
            print("tpc %s <- %s  filled %5d  empty->%g %5d  min %.4f  max %.4f"
                  % (label, HIST_TEMPLATE.format(plane=args.plane, tpc=itpc), len(filled),
                     args.fill, n_empty, min(filled, default=float("nan")),
                     max(filled, default=float("nan"))))

    print("\nWrote %d rows (plane %d) to %s%s" % (nrows, args.plane, args.output,
                                                 " (with header)" if args.header else ""))
    print("Binning: z %d bins [%g, %g] cm, y %d bins [%g, %g] cm"
          % (nz, zedges[0], zedges[-1], ny, yedges[0], yedges[-1]))
    if (ny, nz) != (80, 100):
        print("WARNING: tpc_yz_correction is documented with 80 y bins x 100 z bins")
    print("write_data.py variable list: %s" % ",".join(COLUMNS[1:]))


if __name__ == "__main__":
    main()
