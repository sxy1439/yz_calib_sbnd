// compare_czy_percent_diff.C
//
// For two YZ-map ROOT files (output of yz_median_planes_sce.C), compare the
// correction-factor maps CzyHist_<plane>_<tpc> bin by bin and histogram the
// percentage difference:
//
//     % diff = (Czy_new - Czy_ref) / Czy_ref * 100
//
// and the signed difference (no units, Czy is a ratio):
//
//     diff = Czy_new - Czy_ref
//
// One entry per YZ pixel. Pixels that are empty (0) in either file are skipped.
// One histogram per plane and TPC (6 for % diff, 6 for diff). Each set of 6 shares its binning.
//
// Usage (from src/):
//   root -l -b -q 'compare_czy_percent_diff.C("<ref file>", "<new file>")'
//
// Optional arguments:
//   label_ref, label_new : legend/title text (default: taken from the "(...)" in the file names)
//   range                : x-axis goes from -range to +range in %. <= 0 means automatic.
//   nbins                : number of bins (default 100)
//   out_dir              : where the ROOT file and PNGs go
//   range_diff           : like range, but for the Czy_new - Czy_ref plots. <= 0 means automatic.
//
// Output:
//   <out_dir>/czy_pdiff_<new>_vs_<ref>.root          12 histograms (pdiff_* and diff_*)
//   <out_dir>/czy_pdiff_<new>_vs_<ref>/*.png          one PNG per plane/TPC + one summary PNG, for each set

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

// Text inside the last "(...)" of a file name, e.g. "1000usEast1000usWest".
TString pdiff_tag_from_file(const TString& path)
{
  TString base = gSystem->BaseName(path);
  int open  = base.Last('(');
  int close = base.Last(')');
  if (open >= 0 && close > open) return base(open + 1, close - open - 1);
  base.ReplaceAll(".root", "");
  return base;
}

// "1000usEast1000usWest" -> "1000 #mus East / 1000 #mus West", or "1000 #mus" if East == West.
TString pdiff_label_from_tag(TString tag)
{
  int iw = tag.Index("West");
  int ie = tag.Index("East");
  if (ie > 0 && iw > ie) {
    TString east = tag(0, ie);
    TString west = tag(ie + 4, iw - ie - 4);
    east.ReplaceAll("us", " #mus"); west.ReplaceAll("us", " #mus");
    east.ReplaceAll("ms", " ms");   west.ReplaceAll("ms", " ms");
    if (!east.Contains(" ")) east.ReplaceAll("s", " s");
    if (!west.Contains(" ")) west.ReplaceAll("s", " s");
    if (east == west) return east;
    return east + " East / " + west + " West";
  }
  return tag;
}

const char* pdiff_plane_name(int plane)
{
  if (plane == 0) return "Plane 0 (Induction 1)";
  if (plane == 1) return "Plane 1 (Induction 2)";
  if (plane == 2) return "Plane 2 (Collection)";
  return "Plane ?";
}

const char* pdiff_tpc_name(int tpc) { return tpc == 0 ? "East TPC" : "West TPC"; }

// Percent differences for one plane/TPC. Skips pixels empty in either map.
std::vector<double> pdiff_values(TH2* h_ref, TH2* h_new, int& n_skipped)
{
  std::vector<double> v;
  n_skipped = 0;
  for (int ix = 1; ix <= h_ref->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h_ref->GetNbinsY(); iy++) {
      double ref = h_ref->GetBinContent(ix, iy);
      double nw  = h_new->GetBinContent(ix, iy);
      if (ref == 0 || nw == 0 || !std::isfinite(ref) || !std::isfinite(nw)) { n_skipped++; continue; }
      v.push_back((nw - ref) / ref * 100.0);
    }
  }
  return v;
}

// Signed differences Czy_new - Czy_ref. Skips the same pixels as pdiff_values.
std::vector<double> diff_values(TH2* h_ref, TH2* h_new)
{
  std::vector<double> v;
  for (int ix = 1; ix <= h_ref->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h_ref->GetNbinsY(); iy++) {
      double ref = h_ref->GetBinContent(ix, iy);
      double nw  = h_new->GetBinContent(ix, iy);
      if (ref == 0 || nw == 0 || !std::isfinite(ref) || !std::isfinite(nw)) continue;
      v.push_back(nw - ref);
    }
  }
  return v;
}

// A "nice" round number >= x (1, 2, 5, 10, 20, 50, ...).
double pdiff_nice_ceil(double x)
{
  if (x <= 0) return 1.0;
  double p = std::pow(10.0, std::floor(std::log10(x)));
  for (double m : {1.0, 2.0, 5.0, 10.0}) if (m * p >= x) return m * p;
  return 10.0 * p;
}

void compare_czy_percent_diff(const char* ref_file,
                              const char* new_file,
                              const char* label_ref = "",
                              const char* label_new = "",
                              double range = -1,
                              int nbins = 100,
                              const char* out_dir = "../output_files/run_2_sample",
                              double range_diff = -1)
{
  TFile* f_ref = TFile::Open(ref_file, "READ");
  TFile* f_new = TFile::Open(new_file, "READ");
  if (!f_ref || f_ref->IsZombie()) { std::cerr << "ERROR: cannot open " << ref_file << std::endl; return; }
  if (!f_new || f_new->IsZombie()) { std::cerr << "ERROR: cannot open " << new_file << std::endl; return; }

  TString tag_ref = pdiff_tag_from_file(ref_file);
  TString tag_new = pdiff_tag_from_file(new_file);
  TString lab_ref = TString(label_ref).Length() ? TString(label_ref) : pdiff_label_from_tag(tag_ref);
  TString lab_new = TString(label_new).Length() ? TString(label_new) : pdiff_label_from_tag(tag_new);

  TString stem    = Form("czy_pdiff_%s_vs_%s", tag_new.Data(), tag_ref.Data());
  TString out_png = Form("%s/%s", out_dir, stem.Data());
  TString out_root = Form("%s/%s.root", out_dir, stem.Data());
  gSystem->mkdir(out_png, true);

  // 1) Collect percent differences for all 6 maps.
  std::vector<double> vals[3][2];
  std::vector<double> dvals[3][2];
  int skipped[3][2] = {{0}};
  int total[3][2]   = {{0}};
  std::vector<double> all_abs;

  for (int l = 0; l < 3; l++) for (int k = 0; k < 2; k++) {
    TString name = Form("CzyHist_%i_%i", l, k);
    TH2* h_ref = (TH2*)f_ref->Get(name);
    TH2* h_new = (TH2*)f_new->Get(name);
    if (!h_ref || !h_new) { std::cout << "Skipping missing " << name << std::endl; continue; }
    if (h_ref->GetNbinsX() != h_new->GetNbinsX() || h_ref->GetNbinsY() != h_new->GetNbinsY()) {
      std::cerr << "Binning mismatch for " << name << ", skipping" << std::endl; continue;
    }
    total[l][k] = h_ref->GetNbinsX() * h_ref->GetNbinsY();
    vals[l][k] = pdiff_values(h_ref, h_new, skipped[l][k]);
    dvals[l][k] = diff_values(h_ref, h_new);
    for (double x : vals[l][k]) all_abs.push_back(std::fabs(x));
  }

  if (all_abs.empty()) { std::cerr << "ERROR: no overlapping filled pixels found." << std::endl; return; }

  // 2) Binning: same for every histogram so they can be compared directly.
  //    Automatic: covers 99% of all |% diff| values, rounded up to a nice number.
  if (range <= 0) {
    std::sort(all_abs.begin(), all_abs.end());
    double q99 = all_abs[(size_t)std::floor(0.99 * (all_abs.size() - 1))];
    range = pdiff_nice_ceil(1.1 * q99);
  }
  std::cout << "x-axis range: -" << range << " % to +" << range << " %, " << nbins
            << " bins (" << 2 * range / nbins << " % per bin)" << std::endl;

  // 3) Fill, draw, save.
  gStyle->SetOptStat(0);
  TFile* f_out = new TFile(out_root, "RECREATE");
  TH1D* hists[3][2] = {{nullptr}};
  TString xtitle = Form("(C_{zy}[%s] - C_{zy}[%s]) / C_{zy}[%s] [%%]", lab_new.Data(), lab_ref.Data(), lab_ref.Data());

  printf("\n%-24s %-9s %8s %8s %10s %10s %8s %8s\n", "Plane", "TPC", "pixels", "skipped", "mean[%]", "RMS[%]", "under", "over");
  for (int l = 0; l < 3; l++) for (int k = 0; k < 2; k++) {
    if (total[l][k] == 0) continue;
    TH1D* h = new TH1D(Form("pdiff_CzyHist_%i_%i", l, k),
                       Form("%s, %s;%s;Number of YZ pixels", pdiff_plane_name(l), pdiff_tpc_name(k), xtitle.Data()),
                       nbins, -range, range);
    for (double x : vals[l][k]) h->Fill(x);
    h->SetLineWidth(2);
    h->SetLineColor(kAzure + 2);
    h->SetFillColorAlpha(kAzure + 2, 0.25);
    f_out->cd();
    h->Write();
    hists[l][k] = h;

    // mean/RMS from the raw values, so out-of-range pixels are included
    double mean = TMath::Mean(vals[l][k].begin(), vals[l][k].end());
    double rms  = TMath::RMS(vals[l][k].begin(), vals[l][k].end());
    printf("%-24s %-9s %8zu %8d %10.4f %10.4f %8.0f %8.0f\n", pdiff_plane_name(l), pdiff_tpc_name(k),
           vals[l][k].size(), skipped[l][k], mean, rms, h->GetBinContent(0), h->GetBinContent(nbins + 1));

    TCanvas c("c_pdiff", "", 1100, 850);
    c.SetLeftMargin(0.13); c.SetBottomMargin(0.13); c.SetRightMargin(0.05);
    h->GetXaxis()->SetTitleSize(0.042); h->GetYaxis()->SetTitleSize(0.045);
    h->SetMaximum(1.45 * h->GetMaximum());  // headroom for the text (same as the diff plots)
    h->Draw("HIST");
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
    t.DrawLatex(0.16, 0.85, Form("%s vs %s", lab_new.Data(), lab_ref.Data()));
    t.DrawLatex(0.16, 0.80, Form("Pixels: %zu (%d empty skipped)", vals[l][k].size(), skipped[l][k]));
    t.DrawLatex(0.16, 0.75, Form("Mean = %.3g %%,  RMS = %.3g %%", mean, rms));
    t.DrawLatex(0.16, 0.70, Form("Under/overflow: %.0f / %.0f", h->GetBinContent(0), h->GetBinContent(nbins + 1)));
    c.SaveAs(Form("%s/pdiff_CzyHist_%i_%i.png", out_png.Data(), l, k));
  }

  // Summary: 3 rows (planes) x 2 columns (East, West).
  TCanvas cs("c_pdiff_summary", "", 1400, 1600);
  cs.Divide(2, 3);
  for (int l = 0; l < 3; l++) for (int k = 0; k < 2; k++) {
    if (!hists[l][k]) continue;
    cs.cd(l * 2 + k + 1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.15);
    hists[l][k]->Draw("HIST");
  }
  cs.SaveAs(Form("%s/pdiff_CzyHist_summary.png", out_png.Data()));

  // 4) Same again for the signed difference Czy_new - Czy_ref.
  if (range_diff <= 0) {
    std::vector<double> all_abs_diff;
    for (int l = 0; l < 3; l++) for (int k = 0; k < 2; k++)
      for (double x : dvals[l][k]) all_abs_diff.push_back(std::fabs(x));
    std::sort(all_abs_diff.begin(), all_abs_diff.end());
    double q99 = all_abs_diff[(size_t)std::floor(0.99 * (all_abs_diff.size() - 1))];
    range_diff = pdiff_nice_ceil(1.1 * q99);
  }
  std::cout << "\nDifference x-axis range: -" << range_diff << " to +" << range_diff << ", " << nbins
            << " bins (" << 2 * range_diff / nbins << " per bin)" << std::endl;

  TH1D* dhists[3][2] = {{nullptr}};
  TString dxtitle = Form("C_{zy}[%s] - C_{zy}[%s]", lab_new.Data(), lab_ref.Data());

  printf("\n%-24s %-9s %8s %8s %12s %12s %8s %8s\n", "Plane", "TPC", "pixels", "skipped", "mean", "RMS", "under", "over");
  for (int l = 0; l < 3; l++) for (int k = 0; k < 2; k++) {
    if (total[l][k] == 0) continue;
    TH1D* h = new TH1D(Form("diff_CzyHist_%i_%i", l, k),
                       Form("%s, %s;%s;Number of YZ pixels", pdiff_plane_name(l), pdiff_tpc_name(k), dxtitle.Data()),
                       nbins, -range_diff, range_diff);
    for (double x : dvals[l][k]) h->Fill(x);
    h->SetLineWidth(2);
    h->SetLineColor(kAzure + 2);
    h->SetFillColorAlpha(kAzure + 2, 0.25);
    f_out->cd();
    h->Write();
    dhists[l][k] = h;

    double mean = TMath::Mean(dvals[l][k].begin(), dvals[l][k].end());
    double rms  = TMath::RMS(dvals[l][k].begin(), dvals[l][k].end());
    printf("%-24s %-9s %8zu %8d %12.3e %12.3e %8.0f %8.0f\n", pdiff_plane_name(l), pdiff_tpc_name(k),
           dvals[l][k].size(), skipped[l][k], mean, rms, h->GetBinContent(0), h->GetBinContent(nbins + 1));

    TCanvas c("c_diff", "", 1100, 850);
    c.SetLeftMargin(0.13); c.SetBottomMargin(0.13); c.SetRightMargin(0.05);
    h->GetXaxis()->SetTitleSize(0.042); h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetMaxDigits(5); h->GetXaxis()->SetNdivisions(505);  // plain decimals, no x10^-3
    h->SetMaximum(1.45 * h->GetMaximum());  // headroom for the text
    h->Draw("HIST");
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032);
    t.DrawLatex(0.16, 0.85, Form("%s vs %s", lab_new.Data(), lab_ref.Data()));
    t.DrawLatex(0.16, 0.80, Form("Pixels: %zu (%d empty skipped)", dvals[l][k].size(), skipped[l][k]));
    t.DrawLatex(0.16, 0.75, Form("Mean = %.3g,  RMS = %.3g", mean, rms));
    t.DrawLatex(0.16, 0.70, Form("Under/overflow: %.0f / %.0f", h->GetBinContent(0), h->GetBinContent(nbins + 1)));
    c.SaveAs(Form("%s/diff_CzyHist_%i_%i.png", out_png.Data(), l, k));
  }

  TCanvas cd("c_diff_summary", "", 1400, 1600);
  cd.Divide(2, 3);
  for (int l = 0; l < 3; l++) for (int k = 0; k < 2; k++) {
    if (!dhists[l][k]) continue;
    cd.cd(l * 2 + k + 1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.15);
    dhists[l][k]->Draw("HIST");
  }
  cd.SaveAs(Form("%s/diff_CzyHist_summary.png", out_png.Data()));

  f_out->Close();
  f_ref->Close();
  f_new->Close();
  std::cout << "\nWrote " << out_root << "\n  and PNGs in " << out_png << std::endl;
}
