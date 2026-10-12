// VERSION 1
// overlay_dqdx_lifetime.C
// //
// // Overlay dqdxHist_<plane>_<tpc> histograms from multiple ROOT files.
// // Intended for comparing yz_median_planes.C outputs produced with different
// // electron lifetime assumptions.
// //
// // Produces one PNG per plane/TPC:
// //   dqdxHist_0_0_overlay.png
// //   dqdxHist_0_1_overlay.png
// //   ...
// //   dqdxHist_2_1_overlay.png
// //
// // Set normalize=true to compare shapes.
// // Set normalize=false to compare raw counts.

// #include <TFile.h>
// #include <TH1.h>
// #include <TCanvas.h>
// #include <TLegend.h>
// #include <TPaveText.h>
// #include <TString.h>
// #include <TSystem.h>
// #include <TStyle.h>
// #include <TROOT.h>

// #include <iostream>
// #include <vector>
// #include <algorithm>


// static TString sanitize(TString s)
// {
//   s.ReplaceAll("/", "__");
//   s.ReplaceAll(" ", "_");
//   s.ReplaceAll("(", "");
//   s.ReplaceAll(")", "");
//   s.ReplaceAll(",", "_");
//   return s;
// }

// const char* get_plane_label(const TString& hname)
// {
//   if (hname.Contains("dqdxHist_0_")) return "U plane";
//   if (hname.Contains("dqdxHist_1_")) return "V plane";
//   if (hname.Contains("dqdxHist_2_")) return "Y plane";
//   return "Unknown plane";
// }

// const char* get_tpc_label(const TString& hname)
// {
//   if (hname.EndsWith("_0")) return "East TPC";
//   if (hname.EndsWith("_1")) return "West TPC";
//   return "Unknown TPC";
// }

// bool is_east_tpc(const TString& hname)
// {
//   return hname.EndsWith("_0");
// }



// double hist_median(TH1* h)
// {
//   if (!h || h->GetEntries() <= 0) return 0.0;

//   double q = 0.5;
//   double median = 0.0;
//   h->GetQuantiles(1, &median, &q);

//   return median;
// }


// void format_hist_1d(TH1* h)
// {
//   if (!h) return;

//   h->SetStats(0);
//   h->SetLineWidth(2);

//   h->GetXaxis()->SetTitleSize(0.045);
//   h->GetYaxis()->SetTitleSize(0.045);

//   h->GetXaxis()->SetLabelSize(0.04);
//   h->GetYaxis()->SetLabelSize(0.04);

//   h->GetXaxis()->SetTitleOffset(1.10);
//   h->GetYaxis()->SetTitleOffset(1.30);
// }


// void overlay_one_dqdx_hist(const std::vector<TFile*>& files,
//                            const std::vector<TString>& labels,
//                            const TString& histname,
//                            const TString& outname,
//                            bool normalize = true)
// {
//   if (files.size() != labels.size()) {
//     std::cerr << "ERROR: files and labels have different sizes.\n";
//     return;
//   }

//   std::vector<TH1*> hists;
//   std::vector<double> medians;
//   std::vector<double> means;
//   std::vector<double> stddevs;
//   std::vector<double> entries;

//   // Colors for up to several lifetime variations.
//   std::vector<int> colors = {
//     kBlack,
//     kRed + 1,
//     kBlue + 1,
//     kGreen + 2,
//     kMagenta + 1,
//     kOrange + 7,
//     kCyan + 2,
//     kViolet + 1
//   };

//   for (size_t i = 0; i < files.size(); ++i) {
//     if (!files[i]) continue;

//     TH1* h = (TH1*)files[i]->Get(histname);
//     if (!h) {
//       std::cerr << "WARNING: could not find " << histname
//                 << " in file " << files[i]->GetName() << "\n";
//       continue;
//     }

//     TH1* hc = (TH1*)h->Clone(Form("%s_clone_%zu", histname.Data(), i));
//     hc->SetDirectory(nullptr);

//     // Store stats before optional normalization.
//     entries.push_back(hc->GetEntries());
//     means.push_back(hc->GetMean());
//     stddevs.push_back(hc->GetStdDev());
//     medians.push_back(hist_median(hc));

//     if (normalize) {
//       double integral = hc->Integral();
//       if (integral > 0) hc->Scale(1.0 / integral);
//     }

//     format_hist_1d(hc);

//     int color = colors[i % colors.size()];
//     hc->SetLineColor(color);
//     hc->SetMarkerColor(color);

//     hists.push_back(hc);
//   }

//   if (hists.empty()) {
//     std::cerr << "WARNING: no histograms found for " << histname << "\n";
//     return;
//   }

//   TCanvas* c = new TCanvas("c_overlay_dqdx", "c_overlay_dqdx", 1400, 850);
//   c->cd();

//   gStyle->SetOptStat(0);

//   gPad->SetLeftMargin(0.12);
//   gPad->SetBottomMargin(0.12);
//   gPad->SetRightMargin(0.34);
//   gPad->SetTopMargin(0.08);

//   double ymax = 0.0;
//   for (auto* h : hists) {
//     ymax = std::max(ymax, h->GetMaximum());
//   }

//   hists[0]->SetMaximum(1.18 * ymax);
//   // hists[0]->SetTitle("dqdx"); I replaced "dqdx" with plot_title
//   TString plot_title = Form("dqdx - %s, %s",
//                           get_plane_label(histname),
//                           get_tpc_label(histname));
//   hists[0]->SetTitle(plot_title);

//   hists[0]->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");

//   if (normalize) {
//     hists[0]->GetYaxis()->SetTitle("Fraction of hits");
//   }
//   else {
//     hists[0]->GetYaxis()->SetTitle("Number of hits");
//   }

//   hists[0]->Draw("HIST");

//   for (size_t i = 1; i < hists.size(); ++i) {
//     hists[i]->Draw("HIST SAME");
//   }

//   TLegend* leg = new TLegend(0.69, 0.62, 0.98, 0.90);
//   leg->SetBorderSize(1);
//   leg->SetFillStyle(1001);

//   for (size_t i = 0; i < hists.size(); ++i) {
//     leg->AddEntry(hists[i], labels[i], "l");
//   }

//   leg->Draw();

//   // Compact stats box.
//   // Shows median mainly because electron lifetime changes should shift dQ/dx.
//   TPaveText* stats = new TPaveText(0.69, 0.22, 0.98, 0.58, "NDC");

//   stats->SetBorderSize(1);
//   stats->SetFillStyle(1001);
//   stats->SetTextAlign(12);
//   stats->SetTextSize(0.022);

//   for (size_t i = 0; i < hists.size(); ++i) {
//     stats->AddText(Form("%s: med=%.1f, #sigma=%.1f",
//                         labels[i].Data(),
//                         medians[i],
//                         stddevs[i]));
//   }

//   stats->Draw();

//   c->Modified();
//   c->Update();
//   c->SaveAs(outname);

//   delete stats;
//   delete leg;

//   for (auto* h : hists) {
//     delete h;
//   }

//   delete c;
// }


// void overlay_dqdx_lifetime(const std::vector<TString>& rootfiles,
//                            const std::vector<TString>& labels,
//                            const char* outdir = "dqdx_lifetime_overlay",
//                            bool normalize = true)
// {
//   if (rootfiles.size() != labels.size()) {
//     std::cerr << "ERROR: rootfiles and labels must have the same size.\n";
//     return;
//   }

//   gSystem->mkdir(outdir, true);

//   std::vector<TFile*> files;

//   for (const auto& fname : rootfiles) {
//     TFile* f = TFile::Open(fname, "READ");

//     if (!f || f->IsZombie()) {
//       std::cerr << "ERROR: could not open " << fname << "\n";
//       if (f) delete f;
//       continue;
//     }

//     files.push_back(f);
//   }

//   if (files.empty()) {
//     std::cerr << "ERROR: no valid ROOT files opened.\n";
//     return;
//   }

//   for (int plane = 0; plane < 3; ++plane) {
//     for (int tpc = 0; tpc < 2; ++tpc) {

//       TString histname = Form("dqdxHist_%d_%d", plane, tpc);

//       TString mode = normalize ? "shape" : "counts";

//       TString outname = TString::Format(
//         "%s/%s_%s_overlay.png",
//         outdir,
//         sanitize(histname).Data(),
//         mode.Data()
//       );

//       overlay_one_dqdx_hist(files, labels, histname, outname, normalize);
//     }
//   }

//   for (auto* f : files) {
//     f->Close();
//     delete f;
//   }

//   std::cout << "Saved dqdxHist overlays to " << outdir << "\n";
// }























// VERSION 2. Difference with version 1 is 
// overlay_dqdx_lifetime.C
//
// Overlay dqdxHist_<plane>_<tpc> histograms from multiple ROOT files.
// Intended for comparing yz_median_planes.C / yz_median_planes_sce.C outputs
// produced with different electron lifetime assumptions.
//
// Histogram naming convention assumed:
//   dqdxHist_<plane>_<tpc>
//
// Plane:
//   0 = U plane
//   1 = V plane
//   2 = Y plane
//
// TPC:
//   0 = East TPC
//   1 = West TPC
//
// For *_0 histograms, the legend uses east_lifetimes.
// For *_1 histograms, the legend uses west_lifetimes.

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <algorithm>


static TString sanitize(TString s)
{
  s.ReplaceAll("/", "__");
  s.ReplaceAll(" ", "_");
  s.ReplaceAll("(", "");
  s.ReplaceAll(")", "");
  s.ReplaceAll(",", "_");
  return s;
}


double hist_median(TH1* h)
{
  if (!h || h->GetEntries() <= 0) return 0.0;

  double q = 0.5;
  double median = 0.0;
  h->GetQuantiles(1, &median, &q);

  return median;
}


const char* get_plane_label(const TString& hname)
{
  if (hname.Contains("dqdxHist_0_")) return "U plane";
  if (hname.Contains("dqdxHist_1_")) return "V plane";
  if (hname.Contains("dqdxHist_2_")) return "Y plane";

  return "Unknown plane";
}


const char* get_tpc_label(const TString& hname)
{
  if (hname.EndsWith("_0")) return "East TPC";
  if (hname.EndsWith("_1")) return "West TPC";

  return "Unknown TPC";
}


bool is_east_hist(const TString& hname)
{
  return hname.EndsWith("_0");
}


TString lifetime_legend_label(const TString& hname,
                              const std::vector<double>& east_lifetimes,
                              const std::vector<double>& west_lifetimes,
                              size_t i)
{
  double tau = is_east_hist(hname) ? east_lifetimes[i] : west_lifetimes[i];

  return Form("#tau_{e} = %.0f ms", tau);
}


void format_hist_1d(TH1* h)
{
  if (!h) return;

  h->SetStats(0);
  h->SetLineWidth(2);

  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);

  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);

  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.30);
}


void overlay_one_dqdx_hist(const std::vector<TFile*>& files,
                           const std::vector<double>& east_lifetimes,
                           const std::vector<double>& west_lifetimes,
                           const TString& histname,
                           const TString& outname,
                           bool normalize = true)
{
  if (files.size() != east_lifetimes.size() ||
      files.size() != west_lifetimes.size()) {
    std::cerr << "ERROR: files, east_lifetimes, and west_lifetimes must have the same size.\n";
    return;
  }

  std::vector<TH1*> hists;
  std::vector<double> medians;
  std::vector<double> stddevs;
  std::vector<double> entries;

  std::vector<int> colors = {
    kBlack,
    kRed + 1,
    kBlue + 1,
    kGreen + 2,
    kMagenta + 1,
    kOrange + 7,
    kCyan + 2,
    kViolet + 1
  };

  for (size_t i = 0; i < files.size(); ++i) {
    if (!files[i]) continue;

    TH1* h = (TH1*)files[i]->Get(histname);

    if (!h) {
      std::cerr << "WARNING: could not find " << histname
                << " in file " << files[i]->GetName() << "\n";
      continue;
    }

    TH1* hc = (TH1*)h->Clone(Form("%s_clone_%zu", histname.Data(), i));
    hc->SetDirectory(nullptr);

    entries.push_back(hc->GetEntries());
    medians.push_back(hist_median(hc));
    stddevs.push_back(hc->GetStdDev());

    if (normalize) {
      double integral = hc->Integral();
      if (integral > 0) hc->Scale(1.0 / integral);
    }

    format_hist_1d(hc);

    int color = colors[i % colors.size()];
    hc->SetLineColor(color);
    hc->SetMarkerColor(color);

    hists.push_back(hc);
  }

  if (hists.empty()) {
    std::cerr << "WARNING: no histograms found for " << histname << "\n";
    return;
  }

  TCanvas* c = new TCanvas("c_overlay_dqdx", "c_overlay_dqdx", 1400, 850);
  c->cd();

  gStyle->SetOptStat(0);

  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.34);
  gPad->SetTopMargin(0.08);

  double ymax = 0.0;
  for (auto* h : hists) {
    ymax = std::max(ymax, h->GetMaximum());
  }

  hists[0]->SetMaximum(1.18 * ymax);

  TString plot_title = Form("dqdx - %s, %s",
                            get_plane_label(histname),
                            get_tpc_label(histname));

  hists[0]->SetTitle(plot_title);
  hists[0]->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");

  if (normalize) {
    hists[0]->GetYaxis()->SetTitle("Fraction of hits");
  }
  else {
    hists[0]->GetYaxis()->SetTitle("Number of hits");
  }

  hists[0]->Draw("HIST");

  for (size_t i = 1; i < hists.size(); ++i) {
    hists[i]->Draw("HIST SAME");
  }

  TLegend* leg = new TLegend(0.69, 0.62, 0.98, 0.90);
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.030);

  for (size_t i = 0; i < hists.size(); ++i) {
    TString leg_label = lifetime_legend_label(histname,
                                              east_lifetimes,
                                              west_lifetimes,
                                              i);

    leg->AddEntry(hists[i], leg_label, "l");
  }

  leg->Draw();

  TPaveText* stats = new TPaveText(0.69, 0.22, 0.98, 0.58, "NDC");
  stats->SetBorderSize(1);
  stats->SetFillStyle(1001);
  stats->SetTextAlign(12);
  stats->SetTextSize(0.030);

  for (size_t i = 0; i < hists.size(); ++i) {
    TString stat_label = lifetime_legend_label(histname,
                                               east_lifetimes,
                                               west_lifetimes,
                                               i);

    stats->AddText(Form("%s:", stat_label.Data()));
    stats->AddText(Form("  med=%.1f, #sigma=%.1f",
                        medians[i],
                        stddevs[i]));
  }

  stats->Draw();

  c->Modified();
  c->Update();
  c->SaveAs(outname);

  delete stats;
  delete leg;

  for (auto* h : hists) {
    delete h;
  }

  delete c;
}


void overlay_dqdx_lifetime(const std::vector<TString>& rootfiles,
                           const std::vector<double>& east_lifetimes,
                           const std::vector<double>& west_lifetimes,
                           const char* outdir = "dqdx_lifetime_overlay",
                           bool normalize = true)
{
  if (rootfiles.size() != east_lifetimes.size() ||
      rootfiles.size() != west_lifetimes.size()) {
    std::cerr << "ERROR: rootfiles, east_lifetimes, and west_lifetimes must have the same size.\n";
    return;
  }

  gSystem->mkdir(outdir, true);

  std::vector<TFile*> files;

  for (const auto& fname : rootfiles) {
    TFile* f = TFile::Open(fname, "READ");

    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: could not open " << fname << "\n";
      if (f) delete f;
      continue;
    }

    files.push_back(f);
  }

  if (files.empty()) {
    std::cerr << "ERROR: no valid ROOT files opened.\n";
    return;
  }

  for (int plane = 0; plane < 3; ++plane) {
    for (int tpc = 0; tpc < 2; ++tpc) {

      TString histname = Form("dqdxHist_%d_%d", plane, tpc);

      TString mode = normalize ? "shape" : "counts";

      TString outname = TString::Format(
        "%s/%s_%s_overlay.png",
        outdir,
        sanitize(histname).Data(),
        mode.Data()
      );

      overlay_one_dqdx_hist(files,
                            east_lifetimes,
                            west_lifetimes,
                            histname,
                            outname,
                            normalize);
    }
  }

  for (auto* f : files) {
    f->Close();
    delete f;
  }

  std::cout << "Saved dqdxHist overlays to " << outdir << "\n";
}