#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>

#include <iostream>
#include <algorithm>

const int nplanes = 3;
const int nxsections = 20;

const double xslice_low = -200.0;
const double xslice_width = 20.0;

const char* plane_label(int plane)
{
  if (plane == 0) return "U";
  if (plane == 1) return "V";
  if (plane == 2) return "Y";
  return "Unknown";
}

const char* tpc_label(int tpc)
{
  if (tpc == 0) return "East";
  if (tpc == 1) return "West";
  return "Unknown";
}

void NormaliseHistogram(TH1* h)
{
  if (!h) return;

  double integral = h->Integral();

  if (integral > 0.0) {
    h->Scale(1.0 / integral);
  }
}

void overlay_dqdx_xslice_original_vs_dentforce(
    const char* original_file,
    const char* dentforce_file,
    const char* out_dir = "../plots/dqdx_xslice_original_vs_dentforce",
    bool only_physical_tpc_xslices = true)
{
  gStyle->SetOptStat(0);

  gSystem->mkdir(out_dir, true);

  TFile* f_orig = TFile::Open(original_file, "READ");
  TFile* f_dent = TFile::Open(dentforce_file, "READ");

  if (!f_orig || f_orig->IsZombie()) {
    std::cerr << "ERROR: Could not open original file: "
              << original_file << std::endl;
    return;
  }

  if (!f_dent || f_dent->IsZombie()) {
    std::cerr << "ERROR: Could not open DENTFORCE file: "
              << dentforce_file << std::endl;
    return;
  }

  for (int plane = 0; plane < nplanes; plane++) {
    for (int tpc = 0; tpc < 2; tpc++) {
      for (int xs = 0; xs < nxsections; xs++) {

        // For full -200 to 200 cm slicing:
        // xsec0-9 are negative x / East TPC
        // xsec10-19 are positive x / West TPC
        if (only_physical_tpc_xslices) {
          if (tpc == 0 && xs >= 10) continue;
          if (tpc == 1 && xs < 10) continue;
        }

        double xlow = xslice_low + xs * xslice_width;
        double xhigh = xlow + xslice_width;

        TString hname = Form("dqdxHist_%d_%d_xsec%d",
                             plane, tpc, xs);

        TH1* h_orig_raw = (TH1*)f_orig->Get(hname);
        TH1* h_dent_raw = (TH1*)f_dent->Get(hname);

        if (!h_orig_raw || !h_dent_raw) {
          std::cout << "Skipping missing histogram: "
                    << hname << std::endl;
          continue;
        }

        if (h_orig_raw->GetEntries() <= 0 ||
            h_dent_raw->GetEntries() <= 0) {
          std::cout << "Skipping empty histogram: "
                    << hname << std::endl;
          continue;
        }

        TH1* h_orig = (TH1*)h_orig_raw->Clone(Form("%s_original_clone",
                                                   hname.Data()));
        TH1* h_dent = (TH1*)h_dent_raw->Clone(Form("%s_dentforce_clone",
                                                   hname.Data()));

        h_orig->SetDirectory(0);
        h_dent->SetDirectory(0);

        NormaliseHistogram(h_orig);
        NormaliseHistogram(h_dent);

        h_orig->SetLineColor(kBlue + 1);
        h_orig->SetLineWidth(2);

        h_dent->SetLineColor(kRed + 1);
        h_dent->SetLineWidth(2);

        h_orig->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");
        h_orig->GetYaxis()->SetTitle("Fraction of hits");

        TString title = Form("%s plane, %s TPC, xsec%d: %.0f < X < %.0f cm",
                             plane_label(plane),
                             tpc_label(tpc),
                             xs,
                             xlow,
                             xhigh);

        h_orig->SetTitle(title);

        double ymax = std::max(h_orig->GetMaximum(),
                               h_dent->GetMaximum());

        h_orig->SetMaximum(1.25 * ymax);

        TCanvas* c = new TCanvas(Form("c_%s", hname.Data()),
                                 Form("c_%s", hname.Data()),
                                 900, 700);

        h_orig->Draw("HIST");
        h_dent->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.58, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        leg->AddEntry(h_orig,
                      Form("Original, entries = %.0f",
                           h_orig_raw->GetEntries()),
                      "l");

        leg->AddEntry(h_dent,
                      Form("DENTFORCE, entries = %.0f",
                           h_dent_raw->GetEntries()),
                      "l");

        leg->Draw();

        TString out_png = Form("%s/dqdx_overlay_plane%d_tpc%d_xsec%d.png",
                               out_dir,
                               plane,
                               tpc,
                               xs);

        c->SaveAs(out_png);

        std::cout << "Saved: " << out_png << std::endl;

        delete c;
        delete h_orig;
        delete h_dent;
      }
    }
  }

  f_orig->Close();
  f_dent->Close();

  std::cout << "Done making overlaid, area-normalised dQ/dx plots."
            << std::endl;
}