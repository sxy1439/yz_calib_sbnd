#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <iostream>

#define SBNDSTYLE_ENABLE_AUTOMATICALLY 0
#include "../plots/SBNDStyle.h"


/*
This macro is used to make the median dQ/dx YZ maps in x-slices by running in the relevant ROOT file
created by yz_median_planes_sce_xslices.C to produce PNGs.
It is based on yz_median_planes_sce.C, but with additional histograms
for each x-slice. The number of x-slices and their boundaries can be
adjusted in the GetXSectionFullRange() function.
*/

const int nxsections = 8;
const double xslice_low = -20.0;
const double xslice_width = 20.0;

void DrawLabel(const char* text="text",
               float tscale=0.8,
               double xloc=0.9,
               double yloc=0.95,
               int tcolor=kBlack,
               int align=12)
{
  TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  label->SetTextColor(tcolor);
  label->SetTextFont(42);
  label->SetNDC();
  label->SetTextSize(2.0/30.0);
  label->SetTextAlign(align);
  label->Draw();
}

void SetHistSBND(TH2F *h,
                 const char* xlabel,
                 const char* ylabel,
                 const char* zlabel)
{
  h->SetTitle("");

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();

  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  h->GetZaxis()->SetTitle(zlabel);

  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetZaxis()->SetTitleOffset(1.17);

  h->GetXaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetZaxis()->SetTitleSize(0.055);

  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);
}

void make_xslice_gif_frames(const char* infile,
                            const char* outdir = "../plots/xslice_gif_frames",
                            int plane = 2,
                            const char* hist_prefix = "CzyHist")
{
  gSystem->mkdir(outdir, true);

  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Could not open " << infile << std::endl;
    return;
  }

  sbndstyle::SetSBNDStyle();
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  sbndstyle::SeaPalette();

  int frame = 0;

  for (int xs = 0; xs < nxsections; xs++) {
    int tpc = (xs < 10) ? 0 : 1;

    TString hname = Form("%s_%i_%i_xsec%i", hist_prefix, plane, tpc, xs);
    TH2F* h = (TH2F*)f->Get(hname);

    if (!h) {
      std::cout << "Skipping missing histogram: " << hname << std::endl;
      continue;
    }

    if (h->GetEntries() == 0) {
      std::cout << "Skipping empty histogram: " << hname << std::endl;
      continue;
    }

    double x1 = xslice_low + xs * xslice_width;
    double x2 = x1 + xslice_width;

    TCanvas* c = new TCanvas("c", "c", 1100, 850);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.19);
    c->SetTopMargin(0.08);

    if (TString(hist_prefix) == "CzyHist") {
      SetHistSBND(h,
                  "Z Coordinate [cm]",
                  "Y Coordinate [cm]",
                  "YZ correction factor");

      h->GetZaxis()->SetRangeUser(0.8, 1.12);
    }
    else if (TString(hist_prefix) == "zy") {
      SetHistSBND(h,
                  "Z Coordinate [cm]",
                  "Y Coordinate [cm]",
                  "Median dQ/dx [ADC/cm]");

      h->GetZaxis()->SetRangeUser(500, 1300);
    }
    else if (TString(hist_prefix) == "zynhits") {
      SetHistSBND(h,
                  "Z Coordinate [cm]",
                  "Y Coordinate [cm]",
                  "number of hits");
    }

    h->SetStats(0);
    h->Draw("COLZ");

    const char* tpc_label = (tpc == 0) ? "East TPC" : "West TPC";

    DrawLabel(tpc_label, 0.7, 0.15, 0.95, kBlack, 12);
    DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);

    DrawLabel(Form("%.0f < X < %.0f cm", x1, x2),
              0.7, 0.50, 0.95, kBlack, 22);

    TString outname = Form("%s/frame_%03d.png", outdir, frame);
    c->SaveAs(outname);

    std::cout << "Saved " << outname << " from " << hname << std::endl;

    delete c;
    frame++;
  }

  f->Close();

  std::cout << "Saved " << frame << " frames to " << outdir << std::endl;
}