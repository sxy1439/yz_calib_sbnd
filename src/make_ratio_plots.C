#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <iostream>

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

void SetSimpleStyle()
{
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.19);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.08);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.1, "Y");
  gStyle->SetTitleOffset(1.17, "Z");

  gStyle->SetNumberContours(255);
  gROOT->ForceStyle();
}

void SetHistRatioStyle(TH2* h,
                       const char* ztitle)
{
  h->SetTitle("");

  h->GetXaxis()->SetTitle("Z Coordinate [cm]");
  h->GetYaxis()->SetTitle("Y Coordinate [cm]");
  h->GetZaxis()->SetTitle(ztitle);

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();

  h->GetXaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetZaxis()->SetTitleSize(0.055);

  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);

  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetZaxis()->SetTitleOffset(1.17);

  h->SetStats(0);
}

const char* GetTPCLabel(const TString& name)
{
  if (name.Contains("_0")) return "East TPC";
  if (name.Contains("_1")) return "West TPC";
  return "";
}

bool keep_ratio_hist(const TString& name)
{
  // Keep only the histogram types you care about.
  // Start with correction factor maps.
  if (name.BeginsWith("CzyHist_")) return true;

  // Uncomment these if you also want ratio maps for median dQ/dx and hit counts.
  // if (name.BeginsWith("zy_")) return true;
  // if (name.BeginsWith("zynhits_")) return true;

  return false;
}

TH2D* MakeSafeRatio(TH2* h_num, TH2* h_den, const TString& ratio_name)
{
  if (!h_num || !h_den) return nullptr;

  if (h_num->GetNbinsX() != h_den->GetNbinsX() ||
      h_num->GetNbinsY() != h_den->GetNbinsY()) {
    std::cerr << "Binning mismatch for " << ratio_name << std::endl;
    return nullptr;
  }

  TH2D* h_ratio = (TH2D*)h_num->Clone(ratio_name);
  h_ratio->Reset();

  for (int ix = 1; ix <= h_num->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h_num->GetNbinsY(); iy++) {
      double num = h_num->GetBinContent(ix, iy);
      double den = h_den->GetBinContent(ix, iy);

      // Avoid meaningless divisions.
      // If denominator is zero or empty, leave ratio bin empty.
      if (den == 0) continue;

      h_ratio->SetBinContent(ix, iy, num / den);
    }
  }

  return h_ratio;
}

void save_ratio_png(TH2* h,
                    const TString& outdir,
                    const TString& name)
{
  SetSimpleStyle();

  TCanvas* c = new TCanvas("c", "c", 1100, 850);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  c->SetRightMargin(0.19);
  c->SetTopMargin(0.08);

  TString ztitle = "DENTFORCE / original";

  if (name.BeginsWith("ratio_CzyHist_")) {
    ztitle = "Correction factor ratio";
    h->GetZaxis()->SetRangeUser(0.9, 1.1);
  }
  else if (name.BeginsWith("ratio_zy_")) {
    ztitle = "Median dQ/dx ratio";
    h->GetZaxis()->SetRangeUser(0.8, 1.2);
  }
  else if (name.BeginsWith("ratio_zynhits_")) {
    ztitle = "Hit count ratio";
    h->GetZaxis()->SetRangeUser(0.0, 1.2);
  }

  SetHistRatioStyle(h, ztitle);
  h->Draw("COLZ");

  DrawLabel(GetTPCLabel(name), 0.7, 0.15, 0.95, kBlack, 12);
  DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
  DrawLabel("DENTFORCE / Original", 0.65, 0.50, 0.95, kBlack, 22);

  TString pngname = Form("%s/%s.png", outdir.Data(), name.Data());
  c->SaveAs(pngname);

  delete c;
}

void make_ratio_plots(const char* original_file,
                      const char* dentforce_file,
                      const char* out_root_file = "../output_files/ratio_dentforce_over_original.root",
                      const char* out_png_dir = "../plots/ratio_dentforce_over_original")
{
  gSystem->mkdir(out_png_dir, true);

  TFile* f_orig = TFile::Open(original_file, "READ");
  TFile* f_dent = TFile::Open(dentforce_file, "READ");

  if (!f_orig || f_orig->IsZombie()) {
    std::cerr << "ERROR: Could not open original file: " << original_file << std::endl;
    return;
  }

  if (!f_dent || f_dent->IsZombie()) {
    std::cerr << "ERROR: Could not open DENTFORCE file: " << dentforce_file << std::endl;
    return;
  }

  TFile* f_out = new TFile(out_root_file, "RECREATE");

  TIter next(f_orig->GetListOfKeys());
  TKey* key;

  int nmade = 0;

  while ((key = (TKey*)next())) {
    TString name = key->GetName();

    if (!keep_ratio_hist(name)) continue;

    TObject* obj_orig = key->ReadObj();

    if (!obj_orig->InheritsFrom(TH2::Class())) {
      delete obj_orig;
      continue;
    }

    TH2* h_orig = (TH2*)obj_orig;
    TH2* h_dent = (TH2*)f_dent->Get(name);

    if (!h_dent) {
      std::cout << "Skipping missing DENTFORCE hist: " << name << std::endl;
      delete obj_orig;
      continue;
    }

    TString ratio_name = "ratio_" + name;

    TH2D* h_ratio = MakeSafeRatio(h_dent, h_orig, ratio_name);

    if (!h_ratio) {
      delete obj_orig;
      continue;
    }

    f_out->cd();
    h_ratio->Write();

    save_ratio_png(h_ratio, out_png_dir, ratio_name);

    std::cout << "Made ratio: " << ratio_name << std::endl;

    delete h_ratio;
    delete obj_orig;

    nmade++;
  }

  f_out->Close();
  f_orig->Close();
  f_dent->Close();

  std::cout << "Done. Made " << nmade << " ratio plots." << std::endl;
}