// // Version 3 

// #include <TFile.h>
// #include <TKey.h>
// #include <TClass.h>
// #include <TSystem.h>
// #include <TCanvas.h>
// #include <TString.h>
// #include <TDirectory.h>
// #include <TObject.h>
// #include <TH1.h>
// #include <TH2.h>
// #include <TStyle.h>
// #include <TPaveStats.h>
// #include <TPaletteAxis.h>
// #include <TList.h>
// #include <iostream>
// #define SBNDSTYLE_ENABLE_AUTOMATICALLY 0
// #include "../plots/SBNDStyle.h"
// #include "TLatex.h"

// void SetHistSBND(TH2F *h,
//                  const char* title,
//                  const char* xlabel,
//                  const char* ylabel,
//                  const char* zlabel)
// {
//   h->SetTitle(title);

//   h->GetXaxis()->CenterTitle();
//   h->GetYaxis()->CenterTitle();
//   h->GetZaxis()->CenterTitle();

//   h->GetXaxis()->SetTitle(xlabel);
//   h->GetYaxis()->SetTitle(ylabel);
//   h->GetZaxis()->SetTitle(zlabel);

//   h->GetXaxis()->SetTitleOffset(1.1);
//   h->GetYaxis()->SetTitleOffset(1.1);
//   h->GetZaxis()->SetTitleOffset(1.17);

//   h->GetXaxis()->SetTitleSize(0.055);
//   h->GetYaxis()->SetTitleSize(0.055);
//   h->GetZaxis()->SetTitleSize(0.055);

//   h->GetXaxis()->SetLabelSize(0.05);
//   h->GetYaxis()->SetLabelSize(0.05);
//   h->GetZaxis()->SetLabelSize(0.05);
// }

// void DrawLabel(const char* text="text",
//                float tscale=0.8,
//                double xloc=0.9,
//                double yloc=0.95,
//                int tcolor=kGray+1,
//                int align=12)
// {
//   TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
//   label->SetTextColor(tcolor);
//   label->SetTextFont(42);
//   label->SetNDC();
//   label->SetTextSize(2.0/30.0);
//   label->SetTextAlign(align);
//   label->Draw();
// }


// const char* get_tpc_label(const TString& name) {
//   if (name.Contains("_0_xsec")) return "East TPC";
//   if (name.Contains("_1_xsec")) return "West TPC";
//   return "";
// }

// bool keep_plot(const TString& name) {
//   static std::vector<TString> keep = {


//         // Plane 0 / U
//     "CzyHist_0_0",
//     "CzyHist_0_1",

//     "zy_0_0",
//     "zy_0_1",

//     "zynhits_0_0",
//     "zynhits_0_1",

//     "dqdxHist_0_0",
//     "dqdxHist_0_1",

//     // Plane 1 useful combinations
//     "CzyHist_1_0",
//     "CzyHist_1_1",

//     "zy_1_0",
//     "zy_1_1",

//     "zynhits_1_0",
//     "zynhits_1_1",

//     "dqdxHist_1_0",
//     "dqdxHist_1_1",

//     // Plane 2 useful combinations
//     "CzyHist_2_0",
//     "CzyHist_2_1",

//     "zy_2_0",
//     "zy_2_1",

//     "zynhits_2_0",
//     "zynhits_2_1",

//     "dqdxHist_2_0",
//     "dqdxHist_2_1",

//   };

//   for (const auto& k : keep) {
//     if (name == k) return true;
//   }

//   return false;
// }

// static TString sanitize(TString s) {
//   s.ReplaceAll("/", "__");
//   s.ReplaceAll(" ", "_");
//   return s;
// }


// void save_histogram(TObject* obj, const TString& outdir, int index) {
//   if (!obj) return;

//   TString name = obj->GetName();
//   TString out = TString::Format("%s/%03d_%s.png", outdir.Data(), index, sanitize(name).Data());

//   TCanvas* c = new TCanvas("c", "c", 1100, 850);
//   c->cd();

//   sbndstyle::SetSBNDStyle();
//   gROOT->ForceStyle();
//   gStyle->SetOptStat(0);
//   sbndstyle::SeaPalette();

//   c->SetLeftMargin(0.15);
//   c->SetBottomMargin(0.15);
//   c->SetRightMargin(0.19);
//   c->SetTopMargin(0.08);

//   if (obj->InheritsFrom(TH2::Class())) {
//   TH2F* h2 = (TH2F*)obj;
//   TString name = h2->GetName();

//   sbndstyle::SetSBNDStyle();
//   gROOT->ForceStyle();
//   gStyle->SetOptStat(0);
//   sbndstyle::SeaPalette();

//   c->SetLeftMargin(0.15);
//   c->SetBottomMargin(0.15);
//   c->SetRightMargin(0.19);
//   c->SetTopMargin(0.08);

//   // Decide labels and z range based on histogram type
//   if (name.BeginsWith("CzyHist")) {
//     SetHistSBND(h2, "",
//                 "Z Coordinate [cm]",
//                 "Y Coordinate [cm]",
//                 "YZ correction factor");

//     h2->GetZaxis()->SetRangeUser(0.8, 1.12);
//   }
//   else if (name.BeginsWith("zy_")) {
//     SetHistSBND(h2, "",
//                 "Z Coordinate [cm]",
//                 "Y Coordinate [cm]",
//                 "Median dQ/dx [ADC/cm]");

//     h2->GetZaxis()->SetRangeUser(500, 1300);
//   }
//   else if (name.BeginsWith("zynhits")) {
//     SetHistSBND(h2, "",
//                 "Z Coordinate [cm]",
//                 "Y Coordinate [cm]",
//                 "number of hits");
//   }
//   else {
//     SetHistSBND(h2, "",
//                 "Z Coordinate [cm]",
//                 "Y Coordinate [cm]",
//                 h2->GetZaxis()->GetTitle());
//   }

//   h2->SetStats(0);
//   h2->Draw("COLZ");

//   DrawLabel(get_tpc_label(name), 0.7, 0.15, 0.95, kBlack, 12);
//   DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);

//   c->Modified();
//   c->Update();
//   c->SaveAs(out);
// }


//   else if (obj->InheritsFrom(TH1::Class())) {
//     TH1* h1 = (TH1*)obj;

//     h1->SetTitle(name);
//     h1->GetXaxis()->SetTitleSize(0.045);
//     h1->GetYaxis()->SetTitleSize(0.045);

//     h1->GetXaxis()->SetLabelSize(0.04);
//     h1->GetYaxis()->SetLabelSize(0.04);

//     h1->GetXaxis()->SetTitleOffset(1.1);
//     h1->GetYaxis()->SetTitleOffset(1.3);

//     h1->Draw();
//     c->Update();

//     TPaveStats* st = (TPaveStats*)h1->FindObject("stats");
//     if (st) {
//       st->SetX1NDC(0.78);
//       st->SetX2NDC(0.98);
//       st->SetY1NDC(0.78);
//       st->SetY2NDC(0.92);
//       st->SetBorderSize(1);
//       st->SetFillStyle(1001);
//     }

//     c->Modified();
//     c->Update();
//     c->SaveAs(out);
//   }

//   delete c;
// }

// void walk_directory(TDirectory* dir, const TString& outdir, int& index) {
//   if (!dir) return;

//   TIter next(dir->GetListOfKeys());
//   TKey* key;

//   while ((key = (TKey*)next())) {
//     TObject* obj = key->ReadObj();
//     if (!obj) continue;

//     if (obj->InheritsFrom(TDirectory::Class())) {
//       walk_directory((TDirectory*)obj, outdir, index);
//       delete obj;
//       continue;
//     }



//     // if (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(TH2::Class())) {
//     //   save_histogram(obj, outdir, index);
//     //   ++index;
//     // }

//     if (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(TH2::Class())) {
//     TString name = obj->GetName();

//     if (!keep_plot(name)) {
//       delete obj;
//       continue;
//     }

//     save_histogram(obj, outdir, index);
//     ++index;
//   }

//     delete obj;
//   }
// }

// void export_hists_to_png(const char* rootfile,
//                                    const char* outdir = "png_out")
// {
//   gSystem->mkdir(outdir, true);

//   TFile f(rootfile, "READ");
//   if (f.IsZombie()) {
//     std::cerr << "ERROR: cannot open " << rootfile << "\n";
//     return;
//   }

//   gStyle->SetOptStat(0);

//   int index = 0;
//   walk_directory(&f, outdir, index);

//   f.Close();

//   std::cout << "Exported " << index << " histograms to " << outdir << "\n";
// }






































// Version 5
// Export selected histograms from a ROOT file to PNG.
// TH2 maps use SBND style.
// dqdxHist TH1 histograms use a clean ROOT-like style with title, axes, and stats box.

#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TString.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaletteAxis.h>
#include <TList.h>
#include <TROOT.h>
#include <TLatex.h>

#include <iostream>
#include <vector>

#define SBNDSTYLE_ENABLE_AUTOMATICALLY 0
#include "../plots/SBNDStyle.h"


void SetHistSBND(TH2F *h,
                 const char* title,
                 const char* xlabel,
                 const char* ylabel,
                 const char* zlabel)
{
  h->SetTitle(title);

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


void DrawLabel(const char* text = "text",
               float tscale = 0.8,
               double xloc = 0.9,
               double yloc = 0.95,
               int tcolor = kGray + 1,
               int align = 12)
{
  TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  label->SetTextColor(tcolor);
  label->SetTextFont(42);
  label->SetNDC();
  label->SetTextSize(2.0 / 30.0);
  label->SetTextAlign(align);
  label->Draw();
}


const char* get_tpc_label(const TString& name)
{
  if (name.Contains("_0_xsec")) return "East TPC";
  if (name.Contains("_1_xsec")) return "West TPC";

  if (name.EndsWith("_0")) return "East TPC";
  if (name.EndsWith("_1")) return "West TPC";

  return "";
}


bool keep_plot(const TString& name)
{
  static std::vector<TString> keep = {

    // Plane 0 / U
    "CzyHist_0_0",
    "CzyHist_0_1",

    "zy_0_0",
    "zy_0_1",

    "zynhits_0_0",
    "zynhits_0_1",

    "zydqdx_0_0",
    "zydqdx_0_1",

    "dqdxHist_0_0",
    "dqdxHist_0_1",

    // Plane 1 / V
    "CzyHist_1_0",
    "CzyHist_1_1",

    "zy_1_0",
    "zy_1_1",

    "zynhits_1_0",
    "zynhits_1_1",

    "zydqdx_1_0",
    "zydqdx_1_1",

    "dqdxHist_1_0",
    "dqdxHist_1_1",

    // Plane 2 / Collection
    "CzyHist_2_0",
    "CzyHist_2_1",

    "zy_2_0",
    "zy_2_1",

    "zynhits_2_0",
    "zynhits_2_1",

    "zydqdx_2_0",
    "zydqdx_2_1",

    "dqdxHist_2_0",
    "dqdxHist_2_1"
  };

  for (const auto& k : keep) {
    if (name == k) return true;
  }

  return false;
}


static TString sanitize(TString s)
{
  s.ReplaceAll("/", "__");
  s.ReplaceAll(" ", "_");
  return s;
}


void draw_dqdx_histogram(TH1* h1, const TString& out)
{
  if (!h1) return;

  TCanvas* c = new TCanvas("c_dqdx", "c_dqdx", 1100, 850);
  c->cd();

  // Use a clean ROOT-like style, not SBND style.
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1110);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatFont(42);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatColor(kWhite);

  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.08);

  h1->SetStats(1);

  // Force clean labels.
  h1->SetTitle("");
  h1->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");
  h1->GetYaxis()->SetTitle("number of hits");

  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetTitleSize(0.045);

  h1->GetXaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetLabelSize(0.04);

  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleOffset(1.45);

  h1->SetLineColor(kBlue);
  h1->SetLineWidth(1);

  h1->Draw("HIST");

  // // Calculate median from the histogram (CHANGE PENDING)
  // double q = 0.5;
  // double median = 0.0;
  // h1->GetQuantiles(1, &median, &q);

  // Manually draw a centred title.
  TLatex title;
  title.SetNDC();
  title.SetTextFont(42);
  title.SetTextSize(0.05);
  title.SetTextAlign(22);
  title.DrawLatex(0.50, 0.96, "dqdx");

  c->Modified();
  c->Update();

  TPaveStats* st = (TPaveStats*)h1->FindObject("stats");
  if (st) {
    st->SetX1NDC(0.78);
    st->SetX2NDC(0.98);
    st->SetY1NDC(0.78);
    st->SetY2NDC(0.92);
    st->SetBorderSize(1);
    st->SetFillColor(kWhite);
    st->SetFillStyle(1001);
  }

  // // Display median from the histogram below stats box (CHANGE PENDING)
  // TLatex med_label;
  // med_label.SetNDC();
  // med_label.SetTextFont(42);
  // med_label.SetTextSize(0.025);
  // med_label.SetTextAlign(12);
  // med_label.DrawLatex(0.79, 0.745, Form("Median      %.1f", median));



  c->Modified();
  c->Update();
  c->SaveAs(out);

  delete c;
}


void draw_generic_th1(TH1* h1, const TString& name, const TString& out)
{
  if (!h1) return;

  TCanvas* c = new TCanvas("c_th1", "c_th1", 1100, 850);
  c->cd();

  gStyle->SetOptStat(0);

  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.10);
  c->SetTopMargin(0.08);

  h1->SetStats(0);
  h1->SetTitle(name);

  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetTitleSize(0.045);

  h1->GetXaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetLabelSize(0.04);

  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleOffset(1.3);

  h1->Draw("HIST");

  c->Modified();
  c->Update();
  c->SaveAs(out);

  delete c;
}


void draw_th2_sbnd(TH2F* h2, const TString& name, const TString& out)
{
  if (!h2) return;

  TCanvas* c = new TCanvas("c_th2", "c_th2", 1100, 850);
  c->cd();

  sbndstyle::SetSBNDStyle();
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  sbndstyle::SeaPalette();

  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  c->SetRightMargin(0.19);
  c->SetTopMargin(0.08);

  if (name.BeginsWith("CzyHist")) {
    SetHistSBND(h2, "",
                "Z Coordinate [cm]",
                "Y Coordinate [cm]",
                "YZ correction factor");

    h2->GetZaxis()->SetRangeUser(0.8, 1.12);
  }
  else if (name.BeginsWith("zy_")) {
    SetHistSBND(h2, "",
                "Z Coordinate [cm]",
                "Y Coordinate [cm]",
                "Median dQ/dx [ADC/cm]");

    h2->GetZaxis()->SetRangeUser(500, 1300);
  }
  else if (name.BeginsWith("zynhits")) {
    SetHistSBND(h2, "",
                "Z Coordinate [cm]",
                "Y Coordinate [cm]",
                "number of hits");
  }
  else if (name.BeginsWith("zydqdx")) {
    SetHistSBND(h2, "",
                "Z Coordinate [cm]",
                "Y Coordinate [cm]",
                "dQ/dx [ADC/cm]");
  }
  else {
    SetHistSBND(h2, "",
                "Z Coordinate [cm]",
                "Y Coordinate [cm]",
                h2->GetZaxis()->GetTitle());
  }

  h2->SetStats(0);
  h2->Draw("COLZ");

  DrawLabel(get_tpc_label(name), 0.7, 0.15, 0.95, kBlack, 12);
  DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);

  c->Modified();
  c->Update();
  c->SaveAs(out);

  delete c;
}


void save_histogram(TObject* obj, const TString& outdir, int index)
{
  if (!obj) return;

  TString name = obj->GetName();
  TString out = TString::Format(
      "%s/%03d_%s.png",
      outdir.Data(),
      index,
      sanitize(name).Data()
  );

  if (obj->InheritsFrom(TH2::Class())) {
    TH2F* h2 = (TH2F*)obj;
    draw_th2_sbnd(h2, name, out);
    return;
  }

  if (obj->InheritsFrom(TH1::Class())) {
    TH1* h1 = (TH1*)obj;

    if (name.BeginsWith("dqdxHist")) {
      draw_dqdx_histogram(h1, out);
    }
    else {
      draw_generic_th1(h1, name, out);
    }

    return;
  }
}


void walk_directory(TDirectory* dir, const TString& outdir, int& index)
{
  if (!dir) return;

  TIter next(dir->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (!obj) continue;

    if (obj->InheritsFrom(TDirectory::Class())) {
      walk_directory((TDirectory*)obj, outdir, index);
      delete obj;
      continue;
    }

    if (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(TH2::Class())) {
      TString name = obj->GetName();

      if (!keep_plot(name)) {
        delete obj;
        continue;
      }

      save_histogram(obj, outdir, index);
      ++index;
    }

    delete obj;
  }
}


void export_hists_to_png(const char* rootfile,
                         const char* outdir = "png_out")
{
  gSystem->mkdir(outdir, true);

  TFile f(rootfile, "READ");
  if (f.IsZombie()) {
    std::cerr << "ERROR: cannot open " << rootfile << "\n";
    return;
  }

  int index = 0;
  walk_directory(&f, outdir, index);

  f.Close();

  std::cout << "Exported " << index << " histograms to " << outdir << "\n";
}