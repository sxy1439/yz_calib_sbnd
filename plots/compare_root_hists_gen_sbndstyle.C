#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TPaveStats.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TLatex.h"

#define SBNDSTYLE_ENABLE_AUTOMATICALLY 0
#include "SBNDStyle.h"
#include <TLegend.h>
#include <TPaveText.h>
#include <iostream>
#include <algorithm>

static TString sanitize(TString s) {
  s.ReplaceAll("/", "__");
  s.ReplaceAll(" ", "_");
  return s;
}

// Use histogram's natural bin-content range, not display state
double natural_min(TH2* h) {
  if (!h) return 0.0;
  return h->GetBinContent(h->GetMinimumBin());
}

double natural_max(TH2* h) {
  if (!h) return 0.0;
  return h->GetBinContent(h->GetMaximumBin());
}

void apply_sbnd_plot_style() {
  sbndstyle::SetSBNDStyle();
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  sbndstyle::SeaPalette();
}

void draw_compare_label(const char* text="text", float tscale=0.8,
                        double xloc=0.9, double yloc=0.95,
                        int tcolor=kBlack, int align=12) {
  TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  label->SetTextColor(tcolor);
  label->SetTextFont(42);
  label->SetNDC();
  label->SetTextSize(2.0/30.0);
  label->SetTextAlign(align);
  label->Draw();
}

void format_pad_2d(TH2* h) {
  if (!h) return;

  // Match the pad geometry used in plot_yz_planes.C
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.19);
  gPad->SetTopMargin(0.08);

  h->SetStats(0);
  h->SetTitle("");

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();

  h->GetXaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetZaxis()->SetTitleSize(0.055);

  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetLabelSize(0.05);

  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.10);
  h->GetZaxis()->SetTitleOffset(1.17);
}


void apply_2d_axis_labels(TH2* h) {
  if (!h) return;

  // Match plot_yz_planes.C: no ROOT histogram title above the pad.
  h->SetTitle("");
  h->GetXaxis()->SetTitle("Z Coordinate [cm]");
  h->GetYaxis()->SetTitle("Y Coordinate [cm]");

  TString name = h->GetName();
  if (name.Contains("CzyHist")) {
    h->GetZaxis()->SetTitle("YZ correction factor");
  }
  else if (name.Contains("zynhits")) {
    h->GetZaxis()->SetTitle("Number of hits");
  }
  else if (name.Contains("zy") || name.Contains("dqdx")) {
    h->GetZaxis()->SetTitle("Median dQ/dx [ADC/cm]");
  }
}

void draw_sample_label(const char* filename) {
  TString fname = filename;
  if (fname.Contains("data") || fname.Contains("Data")) {
    draw_compare_label("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
  }
  else {
    draw_compare_label("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
  }
}

void adjust_palette(TH2* h) {
  if (!h) return;
  gPad->Update();

  TPaletteAxis* palette =
    (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
  if (palette) {
    // Slightly inside the right margin, matching the cleaner plot_yz_planes look.
    palette->SetX1NDC(0.83);
    palette->SetX2NDC(0.87);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.92);
  }

  gPad->Modified();
  gPad->Update();
}

void format_hist_1d(TH1* h) {
  if (!h) return;

  h->SetStats(0);
  h->SetTitle("");

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

  h->GetXaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleSize(0.055);

  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);

  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.10);

  h->SetLineWidth(2);
}

void compare_histograms_2d(TH2* h_before, TH2* h_after, const TString& outname, const char* gen1_file, const char* gen2_file) {
  if (!h_before || !h_after) return;

  TH2* hb = (TH2*)h_before->Clone(TString(h_before->GetName()) + "_gen1");
  TH2* ha = (TH2*)h_after->Clone(TString(h_after->GetName()) + "_gen2");

  apply_2d_axis_labels(hb);
  apply_2d_axis_labels(ha);

  // Force both plots to use Gen 1 natural z range
  double zmin = natural_min(hb);
  double zmax = natural_max(hb);

  if (zmin == zmax) {
    zmin -= 1.0;
    zmax += 1.0;
  }

  hb->SetMinimum(zmin);
  hb->SetMaximum(zmax);
  ha->SetMinimum(zmin);
  ha->SetMaximum(zmax);

  apply_sbnd_plot_style();

  TCanvas* c = new TCanvas("c_compare_2d", "comparison_2d", 1800, 800);
  c->Divide(2, 1);

  c->cd(1);
  format_pad_2d(hb);
  hb->Draw("COLZ");
  draw_compare_label("Gen 1", 0.7, 0.15, 0.95, kBlack, 12);
  draw_sample_label(gen1_file);
  adjust_palette(hb);

  c->cd(2);
  format_pad_2d(ha);
  ha->Draw("COLZ");
  draw_compare_label("Gen 2", 0.7, 0.15, 0.95, kBlack, 12);
  draw_sample_label(gen2_file);
  adjust_palette(ha);

  c->SaveAs(outname);

  delete hb;
  delete ha;
  delete c;
}

void draw_histograms_1d_counts(TH1* h_before, TH1* h_after, const TString& outname) {
  if (!h_before || !h_after) return;

  TH1* hb = (TH1*)h_before->Clone(TString(h_before->GetName()) + "_gen1_counts");
  TH1* ha = (TH1*)h_after->Clone(TString(h_after->GetName()) + "_gen2_counts");

  TCanvas* c = new TCanvas("c_compare_1d_counts", "comparison_1d_counts", 1100, 850);
  c->cd();

  apply_sbnd_plot_style();

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.08);

  format_hist_1d(hb);
  format_hist_1d(ha);

  hb->SetLineColor(kBlack);
  ha->SetLineColor(sbndstyle::colors::kOkabeItoVermilion);

  hb->SetMarkerColor(kBlack);
  ha->SetMarkerColor(sbndstyle::colors::kOkabeItoVermilion);

  hb->GetYaxis()->SetTitle("Number of hits");

  double ymax = std::max(hb->GetMaximum(), ha->GetMaximum());
  hb->SetMaximum(1.15 * ymax);

  hb->Draw("HIST");
  ha->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.68, 0.78, 0.90, 0.90);
  leg->AddEntry(hb, "Gen 1", "l");
  leg->AddEntry(ha, "Gen 2", "l");
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->Draw();

  TPaveText* stats = new TPaveText(0.16, 0.74, 0.46, 0.84, "NDC");
  stats->SetBorderSize(1);
  stats->SetFillStyle(1001);
  stats->SetTextAlign(12);
  stats->SetTextSize(0.022);

  stats->AddText(Form("Gen 1: N=%.3g, #mu=%.1f, #sigma=%.1f",
                      hb->GetEntries(), hb->GetMean(), hb->GetStdDev()));
  stats->AddText(Form("Gen 2: N=%.3g, #mu=%.1f, #sigma=%.1f",
                      ha->GetEntries(), ha->GetMean(), ha->GetStdDev()));
  stats->Draw();

  c->Modified();
  c->Update();
  c->SaveAs(outname);

  delete stats;
  delete leg;
  delete hb;
  delete ha;
  delete c;
}

void draw_histograms_1d_shape(TH1* h_before, TH1* h_after, const TString& outname) {
  if (!h_before || !h_after) return;

  TH1* hb = (TH1*)h_before->Clone(TString(h_before->GetName()) + "_gen1_shape");
  TH1* ha = (TH1*)h_after->Clone(TString(h_after->GetName()) + "_gen2_shape");

  TCanvas* c = new TCanvas("c_compare_1d_shape", "comparison_1d_shape", 1100, 850);
  c->cd();

  apply_sbnd_plot_style();

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.08);

  format_hist_1d(hb);
  format_hist_1d(ha);

  hb->SetLineColor(kBlack);
  ha->SetLineColor(sbndstyle::colors::kOkabeItoVermilion);

  hb->SetMarkerColor(kBlack);
  ha->SetMarkerColor(sbndstyle::colors::kOkabeItoVermilion);

  double int_before = hb->Integral();
  double int_after  = ha->Integral();

  if (int_before > 0) hb->Scale(1.0 / int_before);
  if (int_after  > 0) ha->Scale(1.0 / int_after);

  hb->GetYaxis()->SetTitle("Fraction of hits");

  double ymax = std::max(hb->GetMaximum(), ha->GetMaximum());
  hb->SetMaximum(1.15 * ymax);

  hb->Draw("HIST");
  ha->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.68, 0.78, 0.90, 0.90);
  leg->AddEntry(hb, "Gen 1", "l");
  leg->AddEntry(ha, "Gen 2", "l");
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->Draw();

  TPaveText* stats = new TPaveText(0.16, 0.74, 0.46, 0.84, "NDC");
  stats->SetBorderSize(1);
  stats->SetFillStyle(1001);
  stats->SetTextAlign(12);
  stats->SetTextSize(0.022);

  stats->AddText(Form("Gen 1: #mu=%.1f, #sigma=%.1f",
                      hb->GetMean(), hb->GetStdDev()));
  stats->AddText(Form("Gen 2: #mu=%.1f, #sigma=%.1f",
                      ha->GetMean(), ha->GetStdDev()));
  stats->Draw();

  c->Modified();
  c->Update();
  c->SaveAs(outname);

  delete stats;
  delete leg;
  delete hb;
  delete ha;
  delete c;
}









void walk_and_compare(TDirectory* d_before,
                      TDirectory* d_after,
                      const TString& outdir,
                      int& index,
                      const char* gen1_file,
                      const char* gen2_file)
{
  if (!d_before || !d_after) return;

  TIter next(d_before->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)next())) {
    TObject* obj_before = key->ReadObj();
    if (!obj_before) continue;

    TString name = obj_before->GetName();
    TObject* obj_after = d_after->Get(name);

    if (!obj_after) {
      delete obj_before;
      continue;
    }

    if (obj_before->InheritsFrom(TDirectory::Class())) {
      if (obj_after->InheritsFrom(TDirectory::Class())) {
        walk_and_compare((TDirectory*)obj_before,
                         (TDirectory*)obj_after,
                         outdir,
                         index,
                         gen1_file,
                         gen2_file);
      }
      delete obj_before;
      continue;
    }

    if (obj_before->InheritsFrom(TH2::Class()) &&
        obj_after->InheritsFrom(TH2::Class())) {

      TString out = TString::Format("%s/%03d_%s.png",
                                    outdir.Data(),
                                    index,
                                    sanitize(name).Data());

      compare_histograms_2d((TH2*)obj_before, (TH2*)obj_after, out, gen1_file, gen2_file);
      ++index;
    }
    else if (obj_before->InheritsFrom(TH1::Class()) &&
             obj_after->InheritsFrom(TH1::Class()) &&
             !obj_before->InheritsFrom(TH2::Class()) &&
             !obj_after->InheritsFrom(TH2::Class())) {

      TString out_counts = TString::Format("%s/%03d_%s_counts.png",
                                           outdir.Data(),
                                           index,
                                           sanitize(name).Data());
      ++index;

      TString out_shape = TString::Format("%s/%03d_%s_shape.png",
                                          outdir.Data(),
                                          index,
                                          sanitize(name).Data());
      ++index;

      draw_histograms_1d_counts((TH1*)obj_before, (TH1*)obj_after, out_counts);
      draw_histograms_1d_shape((TH1*)obj_before, (TH1*)obj_after, out_shape);
    }

    delete obj_before;
  }
}


void compare_root_files(const char* before_file,
                        const char* after_file,
                        const char* outdir = "comparison_pngs")
{
  gSystem->mkdir(outdir, true);
  apply_sbnd_plot_style();

  TFile f_before(before_file, "READ");
  TFile f_after(after_file, "READ");

  if (f_before.IsZombie() || f_after.IsZombie()) {
    std::cerr << "ERROR: could not open one or both ROOT files.\n";
    return;
  }

  int index = 0;
  walk_and_compare(&f_before, &f_after, outdir, index, before_file, after_file);

  std::cout << "Saved " << index << " comparison groups to " << outdir << "\n";

  f_before.Close();
  f_after.Close();
}
