// // VERSION 2
// #include <TPaveText.h>
// #include <TText.h>
// #include <TFile.h>
// #include <TKey.h>
// #include <TDirectory.h>
// #include <TH1.h>
// #include <TH2.h>
// #include <TCanvas.h>
// #include <TString.h>
// #include <TSystem.h>
// #include <TPaveStats.h>
// #include <TPaletteAxis.h>
// #include <TStyle.h>
// #include <TLegend.h>
// #include <TROOT.h>
// #include <iostream>
// #include <algorithm>

// static TString sanitize(TString s) {
//   s.ReplaceAll("/", "__");
//   s.ReplaceAll(" ", "_");
//   return s;
// }

// // Use histogram's natural bin-content range, not display state
// double natural_min(TH2* h) {
//   if (!h) return 0.0;
//   return h->GetBinContent(h->GetMinimumBin());
// }

// double natural_max(TH2* h) {
//   if (!h) return 0.0;
//   return h->GetBinContent(h->GetMaximumBin());
// }

// void format_pad_2d(TH2* h) {
//   gPad->SetLeftMargin(0.12);
//   gPad->SetBottomMargin(0.12);
//   gPad->SetRightMargin(0.22);
//   gPad->SetTopMargin(0.08);

//   h->GetXaxis()->SetTitleSize(0.045);
//   h->GetYaxis()->SetTitleSize(0.045);
//   h->GetZaxis()->SetTitleSize(0.045);

//   h->GetXaxis()->SetLabelSize(0.04);
//   h->GetYaxis()->SetLabelSize(0.04);
//   h->GetZaxis()->SetLabelSize(0.04);

//   h->GetXaxis()->SetTitleOffset(1.10);
//   h->GetYaxis()->SetTitleOffset(1.30);
//   h->GetZaxis()->SetTitleOffset(1.35);
// }

// void adjust_palette_and_stats(TH2* h) {
//   gPad->Update();

//   TPaletteAxis* palette =
//     (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
//   if (palette) {
//     palette->SetX1NDC(0.80);
//     palette->SetX2NDC(0.84);
//     palette->SetY1NDC(0.14);
//     palette->SetY2NDC(0.88);
//   }

//   TPaveStats* st = (TPaveStats*)h->FindObject("stats");
//   if (st) {
//     st->SetX1NDC(0.80);
//     st->SetX2NDC(0.98);
//     st->SetY1NDC(0.89);
//     st->SetY2NDC(0.98);
//     st->SetBorderSize(1);
//     st->SetFillStyle(1001);
//   }

//   gPad->Modified();
//   gPad->Update();
// }

// void format_hist_1d(TH1* h) {
//   if (!h) return;

//   h->GetXaxis()->SetTitleSize(0.045);
//   h->GetYaxis()->SetTitleSize(0.045);

//   h->GetXaxis()->SetLabelSize(0.04);
//   h->GetYaxis()->SetLabelSize(0.04);

//   h->GetXaxis()->SetTitleOffset(1.10);
//   h->GetYaxis()->SetTitleOffset(1.30);

//   h->SetStats(0);
//   h->SetLineWidth(2);
// }

// void compare_histograms_2d(TH2* h_before, TH2* h_after, const TString& outname) {
//   if (!h_before || !h_after) return;

//   TH2* hb = (TH2*)h_before->Clone(TString(h_before->GetName()) + "_before_clone");
//   TH2* ha = (TH2*)h_after->Clone(TString(h_after->GetName()) + "_after_clone");


//   // Lines to overwrite plot titles with the histogram file names 
//   hb->SetTitle(h_before->GetName());
//   ha->SetTitle(h_after->GetName());


//   // Use ONLY the natural z-range of the BEFORE histogram
//   double zmin = natural_min(hb);
//   double zmax = natural_max(hb);

//   if (zmin == zmax) {
//     zmin -= 1.0;
//     zmax += 1.0;
//   }

//   hb->SetMinimum(zmin);
//   hb->SetMaximum(zmax);
//   ha->SetMinimum(zmin);
//   ha->SetMaximum(zmax);

//   TCanvas* c = new TCanvas("c_compare_2d", "comparison_2d", 1800, 800);
//   c->Divide(2, 1);

//   // LEFT: before
//   c->cd(1);
//   format_pad_2d(hb);
//   hb->Draw("COLZ");
//   adjust_palette_and_stats(hb);

//   // RIGHT: after
//   c->cd(2);
//   format_pad_2d(ha);
//   ha->Draw("COLZ");
//   adjust_palette_and_stats(ha);

//   c->SaveAs(outname);

//   delete hb;
//   delete ha;
//   delete c;
// }


// void compare_histograms_1d(TH1* h_before, TH1* h_after, const TString& outname) {
//   if (!h_before || !h_after) return;

//   TH1* hb = (TH1*)h_before->Clone(TString(h_before->GetName()) + "_before_clone");
//   TH1* ha = (TH1*)h_after->Clone(TString(h_after->GetName()) + "_after_clone");

//   TCanvas* c = new TCanvas("c_compare_1d", "comparison_1d", 1100, 850);
//   c->cd();

//   gPad->SetLeftMargin(0.12);
//   gPad->SetBottomMargin(0.12);
//   gPad->SetRightMargin(0.06);
//   gPad->SetTopMargin(0.08);

//   format_hist_1d(hb);
//   format_hist_1d(ha);

//   hb->SetLineColor(kBlack);
//   ha->SetLineColor(kRed);

//   hb->SetMarkerColor(kBlack);
//   ha->SetMarkerColor(kRed);

//   double ymax = std::max(hb->GetMaximum(), ha->GetMaximum());
//   hb->SetMaximum(1.15 * ymax);

//   hb->Draw("HIST");
//   ha->Draw("HIST SAME");

//   TLegend* leg = new TLegend(0.68, 0.78, 0.90, 0.90);
//   leg->AddEntry(hb, "Run 1", "l");
//   leg->AddEntry(ha, "Run 2", "l");
//   leg->SetBorderSize(1);
//   leg->SetFillStyle(1001);
//   leg->Draw();



//   TPaveText* stats = new TPaveText(0.68, 0.58, 0.90, 0.76, "NDC");
//   stats->SetBorderSize(1);
//   stats->SetFillStyle(1001);
//   stats->SetTextAlign(12);

//   stats->AddText(Form("Run 1:"));
//   stats->AddText(Form("  Entries = %.0f", hb->GetEntries()));
//   stats->AddText(Form("  Mean = %.2f", hb->GetMean()));
//   stats->AddText(Form("  Std Dev = %.2f", hb->GetStdDev()));

//   stats->AddText("");

//   stats->AddText(Form("Run 2:"));
//   stats->AddText(Form("  Entries = %.0f", ha->GetEntries()));
//   stats->AddText(Form("  Mean = %.2f", ha->GetMean()));
//   stats->AddText(Form("  Std Dev = %.2f", ha->GetStdDev()));

//   stats->Draw();

//   c->Modified();
//   c->Update();
//   c->SaveAs(outname);

//   delete stats;
//   delete leg;
//   delete hb;
//   delete ha;
//   delete c;
// }

// void walk_and_compare(TDirectory* d_before,
//                       TDirectory* d_after,
//                       const TString& outdir,
//                       int& index)
// {
//   if (!d_before || !d_after) return;

//   TIter next(d_before->GetListOfKeys());
//   TKey* key;

//   while ((key = (TKey*)next())) {
//     TObject* obj_before = key->ReadObj();
//     if (!obj_before) continue;

//     TString name = obj_before->GetName();
//     TObject* obj_after = d_after->Get(name);

//     if (!obj_after) {
//       delete obj_before;
//       continue;
//     }

//     if (obj_before->InheritsFrom(TDirectory::Class())) {
//       if (obj_after->InheritsFrom(TDirectory::Class())) {
//         walk_and_compare((TDirectory*)obj_before,
//                          (TDirectory*)obj_after,
//                          outdir,
//                          index);
//       }
//       delete obj_before;
//       continue;
//     }

//     TString out = TString::Format("%s/%03d_%s.png",
//                                   outdir.Data(),
//                                   index,
//                                   sanitize(name).Data());

//     if (obj_before->InheritsFrom(TH2::Class()) &&
//         obj_after->InheritsFrom(TH2::Class())) {
//       compare_histograms_2d((TH2*)obj_before, (TH2*)obj_after, out);
//       ++index;
//     }
//     else if (obj_before->InheritsFrom(TH1::Class()) &&
//              obj_after->InheritsFrom(TH1::Class()) &&
//              !obj_before->InheritsFrom(TH2::Class()) &&
//              !obj_after->InheritsFrom(TH2::Class())) {
//       compare_histograms_1d((TH1*)obj_before, (TH1*)obj_after, out);
//       ++index;
//     }

//     delete obj_before;
//   }
// }

// void compare_root_files(const char* before_file,
//                         const char* after_file,
//                         const char* outdir = "comparison_pngs")
// {
//   gSystem->mkdir(outdir, true);
//   gStyle->SetOptStat(1111);

//   TFile f_before(before_file, "READ");
//   TFile f_after(after_file, "READ");

//   if (f_before.IsZombie() || f_after.IsZombie()) {
//     std::cerr << "ERROR: could not open one or both ROOT files.\n";
//     return;
//   }

//   int index = 0;
//   walk_and_compare(&f_before, &f_after, outdir, index);

//   std::cout << "Saved " << index << " comparison plots to " << outdir << "\n";

//   f_before.Close();
//   f_after.Close();
// }























// VERSION 3 Where I implemented plot renormalisation and more detailed stats in the 1D comparison. 
// I also changed the stats box position and size to be more consistent across different histograms, 
// and added a border to it for better visibility. I also added a message at the end to report how many comparisons were saved.

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

void format_pad_2d(TH2* h) {
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.22);
  gPad->SetTopMargin(0.08);

  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetZaxis()->SetTitleSize(0.045);

  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetZaxis()->SetLabelSize(0.04);

  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.30);
  h->GetZaxis()->SetTitleOffset(1.35);
}

void adjust_palette_and_stats(TH2* h) {
  gPad->Update();

  TPaletteAxis* palette =
    (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
  if (palette) {
    palette->SetX1NDC(0.80);
    palette->SetX2NDC(0.84);
    palette->SetY1NDC(0.14);
    palette->SetY2NDC(0.88);
  }

  TPaveStats* st = (TPaveStats*)h->FindObject("stats");
  if (st) {
    st->SetX1NDC(0.80);
    st->SetX2NDC(0.98);
    st->SetY1NDC(0.89);
    st->SetY2NDC(0.98);
    st->SetBorderSize(1);
    st->SetFillStyle(1001);
  }

  gPad->Modified();
  gPad->Update();
}

void format_hist_1d(TH1* h) {
  if (!h) return;

  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);

  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);

  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.30);

  h->SetStats(0);
  h->SetLineWidth(2);
}

void compare_histograms_2d(TH2* h_before, TH2* h_after, const TString& outname) {
  if (!h_before || !h_after) return;

  TH2* hb = (TH2*)h_before->Clone(TString(h_before->GetName()) + "_Old");
  TH2* ha = (TH2*)h_after->Clone(TString(h_after->GetName()) + "_New");

  // Keep original titles from the ROOT file
  // Force both plots to use Run 1 natural z range
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

  TCanvas* c = new TCanvas("c_compare_2d", "comparison_2d", 1800, 800);
  c->Divide(2, 1);

  c->cd(1);
  format_pad_2d(hb);
  hb->Draw("COLZ");
  adjust_palette_and_stats(hb);

  c->cd(2);
  format_pad_2d(ha);
  ha->Draw("COLZ");
  adjust_palette_and_stats(ha);

  c->SaveAs(outname);

  delete hb;
  delete ha;
  delete c;
}

void draw_histograms_1d_counts(TH1* h_before, TH1* h_after, const TString& outname) {
  if (!h_before || !h_after) return;

  TH1* hb = (TH1*)h_before->Clone(TString(h_before->GetName()) + "_old_counts");
  TH1* ha = (TH1*)h_after->Clone(TString(h_after->GetName()) + "_new_counts");

  TCanvas* c = new TCanvas("c_compare_1d_counts", "comparison_1d_counts", 1100, 850);
  c->cd();

  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.08);

  format_hist_1d(hb);
  format_hist_1d(ha);

  hb->SetLineColor(kBlack);
  ha->SetLineColor(kRed);

  hb->SetMarkerColor(kBlack);
  ha->SetMarkerColor(kRed);

  hb->GetYaxis()->SetTitle("Number of hits");

  double ymax = std::max(hb->GetMaximum(), ha->GetMaximum());
  hb->SetMaximum(1.15 * ymax);

  hb->Draw("HIST");
  ha->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.68, 0.78, 0.90, 0.90);
  leg->AddEntry(hb, "Old", "l");
  leg->AddEntry(ha, "New", "l");
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->Draw();

  TPaveText* stats = new TPaveText(0.16, 0.74, 0.46, 0.84, "NDC");
  stats->SetBorderSize(1);
  stats->SetFillStyle(1001);
  stats->SetTextAlign(12);
  stats->SetTextSize(0.022);

  stats->AddText(Form("Old: N=%.3g, #mu=%.1f, #sigma=%.1f",
                      hb->GetEntries(), hb->GetMean(), hb->GetStdDev()));
  stats->AddText(Form("New: N=%.3g, #mu=%.1f, #sigma=%.1f",
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

  TH1* hb = (TH1*)h_before->Clone(TString(h_before->GetName()) + "_old_shape");
  TH1* ha = (TH1*)h_after->Clone(TString(h_after->GetName()) + "_new_shape");

  TCanvas* c = new TCanvas("c_compare_1d_shape", "comparison_1d_shape", 1100, 850);
  c->cd();

  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.08);

  format_hist_1d(hb);
  format_hist_1d(ha);

  hb->SetLineColor(kBlack);
  ha->SetLineColor(kRed);

  hb->SetMarkerColor(kBlack);
  ha->SetMarkerColor(kRed);

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
  leg->AddEntry(hb, "Old", "l");
  leg->AddEntry(ha, "New", "l");
  leg->SetBorderSize(1);
  leg->SetFillStyle(1001);
  leg->Draw();

  TPaveText* stats = new TPaveText(0.16, 0.74, 0.46, 0.84, "NDC");
  stats->SetBorderSize(1);
  stats->SetFillStyle(1001);
  stats->SetTextAlign(12);
  stats->SetTextSize(0.022);

  stats->AddText(Form("Old: #mu=%.1f, #sigma=%.1f",
                      hb->GetMean(), hb->GetStdDev()));
  stats->AddText(Form("New: #mu=%.1f, #sigma=%.1f",
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
                      int& index)
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
                         index);
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

      compare_histograms_2d((TH2*)obj_before, (TH2*)obj_after, out);
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
  gStyle->SetOptStat(1111);

  TFile f_before(before_file, "READ");
  TFile f_after(after_file, "READ");

  if (f_before.IsZombie() || f_after.IsZombie()) {
    std::cerr << "ERROR: could not open one or both ROOT files.\n";
    return;
  }

  int index = 0;
  walk_and_compare(&f_before, &f_after, outdir, index);

  std::cout << "Saved " << index << " comparison groups to " << outdir << "\n";

  f_before.Close();
  f_after.Close();
}