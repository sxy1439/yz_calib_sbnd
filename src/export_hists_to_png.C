
// Version 1 
// // export_hists_to_png.C
// #include <TFile.h>
// #include <TKey.h>
// #include <TClass.h>
// #include <TSystem.h>
// #include <TCanvas.h>
// #include <TString.h>
// #include <TH1.h>
// #include <TH2.h>
// #include <iostream>

// static TString sanitize(TString s) {
//   s.ReplaceAll("/", "__");
//   s.ReplaceAll(" ", "_");
//   return s;
// }

// void export_hists_to_png(const char* rootfile,
//                          const char* outdir = "png_out")
// {
//   gSystem->mkdir(outdir, true);

//   TFile f(rootfile, "READ");
//   if (f.IsZombie()) {
//     std::cerr << "ERROR: cannot open " << rootfile << "\n";
//     return;
//   }

//   TIter next(f.GetListOfKeys());
//   TKey* key;
//   int i = 0;

//   while ((key = (TKey*)next())) {
//     TObject* obj = key->ReadObj();
//     if (!obj) continue;

//     TString name = obj->GetName();
//     TString out  = TString::Format("%s/%03d_%s.png", outdir, i++, sanitize(name).Data());

//     TCanvas c("c", "c", 1100, 850);

//     if (obj->InheritsFrom(TH2::Class())) {
//       obj->Draw("COLZ");
//       c.Modified(); c.Update();
//       c.SaveAs(out);
//       std::cout << "Saved TH2: " << out << "\n";
//     } else if (obj->InheritsFrom(TH1::Class())) {
//       obj->Draw("HIST");
//       c.Modified(); c.Update();
//       c.SaveAs(out);
//       std::cout << "Saved TH1: " << out << "\n";
//     } else {
//       // ignore non-hist objects
//     }
//     delete obj;
//   }

//   std::cout << "Done. PNGs in: " << outdir << "\n";
// }













// void change_max_value(const char* rootfile_run1, const char* rootfile_run2, const char* outdir = "png_out")
// {
//   gSystem->mkdir(outdir, true);
// // This is for the input root file to get the histograms and change the max value of the histogram to 
//   // To do lst: rename all the bits to match the inputs

//   TFile f(rootfile, "READ");
//   if (f.IsZombie()) {
//     std::cerr << "ERROR: cannot open " << rootfile << "\n";
//     return;
//   }
//   TIter next(f.GetListOfKeys());
//   TKey* key;
//   int i = 0;

//   while ((key = (TKey*)next())) {
//     TObject* obj = key->ReadObj();
//     if (!obj) continue;

//     TString name = obj->GetName();
//     TString out  = TString::Format("%s/%03d_%s.png", outdir, i++, sanitize(name).Data());

// // Rename to match the input name also make key 1 and key 2 to istinguish betwwen run1 and run 2 keys
// // This is for the 2nd root file to compare with the 1st root file

//   TFile f(rootfile, "READ");
//   if (f.IsZombie()) {
//     std::cerr << "ERROR: cannot open " << rootfile << "\n";
//     return;
//   }

//   TIter next(f.GetListOfKeys());
//   TKey* key;
//   int i = 0;

//   while ((key = (TKey*)next())) {
//     TObject* obj = key->ReadObj();
//     if (!obj) continue;

//     TString name = obj->GetName();
//     TString out  = TString::Format("%s/%03d_%s.png", outdir, i++, sanitize(name).Data());



// // Make a vector here to store the colz values I get from the histograms and then use that vector to set the max value for the histograms in the 2nd root file
//     TCanvas c("c", "c", 1100, 850);
// // Make an if statement to check for the name of the historgram that it matches the plots I want
// // 



//     if (obj->InheritsFrom(TH2::Class())) {
//       obj->Draw("COLZ"); 
//       //This is where I also get the maxiumum value for colz
//       c.Modified(); c.Update();
//       c.SaveAs(out);
//       std::cout << "Saved TH2: " << out << "\n";
//     } else if (obj->InheritsFrom(TH1::Class())) {
//       obj->Draw("HIST");
//       c.Modified(); c.Update();
//       c.SaveAs(out);
//       std::cout << "Saved TH1: " << out << "\n";
//     } else {
//       // ignore non-hist objects
//     }
// // Open the 2nd historgram and plot it by setting the value of the ith component of the vector from the 1st histogrm using an  loop and export to png as before

//     delete obj;
//   }

//   std::cout << "Done. PNGs in: " << outdir << "\n";
// }



// Version 2





// #include <TFile.h>
// #include <TKey.h>
// #include <TClass.h>
// #include <TSystem.h>
// #include <TCanvas.h>
// #include <TString.h>
// #include <TH1.h>
// #include <TH2.h>
// #include <TStyle.h>
// #include <TPaveStats.h>
// #include <iostream>
// #include <TPaletteAxis.h>
// #include <TPaveStats.h>

// static TString sanitize(TString s) {
//   s.ReplaceAll("/", "__");
//   s.ReplaceAll(" ", "_");
//   return s;
// }

// void export_hists_to_png(const char* rootfile,
//                          const char* outdir = "png_out")
// {
//   gSystem->mkdir(outdir, true);

//   TFile f(rootfile, "READ");
//   if (f.IsZombie()) {
//     std::cerr << "ERROR: cannot open " << rootfile << "\n";
//     return;
//   }

//   gStyle->SetOptStat(1111);

//   TIter next(f.GetListOfKeys());
//   TKey* key;
//   int i = 0;

//   while ((key = (TKey*)next())) {
//     TObject* obj = key->ReadObj();
//     if (!obj) continue;

//     TString name = obj->GetName();
//     TString out  = TString::Format("%s/%03d_%s.png", outdir, i, sanitize(name).Data());

//     TCanvas* c = new TCanvas("c", "c", 1100, 850);
//     c->cd();

//     c->SetLeftMargin(0.12);
//     c->SetBottomMargin(0.12);
//     c->SetRightMargin(0.22);
//     c->SetTopMargin(0.10);

//     if (obj->InheritsFrom(TH2::Class())) {
//       TH2* h2 = (TH2*)obj;

//       h2->SetTitle(name);
//       h2->GetXaxis()->SetTitleSize(0.045);
//       h2->GetYaxis()->SetTitleSize(0.045);
//       h2->GetZaxis()->SetTitleSize(0.045);

//       h2->GetXaxis()->SetLabelSize(0.04);
//       h2->GetYaxis()->SetLabelSize(0.04);
//       h2->GetZaxis()->SetLabelSize(0.04);

//       h2->GetXaxis()->SetTitleOffset(1.1);
//       h2->GetYaxis()->SetTitleOffset(1.3);
//       h2->GetZaxis()->SetTitleOffset(1.3);

//       h2->Draw("COLZ");
//       c->Update();



//       // --- Move palette (color bar) ---
//       TPaletteAxis* palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
//       if (palette) {
//         palette->SetX1NDC(0.82);
//         palette->SetX2NDC(0.86);
//         palette->SetY1NDC(0.15);
//         palette->SetY2NDC(0.90);
//       }


//       TPaveStats* st = (TPaveStats*)h2->FindObject("stats");
//       if (st) {
//         st->SetX1NDC(0.82);
//         st->SetX2NDC(0.98);
//         st->SetY1NDC(0.90);
//         st->SetY2NDC(0.98);
//         st->SetBorderSize(1);
//         st->SetFillStyle(1001);
//       }

//       c->Modified();
//       c->Update();
//       c->SaveAs(out);
//     }
//     else if (obj->InheritsFrom(TH1::Class())) {
//       TH1* h1 = (TH1*)obj;

//       h1->SetTitle(name);
//       h1->GetXaxis()->SetTitleSize(0.045);
//       h1->GetYaxis()->SetTitleSize(0.045);

//       h1->GetXaxis()->SetLabelSize(0.04);
//       h1->GetYaxis()->SetLabelSize(0.04);

//       h1->GetXaxis()->SetTitleOffset(1.1);
//       h1->GetYaxis()->SetTitleOffset(1.3);

//       h1->Draw();
//       c->Update();

//       TPaveStats* st = (TPaveStats*)h1->FindObject("stats");
//       if (st) {
//         st->SetX1NDC(0.78);
//         st->SetX2NDC(0.98);
//         st->SetY1NDC(0.75);
//         st->SetY2NDC(0.92);
//         st->SetBorderSize(1);
//         st->SetFillStyle(1001);
//       }

//       c->Modified();
//       c->Update();
//       c->SaveAs(out);
//     }

//     delete c;
//     delete obj;
//     ++i;
//   }

//   f.Close();
// }





// Version 3 

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
#include <iostream>

static TString sanitize(TString s) {
  s.ReplaceAll("/", "__");
  s.ReplaceAll(" ", "_");
  return s;
}

void save_histogram(TObject* obj, const TString& outdir, int index) {
  if (!obj) return;

  TString name = obj->GetName();
  TString out = TString::Format("%s/%03d_%s.png", outdir.Data(), index, sanitize(name).Data());

  TCanvas* c = new TCanvas("c", "c", 1100, 850);
  c->cd();

  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.22);
  c->SetTopMargin(0.08);

  if (obj->InheritsFrom(TH2::Class())) {
    TH2* h2 = (TH2*)obj;

    h2->SetTitle(name);
    h2->GetXaxis()->SetTitleSize(0.045);
    h2->GetYaxis()->SetTitleSize(0.045);
    h2->GetZaxis()->SetTitleSize(0.045);

    h2->GetXaxis()->SetLabelSize(0.04);
    h2->GetYaxis()->SetLabelSize(0.04);
    h2->GetZaxis()->SetLabelSize(0.04);

    h2->GetXaxis()->SetTitleOffset(1.1);
    h2->GetYaxis()->SetTitleOffset(1.3);
    h2->GetZaxis()->SetTitleOffset(1.35);

    h2->Draw("COLZ");
    c->Update();

    TPaletteAxis* palette =
      (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
    if (palette) {
      palette->SetX1NDC(0.80);
      palette->SetX2NDC(0.84);
      palette->SetY1NDC(0.14);
      palette->SetY2NDC(0.88);
    }

    TPaveStats* st = (TPaveStats*)h2->FindObject("stats");
    if (st) {
      st->SetX1NDC(0.80);
      st->SetX2NDC(0.98);
      st->SetY1NDC(0.89);
      st->SetY2NDC(0.98);
      st->SetBorderSize(1);
      st->SetFillStyle(1001);
    }

    c->Modified();
    c->Update();
    c->SaveAs(out);
  }
  else if (obj->InheritsFrom(TH1::Class())) {
    TH1* h1 = (TH1*)obj;

    h1->SetTitle(name);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetTitleSize(0.045);

    h1->GetXaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetLabelSize(0.04);

    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetYaxis()->SetTitleOffset(1.3);

    h1->Draw();
    c->Update();

    TPaveStats* st = (TPaveStats*)h1->FindObject("stats");
    if (st) {
      st->SetX1NDC(0.78);
      st->SetX2NDC(0.98);
      st->SetY1NDC(0.78);
      st->SetY2NDC(0.92);
      st->SetBorderSize(1);
      st->SetFillStyle(1001);
    }

    c->Modified();
    c->Update();
    c->SaveAs(out);
  }

  delete c;
}

void walk_directory(TDirectory* dir, const TString& outdir, int& index) {
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

  gStyle->SetOptStat(1111);

  int index = 0;
  walk_directory(&f, outdir, index);

  f.Close();

  std::cout << "Exported " << index << " histograms to " << outdir << "\n";
}