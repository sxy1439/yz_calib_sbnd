#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TLatex.h"

void SetHist(TH1F *h, const char* title, const char* xlabel, const char* ylabel){
  h->SetTitle(title);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.1);   
}

void SetAxis(TH1F *h, double xmin, double xmax, double ymin, double ymax){
  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetRangeUser(ymin, ymax);
}

void DrawSBNDLabel(const char* text="text", float tscale=0.8, double xloc=0.9, double yloc=0.95){
  //TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{SBND Preliminary}", tscale));
  TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  prelim->SetTextColor(kGray+2);
  prelim->SetTextFont(132);
  prelim->SetNDC(); 
  prelim->SetTextSize(2.0/30.0); 
  prelim->SetTextAlign(32); 
  prelim->Draw();
}

void DrawLabel(const char* text="your text", float tscale=0.8, double xloc=0.9, double yloc=0.95, int tcolor=kGray+1, int align=12){
  TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  label->SetTextColor(tcolor);
  label->SetTextFont(132); // new times roman
  label->SetNDC(); 
  label->SetTextSize(2.0/30.0); 
  label->SetTextAlign(align); 
  label->Draw();
}

void SetLineStyle(TH1F *h, int color){
  h->SetLineColor(color);
  h->SetLineWidth(3);
}



void plot_dqdx(const char* instring){

  TFile *file_nosce = new TFile("/exp/sbnd/app/users/yadav/Calibration/YZ_X_Calib/Median/split_macros/output_files/dQdx_data_dev_nosce.root", "READ");
  //TFile *file_sce = new TFile("/exp/sbnd/app/users/yadav/Calibration/YZ_X_Calib/Median/split_macros/output_files/dQdx_data_dev_sce.root", "READ");
  TFile *file_sce = new TFile("/exp/sbnd/app/users/yadav/Calibration/YZ_X_Calib/Median/split_macros/output_files/dQdx_data_dev.root", "READ");
 
  int nplanes = 3;
  
  TH1F *dqdx_nosce[nplanes][6];
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<6;k++){
      dqdx_nosce[l][k] = (TH1F*)file_nosce->Get(Form("dqdxHist_%i_%i",l,k));
    }
  }

  TH1F *dqdx_sce[nplanes][6];
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<6;k++){
      dqdx_sce[l][k] = (TH1F*)file_sce->Get(Form("dqdxHist_%i_%i",l,k));
    }
  }

  const char* plane_label[3] = {"Induction Plane 1", "Induction Plane 2", "Collection Plane"};

  {
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);
    
    SetLineStyle(dqdx_sce[2][0], kBlack);
    SetLineStyle(dqdx_sce[2][1], kAzure+1);
    SetLineStyle(dqdx_nosce[2][2], kGreen+2);
    SetLineStyle(dqdx_sce[2][2], kViolet-2);
    
    //dqdx_sce[2][0]->GetXaxis()->SetRangeUser(600, 2000);
    dqdx_sce[2][0]->GetXaxis()->SetRangeUser(30E3, 100E3);
    SetHist(dqdx_sce[2][0], "", "dQ/dx  [ electrons / cm ]", "number of hits");
    dqdx_sce[2][0]->Draw("hist");
    dqdx_sce[2][1]->Draw("hist same");
    //dqdx_nosce[2][2]->Draw("hist same");
    dqdx_sce[2][2]->Draw("hist same");
    
    if(string(instring).find("data") != string::npos){
      //DrawSBNDLabel("#bf{SBND Preliminary}", 0.7, 0.85, 0.85);
      DrawSBNDLabel("SBND Preliminary", 0.73, 0.84, 0.85);
      DrawLabel("#bf{Data Run 18115}", 0.7, 0.88, 0.95, kBlack, 32);
    } else {
      DrawSBNDLabel("SBND Simulation", 0.6, 0.85, 0.95);
    }
    
    //DrawLabel(plane_label[2], 0.7, 0.85, 0.85, kBlack, 32);
    //DrawLabel("Cathode-Crossing tracks", 0.7, 0.85, 0.8, kBlack, 32);
    
    TLegend *l = new TLegend(0.6, 0.58, 0.85, 0.8);
    l->SetTextSize(0.035);
    l->AddEntry(dqdx_sce[2][0], "No Correction", "lf");
    l->AddEntry(dqdx_sce[2][1], "+ SCE correction", "lf");
    //l->AddEntry(dqdx_nosce[2][2], "+ YZ correction", "lf");
    l->AddEntry(dqdx_sce[2][2], "+ YZ correction", "lf");
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    
    l->Draw();
    c->SaveAs(Form("plot_dir/plot_dqdx/dqdx_%i.pdf", 2));
  }


  {
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);
    
    SetLineStyle(dqdx_sce[2][3], kBlack);
    SetLineStyle(dqdx_sce[2][4], kAzure+1);
    SetLineStyle(dqdx_nosce[2][5], kGreen+2);
    SetLineStyle(dqdx_sce[2][5], kViolet-2);
    
    //dqdx_sce[2][3]->GetXaxis()->SetRangeUser(600, 2000);
    dqdx_sce[2][3]->GetXaxis()->SetRangeUser(30E3, 100E3);
    SetHist(dqdx_sce[2][3], "", "dQ/dx  [ electrons / cm ]", "Ratio to nominal");
    dqdx_sce[2][3]->Draw("hist");
    dqdx_sce[2][4]->Draw("hist same");
    //dqdx_nosce[2][5]->Draw("hist same");
    dqdx_sce[2][5]->Draw("hist same");
    
    if(string(instring).find("data") != string::npos){
      //DrawSBNDLabel("SBND Preliminary", 0.6, 0.85, 0.95);
      DrawLabel("#bf{Data Run 18115}", 0.7, 0.88, 0.95, kBlack, 32);
    } else {
      DrawSBNDLabel("SBND Simulation", 0.6, 0.85, 0.95);
    }
    
    TLegend *l = new TLegend(0.6, 0.65, 0.85, 0.85);
    l->SetTextSize(0.035);
    l->AddEntry(dqdx_sce[2][3], "No Correction", "lf");
    l->AddEntry(dqdx_sce[2][4], "+ SCE correction", "lf");
    //l->AddEntry(dqdx_nosce[2][5], "+ YZ equalization", "lf");
    l->AddEntry(dqdx_sce[2][5], "+ YZ correction", "lf");
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    
    l->Draw();
    c->SaveAs(Form("plot_dir/plot_dqdx/ratio_%i.pdf", 2));
  }
  
 
}
