#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TLatex.h"

#define SBNDSTYLE_ENABLE_AUTOMATICALLY 0
#include "SBNDStyle.h"

void SetHist(TH2F *h, const char* title, const char* xlabel, const char* ylabel, const char* zlabel){
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

  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->GetZaxis()->SetLabelSize(.05);
}

void SetAxis(TH2F *h, double xmin, double xmax, double ymin, double ymax){
  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetRangeUser(ymin, ymax);
}

void DrawSBNDLabel(const char* text="text", float tscale=0.8, double xloc=0.9, double yloc=0.95){
  //TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{SBND Preliminary}", tscale));
  TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  prelim->SetTextColor(kBlack);
  prelim->SetTextFont(42);
  prelim->SetNDC(); 
  prelim->SetTextSize(2.0/30.0); 
  prelim->SetTextAlign(32); 
  prelim->Draw();
}

void DrawLabel(const char* text="text", float tscale=0.8, double xloc=0.9, double yloc=0.95, int tcolor=kGray+1, int align=12){
  TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  label->SetTextColor(tcolor);
  label->SetTextFont(42); // new times roman
  label->SetNDC(); 
  label->SetTextSize(2.0/30.0); 
  label->SetTextAlign(align); 
  label->Draw();
}


// median dqdx wrt channel numbers

//void plot_yz_planes(const char* inYZOut, const char* data_run) {
void plot_yz_planes(const char* inYZOut, bool HasSCE=false) {

  int nplanes = 3;
  
  TH2F *zyHist[nplanes][2];
  TH2F *CzyHist[nplanes][2];
  TFile* file_YZOut = new TFile(inYZOut, "READ");
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      zyHist[l][k] = (TH2F*)file_YZOut->Get(Form("zy_%i_%i",l,k));
      CzyHist[l][k] = (TH2F*)file_YZOut->Get(Form("CzyHist_%i_%i",l,k));
    }
  }

  const char* plane_label[3] = {"Induction Plane 1", "Induction Plane 2", "Collection Plane"};
  const char* tpc[2] = {"East TPC", "West TPC"};

  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      TCanvas *c = new TCanvas();
      gStyle->SetOptStat(0);

      SetHist(CzyHist[l][k], "", "Z Coordinate [cm]", "Y Coordinate [cm]", "YZ correction factor");
      if(string(inYZOut).find("data") != string::npos){
	CzyHist[l][k]->GetZaxis()->SetRangeUser(0.8, 1.12);
      } else {
	CzyHist[l][k]->GetZaxis()->SetRangeUser(0.9, 1.12);
      }

      sbndstyle::SetSBNDStyle();
      gROOT->ForceStyle();
      gStyle->SetOptStat(0);
      
      CzyHist[l][k]->SetStats(0);

      sbndstyle::SeaPalette();
      //sbndstyle::SymmetricPalette();
      
      CzyHist[l][k]->Draw("colz");
      

      DrawLabel(tpc[k], 0.7, 0.15, 0.95, kBlack, 12);
      
      if(string(inYZOut).find("data") != string::npos){
	//DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
	DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
      } else {
	DrawLabel("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
      }

      c->SetLeftMargin(0.15);
      c->SetBottomMargin(0.15);
      c->SetRightMargin(0.19);
      if(HasSCE){
	if(string(inYZOut).find("data") != string::npos) c->Print(Form("plot_dir/plot_yz_planes/corr_factor/sce_data_yz_cf_%i_%i.pdf", l, k));
	else c->Print(Form("plot_dir/plot_yz_planes/corr_factor/sce_yz_cf_%i_%i.pdf", l, k));
      } else {
	if(string(inYZOut).find("data") != string::npos) c->Print(Form("plot_dir/plot_yz_planes/corr_factor/nosce_data_yz_cf_%i_%i.pdf", l, k));
	else c->Print(Form("plot_dir/plot_yz_planes/corr_factor/nosce_yz_cf_%i_%i.pdf", l, k));
      }
      
    }
  }


  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      TCanvas *c = new TCanvas();
      gStyle->SetOptStat(0);

      SetHist(zyHist[l][k], "", "Z Coordinate [cm]", "Y Coordinate [cm]", "Median dQ/dx [ADC/cm]");
      if(string(inYZOut).find("data") != string::npos){
	zyHist[l][k]->GetZaxis()->SetRangeUser(500, 1300);
      } else {
	zyHist[l][k]->GetZaxis()->SetRangeUser(500, 1300);
      }

      sbndstyle::SetSBNDStyle();
      gROOT->ForceStyle();
      gStyle->SetOptStat(0);
      
      zyHist[l][k]->SetStats(0);

      sbndstyle::SeaPalette();
      //sbndstyle::SymmetricPalette();
      
      zyHist[l][k]->Draw("colz");
      

      DrawLabel(tpc[k], 0.7, 0.15, 0.95, kBlack, 12);
      
      if(string(inYZOut).find("data") != string::npos){
	//DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
	DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
      } else {
	DrawLabel("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
      }

      c->SetLeftMargin(0.15);
      c->SetBottomMargin(0.15);
      c->SetRightMargin(0.19);
      if(HasSCE){
	if(string(inYZOut).find("data") != string::npos) c->Print(Form("plot_dir/plot_yz_planes/median/sce_data_yz_%i_%i.pdf", l, k));
	else c->Print(Form("plot_dir/plot_yz_planes/median/sce_yz_%i_%i.pdf", l, k));
      } else {
	if(string(inYZOut).find("data") != string::npos) c->Print(Form("plot_dir/plot_yz_planes/median/nosce_data_yz_%i_%i.pdf", l, k));
	else c->Print(Form("plot_dir/median/median/nosce_yz_%i_%i.pdf", l, k));
      }

    }
  }




}

