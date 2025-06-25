#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TLatex.h"

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
  h->GetZaxis()->SetTitleOffset(1.25);

  h->GetXaxis()->SetTitleSize(0.047);
  h->GetYaxis()->SetTitleSize(0.047);
  h->GetZaxis()->SetTitleSize(0.047);
}

void SetAxis(TH2F *h, double xmin, double xmax, double ymin, double ymax){
  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetRangeUser(ymin, ymax);
}

void DrawSBNDLabel(const char* text="text", float tscale=0.8, double xloc=0.9, double yloc=0.95){
  //TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{SBND Preliminary}", tscale));
  TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  prelim->SetTextColor(kBlack);
  prelim->SetTextFont(132);
  prelim->SetNDC(); 
  prelim->SetTextSize(2.0/30.0); 
  prelim->SetTextAlign(32); 
  prelim->Draw();
}

void DrawLabel(const char* text="text", float tscale=0.8, double xloc=0.9, double yloc=0.95, int tcolor=kGray+1, int align=12){
  TLatex* label = new TLatex(xloc, yloc, Form("#scale[%f]{%s}", tscale, text));
  label->SetTextColor(tcolor);
  label->SetTextFont(132); // new times roman
  label->SetNDC(); 
  label->SetTextSize(2.0/30.0); 
  label->SetTextAlign(align); 
  label->Draw();
}


// median dqdx wrt channel numbers

void plot_yz_planes(const char* inYZOut) {

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
	//CzyHist[l][k]->GetZaxis()->SetRangeUser(0.9, 1.1);
	CzyHist[l][k]->GetZaxis()->SetRangeUser(0, 2);
      } else {
	CzyHist[l][k]->GetZaxis()->SetRangeUser(0.9, 1.1);
	//CzyHist[l][k]->GetZaxis()->SetRangeUser(0, 2);
      }
      CzyHist[l][k]->GetZaxis()->SetTitleOffset(1.18);
      if(string(inYZOut).find("data") != string::npos) CzyHist[l][k]->GetZaxis()->SetTitleOffset(0.98); 
      CzyHist[l][k]->Draw("colz");
      
      //DrawSBNDLabel(0.6, 0.85, 0.95);
      if(string(inYZOut).find("data") != string::npos){
	//DrawSBNDLabel("SBND Preliminary", 0.6, 0.85, 0.95);
	//DrawLabel("Data Run 18115", 0.6, 0.1, 0.95, kBlack, 12);
	
	if(l==2){	
	  DrawSBNDLabel("SBND Preliminary", 0.6, 0.9, 0.95);
	} else {
	  DrawSBNDLabel("SBND Data", 0.6, 0.9, 0.95);
	}
	
	DrawLabel("Data Run 18115", 0.6, 0.15, 0.95, kBlack, 12);
      } else {
	DrawSBNDLabel("SBND Simulation", 0.6, 0.9, 0.95);
	//DrawLabel("Simulation", 0.6, 0.1, 0.95, kBlack, 12);
      }
      
      DrawLabel(tpc[k], 0.6, 0.5, 0.95, kBlack, 22);
      /*
      if(k==0){
	DrawLabel(plane_label[l], 0.7, 0.78, 0.85, kWhite, 32);
	DrawLabel("Cathode-Crossing tracks", 0.7, 0.78, 0.8, kWhite, 32);
      } else {
	DrawLabel(plane_label[l], 0.7, 0.12, 0.85, kWhite, 12);
	DrawLabel("Cathode-Crossing tracks", 0.7, 0.12, 0.8, kWhite, 12);
      }
      */

      if(k==0){
	DrawLabel(plane_label[l], 0.7, 0.83, 0.85, kWhite, 32);
	DrawLabel("Cathode-Crossing tracks", 0.7, 0.83, 0.8, kWhite, 32);
      } else {
	DrawLabel(plane_label[l], 0.7, 0.17, 0.85, kWhite, 12);
	DrawLabel("Cathode-Crossing tracks", 0.7, 0.17, 0.8, kWhite, 12);
      }

      c->SetLeftMargin(0.15);
      c->SetRightMargin(0.15);
      c->SetBottomMargin(0.15);
      if(string(inYZOut).find("data") != string::npos) c->Print(Form("plot_dir/plot_yz_planes/corr_factor/sce_data_yz_cf_%i_%i.pdf", l, k));
      else c->Print(Form("plot_dir/plot_yz_planes/corr_factor/sce_yz_cf_%i_%i.pdf", l, k));

    }
  }

  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      TCanvas *c = new TCanvas();
      gStyle->SetOptStat(0);

      SetHist(zyHist[l][k], "", "Z Coordinate [cm]", "Y Coordinate [cm]", "Median dQ/dx [ADC/cm]");
      zyHist[l][k]->GetZaxis()->SetRangeUser(400, 1800);
      //zyHist[l][k]->SetMinimum(400);
      zyHist[l][k]->Draw("colz");

      //DrawSBNDLabel(0.6, 0.85, 0.95);
      if(string(inYZOut).find("data") != string::npos){
	DrawSBNDLabel("SBND Data", 0.6, 0.85, 0.95);
	DrawLabel("Data Run 18115", 0.6, 0.1, 0.95, kBlack, 12);
      } else {
	DrawSBNDLabel("SBND Simulation", 0.6, 0.85, 0.95);
	//DrawLabel("Simulation", 0.6, 0.1, 0.95, kBlack, 12);
      }

      DrawLabel(tpc[k], 0.6, 0.45, 0.95, kBlack, 22);
      if(k==0){
	DrawLabel(plane_label[l], 0.7, 0.78, 0.85, kWhite, 32);
	DrawLabel("Cathode-Crossing tracks", 0.7, 0.78, 0.8, kWhite, 32);
      } else {
	DrawLabel(plane_label[l], 0.7, 0.12, 0.85, kWhite, 12);
	DrawLabel("Cathode-Crossing tracks", 0.7, 0.12, 0.8, kWhite, 12);
      }

      c->SetLeftMargin(0.1);
      c->SetRightMargin(0.2);
      
      if(string(inYZOut).find("data") != string::npos) c->Print(Form("plot_dir/plot_yz_planes/median/sce_data_yz_cf_%i_%i.pdf", l, k));
      else c->Print(Form("plot_dir/plot_yz_planes/median/sce_yz_cf_%i_%i.pdf", l, k));

    }
  }



}

