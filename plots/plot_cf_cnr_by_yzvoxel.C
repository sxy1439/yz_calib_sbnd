#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TLatex.h"

#define SBNDSTYLE_ENABLE_AUTOMATICALLY 0
#include "SBNDStyle.h"

void SetHist(TH1F *h, const char* title, const char* xlabel, const char* ylabel){
  h->SetTitle(title);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.1);

  h->GetXaxis()->SetTitleSize(0.055);
  h->GetYaxis()->SetTitleSize(0.055);

  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
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

TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}




void plot_cf_cnr_by_yzvoxel(const char* instring){

  TFile *file_withcnr = new TFile("/exp/sbnd/app/users/yadav/Calibration/YZ_X_Calib/Median/split_macros/output_files/calib_paper_output/yz_data1e20.root", "READ");
  TFile *file_nocnr = new TFile("/exp/sbnd/app/users/yadav/Calibration/YZ_X_Calib/Median/split_macros/output_files/calib_paper_output/sce/yz_data_fallrun1dev.root", "READ");

  int nplanes = 3;
  int ntpc = 2;
  
  TH1F *histyzcf_withcnr[nplanes][ntpc];
  TH1F *histyzcf_nocnr[nplanes][ntpc];
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      histyzcf_withcnr[l][k] = (TH1F*)file_withcnr->Get(Form("histyzvox_cf_%i_%i",l,k));
      histyzcf_nocnr[l][k] = (TH1F*)file_nocnr->Get(Form("histyzvox_cf_%i_%i",l,k));
    }
  }

  const char* plane_label[3] = {"Induction Plane 1", "Induction Plane 2", "Collection Plane"};
  const char* tpc[2] = {"East TPC", "West TPC"};

  
  {
    TCanvas *c = new TCanvas();
    c->Clear();
    TLegend * leg = MakeLegend(0.6, 0.72, 0.82, 0.85);
    leg->Clear();
    
    gStyle->SetOptStat(0);

    for(int l=0;l<nplanes;l++){
      for(int k=0;k<2;k++){
	SetLineStyle(histyzcf_withcnr[l][k], sbndstyle::colors::kOkabeItoOrange);
	SetLineStyle(histyzcf_nocnr[l][k], sbndstyle::colors::kOkabeItoBlueGreen);
	
	double maxy = std::max(histyzcf_withcnr[l][k]->GetMaximum(), histyzcf_nocnr[l][k]->GetMaximum());
	histyzcf_withcnr[l][k]->GetXaxis()->SetRangeUser(0.9, 1.15);
	histyzcf_withcnr[l][k]->GetYaxis()->SetRangeUser(0, 1.05*maxy);
	SetHist(histyzcf_withcnr[l][k], "", "Correction Factor", "Number of YZ voxels");
	
	sbndstyle::SetSBNDStyle();
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	
	histyzcf_withcnr[l][k]->SetStats(0);
	histyzcf_nocnr[l][k]->SetStats(0);
	
	histyzcf_withcnr[l][k]->Draw("hist");
	histyzcf_nocnr[l][k]->Draw("hist same");
	
	leg->Clear();
	leg->AddEntry(histyzcf_withcnr[l][k],"with CNR","lf");
	leg->AddEntry(histyzcf_nocnr[l][k],"without CNR","lf");
	
	// sbndstyle::colors::kOkabeItoBlue
	
	DrawLabel(tpc[k], 0.7, 0.15, 0.95, kBlack, 12);
	
	if(string(instring).find("data") != string::npos){
	  //DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
	  DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
	} else {
	  DrawLabel("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
	}
	
	leg->Draw();
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);
	//c->SetRightMargin(0.19);
	c->SaveAs(Form("plot_dir/plot_cf_cnr_by_yzvoxel/test_cf_per_yzvoxel_%i_%i.pdf", l, k));
      }
    }
  }


  
  {
    TCanvas *c = new TCanvas();
    c->Clear();
    TLegend * leg = MakeLegend(0.58, 0.67, 0.8, 0.85);
    leg->SetTextSize(0.05);
    leg->Clear();
    
    gStyle->SetOptStat(0);

    for(int l=0;l<nplanes;l++){
      for(int k=0;k<2;k++){
	SetLineStyle(histyzcf_withcnr[l][k], sbndstyle::colors::kOkabeItoOrange);
	SetLineStyle(histyzcf_nocnr[l][k], sbndstyle::colors::kOkabeItoBlueGreen);
	
	double maxy = std::max(histyzcf_withcnr[l][k]->GetMaximum(), histyzcf_nocnr[l][k]->GetMaximum());
	histyzcf_withcnr[l][k]->GetXaxis()->SetRangeUser(0.9, 1.15);
	histyzcf_withcnr[l][k]->GetYaxis()->SetRangeUser(0, 1.05*maxy);
	SetHist(histyzcf_withcnr[l][k], "", "Correction Factor", "Number of YZ voxels");
	
	sbndstyle::SetSBNDStyle();
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	
	histyzcf_withcnr[l][k]->SetStats(0);
	histyzcf_nocnr[l][k]->SetStats(0);
	
	histyzcf_withcnr[l][k]->Draw("hist");
	histyzcf_nocnr[l][k]->Draw("hist same");
	
	leg->Clear();
	leg->AddEntry(histyzcf_withcnr[l][k],"with CNR","lf");
	leg->AddEntry(histyzcf_nocnr[l][k],"without CNR","lf");
	
	// sbndstyle::colors::kOkabeItoBlue

	DrawLabel(tpc[k], 0.7, 0.15, 0.95, kBlack, 12);
	
	if(string(instring).find("data") != string::npos){
	  //DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
	  DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
	} else {
	  DrawLabel("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
	}
	
	leg->Draw();
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);
	//c->SetRightMargin(0.19);
	c->SaveAs(Form("plot_dir/plot_cf_cnr_by_yzvoxel/cf_per_yzvoxel_%i_%i.pdf", l, k));
      }
    }
  }


  

 
 
}
