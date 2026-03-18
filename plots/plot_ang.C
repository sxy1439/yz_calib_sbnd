#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TLatex.h"
#include "TLine.h"

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

void DrawSBNDLabel(float tscale=0.8, double xloc=0.9, double yloc=0.95){
  TLatex* prelim = new TLatex(xloc, yloc, Form("#scale[%f]{SBND Preliminary}", tscale));
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

void DrawLine(double x1, double y1, double x2, double y2, int color = kBlack, int style = 1, int width = 2){
  TLine* line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(color);
  line->SetLineStyle(style); // 1 = solid, 2 = dashed, etc.
  line->SetLineWidth(width);
  line->Draw("same");
}



//void plot_ang(const char* inFile, const char* data_run) {
void plot_ang(const char* inFile) {
  
  TFile* file = new TFile(inFile, "READ");
  
  TH2F *thetaHistEast = (TH2F*)file->Get("thetaHistNeg");
  TH2F *thetaHistWest = (TH2F*)file->Get("thetaHistPos");

  const char* tpc[2] = {"East TPC", "West TPC"};
  double maxz = std::max(thetaHistEast->GetMaximum(), thetaHistWest->GetMaximum());


  {
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);
    
    SetHist(thetaHistEast, "", "#theta_{xz} [deg.]", "#theta_{yz} [deg.]", "Avg. dQ/dx [ADC/cm]");
    thetaHistEast->GetZaxis()->SetRangeUser(0, maxz);
    sbndstyle::SetSBNDStyle();
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    
    thetaHistEast->SetStats(0);
    
    thetaHistEast->Draw("colz");

    DrawLabel(tpc[0], 0.7, 0.15, 0.95, kBlack, 12);
    if(string(inFile).find("data") != string::npos){
      //DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
      //DrawLabel("SBND Data, Preliminary", 0.7, 0.85, 0.95, kBlack, 32);
      DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
    } else {
      DrawLabel("SBND Simulation, Preliminary", 0.7, 0.85, 0.95, kBlack, 32);
    }

    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), 110, thetaHistEast->GetXaxis()->GetXmax(), 112, kBlack, 2, 2);
    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), -110, thetaHistEast->GetXaxis()->GetXmax(), -112, kBlack, 2, 2);
    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), 70, thetaHistEast->GetXaxis()->GetXmax(), 72, kBlack, 2, 2);
    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), -70, thetaHistEast->GetXaxis()->GetXmax(), -72, kBlack, 2, 2);

    DrawLine(115, thetaHistEast->GetYaxis()->GetXmin(), 117, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(-115, thetaHistEast->GetYaxis()->GetXmin(), -117, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(65, thetaHistEast->GetYaxis()->GetXmin(), 67, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(-65, thetaHistEast->GetYaxis()->GetXmin(), -67, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.19);
    if(string(inFile).find("data") != string::npos) c->Print("plot_dir/plot_ang/east_ang_data_test.pdf");
    else c->Print("plot_dir/plot_ang/east_ang_test.pdf");
  }
  

  {
    TCanvas *c = new TCanvas();
   
    SetHist(thetaHistEast, "", "#theta_{xz} [deg.]", "#theta_{yz} [deg.]", "Avg. dQ/dx [ADC/cm]");
    thetaHistEast->GetZaxis()->SetRangeUser(0, maxz);
    sbndstyle::SetSBNDStyle();
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    
    thetaHistEast->SetStats(0);
    
    thetaHistEast->Draw("colz");
    //sbndstyle::Preliminary();

    DrawLabel(tpc[0], 0.7, 0.15, 0.95, kBlack, 12);
    
    if(string(inFile).find("data") != string::npos){
      //DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
      //DrawLabel("SBND Data, Preliminary", 0.7, 0.85, 0.95, kBlack, 32);
      DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
    } else {
      DrawLabel("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
    }

    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), 110, thetaHistEast->GetXaxis()->GetXmax(), 112, kBlack, 2, 2);
    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), -110, thetaHistEast->GetXaxis()->GetXmax(), -112, kBlack, 2, 2);
    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), 70, thetaHistEast->GetXaxis()->GetXmax(), 72, kBlack, 2, 2);
    DrawLine(thetaHistEast->GetXaxis()->GetXmin(), -70, thetaHistEast->GetXaxis()->GetXmax(), -72, kBlack, 2, 2);

    DrawLine(115, thetaHistEast->GetYaxis()->GetXmin(), 117, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(-115, thetaHistEast->GetYaxis()->GetXmin(), -117, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(65, thetaHistEast->GetYaxis()->GetXmin(), 67, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(-65, thetaHistEast->GetYaxis()->GetXmin(), -67, thetaHistEast->GetYaxis()->GetXmax(), kBlack, 2, 2);

    sbndstyle::SetSBNDStyle();

    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.19);
    if(string(inFile).find("data") != string::npos) c->Print("plot_dir/plot_ang/east_ang_data.pdf");
    else c->Print("plot_dir/plot_ang/east_ang.pdf");
  }
  
  
  {
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);
    
    SetHist(thetaHistWest, "", "#theta_{xz} [deg.]", "#theta_{yz} [deg.]", "Avg. dQ/dx [ADC/cm]");
    thetaHistWest->GetZaxis()->SetRangeUser(0, maxz);
    sbndstyle::SetSBNDStyle();
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    
    thetaHistWest->SetStats(0);
    
    thetaHistWest->Draw("colz");

    DrawLabel(tpc[1], 0.7, 0.15, 0.95, kBlack, 12);
    if(string(inFile).find("data") != string::npos){
      //DrawLabel(Form("SBND Data Run %s", data_run), 0.7, 0.85, 0.95, kBlack, 32);
      //DrawLabel("SBND Data, Preliminary", 0.7, 0.85, 0.95, kBlack, 32);
      DrawLabel("SBND Data", 0.7, 0.85, 0.95, kBlack, 32);
    } else {
      DrawLabel("SBND Simulation", 0.7, 0.85, 0.95, kBlack, 32);
    }

    DrawLine(thetaHistWest->GetXaxis()->GetXmin(), 110, thetaHistWest->GetXaxis()->GetXmax(), 112, kBlack, 2, 2);
    DrawLine(thetaHistWest->GetXaxis()->GetXmin(), -110, thetaHistWest->GetXaxis()->GetXmax(), -112, kBlack, 2, 2);
    DrawLine(thetaHistWest->GetXaxis()->GetXmin(), 70, thetaHistWest->GetXaxis()->GetXmax(), 72, kBlack, 2, 2);
    DrawLine(thetaHistWest->GetXaxis()->GetXmin(), -70, thetaHistWest->GetXaxis()->GetXmax(), -72, kBlack, 2, 2);

    DrawLine(115, thetaHistWest->GetYaxis()->GetXmin(), 117, thetaHistWest->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(-115, thetaHistWest->GetYaxis()->GetXmin(), -117, thetaHistWest->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(65, thetaHistWest->GetYaxis()->GetXmin(), 67, thetaHistWest->GetYaxis()->GetXmax(), kBlack, 2, 2);
    DrawLine(-65, thetaHistWest->GetYaxis()->GetXmin(), -67, thetaHistWest->GetYaxis()->GetXmax(), kBlack, 2, 2);
    
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.19);
    if(string(inFile).find("data") != string::npos) c->Print("plot_dir/plot_ang/west_ang_data.pdf");
    else c->Print("plot_dir/plot_ang/west_ang.pdf");
  }

  


}

