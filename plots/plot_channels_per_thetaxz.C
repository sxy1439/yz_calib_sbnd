#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>

void SetHist(TH1F *h){

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

  h->GetXaxis()->SetTitle("Channel Number");
  h->GetYaxis()->SetTitle("Median dQ/dx [ADC/cm]");

}

void SetAxis(TH1F *h, double xmin, double xmax, double ymin, double ymax){

  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetRangeUser(ymin, ymax);

}

// median dqdx wrt channel numbers

void plot_channels_per_thetaxz() {
  
  TFile *file_mcp2024b = new TFile("output_files/channel_per_thetaxz_mc2024B_sub123.root", "READ");
  TFile *file_mcp2025av3 = new TFile("output_files/channel_per_thetaxz_mcp2025av3.root", "READ");
  TFile *file_data_dev = new TFile("output_files/channel_per_thetaxz_data_dev.root", "READ");
  TFile *file_data_dec = new TFile("output_files/channel_per_thetaxz_data_runs17742_to_87_sub1.root", "READ");
  
  int nplanes = 3, nthetabin = 19;
  
  TH1F *hmcp2024b[nplanes][nthetabin];
  TH1F *hmcp2025av3[nplanes][nthetabin];
  TH1F *hdata_dev[nplanes][nthetabin];
  TH1F *hdata_dec[nplanes][nthetabin];
  
  for(int i=0; i<nplanes; i++){
    for(int j=0; j<nthetabin; j++){
      hmcp2024b[i][j] = (TH1F*)file_mcp2024b->Get(Form("mediandqdx_%i_%i", i, j));
      hmcp2025av3[i][j] = (TH1F*)file_mcp2025av3->Get(Form("mediandqdx_%i_%i", i, j));
      hdata_dev[i][j] = (TH1F*)file_data_dev->Get(Form("mediandqdx_%i_%i", i, j));
      hdata_dec[i][j] = (TH1F*)file_data_dec->Get(Form("mediandqdx_%i_%i", i, j));
    }
  }

  const char* sample[4] = {"MC 2024B", "MC 2025Av3", "Data Dev", "Data December 2024"};
  
  for(int i=0; i<nplanes; i++){
    for(int j=0; j<nthetabin; j++){
      
      //TCanvas *c = new TCanvas("c", "c", 800, 400);
      TCanvas *c = new TCanvas();
      gStyle->SetOptStat(0);
      TLegend *l = new TLegend(0.2, 0.17, 0.5, 0.35);
      if(i==2) l = new TLegend(0.35, 0.17, 0.65, 0.35);
      l->SetFillStyle(0);
      
      hmcp2024b[i][j]->SetLineColor(kBlack);
      hmcp2024b[i][j]->SetLineWidth(2);
      
      hmcp2025av3[i][j]->SetLineColor(kGreen-6);
      hmcp2025av3[i][j]->SetLineWidth(2);
      
      hdata_dev[i][j]->SetLineColor(kViolet-9);
      hdata_dev[i][j]->SetLineWidth(2);
      
      hdata_dec[i][j]->SetLineColor(kOrange+2);
      hdata_dec[i][j]->SetLineWidth(2);
      
      SetHist(hmcp2024b[i][j]);
      SetHist(hmcp2025av3[i][j]);
      SetHist(hdata_dev[i][j]);
      SetHist(hdata_dec[i][j]);

      hmcp2024b[i][j]->Smooth(2);
      hmcp2025av3[i][j]->Smooth(2);
      hdata_dev[i][j]->Smooth(2);
      hdata_dec[i][j]->Smooth(2);
      
      l->AddEntry(hmcp2024b[i][j], sample[0], "l");
      l->AddEntry(hmcp2025av3[i][j], sample[1], "l");
      l->AddEntry(hdata_dev[i][j], sample[2], "l");
      //l->AddEntry(hdata_dec[i][j], sample[3], "l");
      
      
      if(i==0){
	hmcp2024b[i][j]->SetTitle("Induction (Plane 0)");
	
	SetAxis(hmcp2024b[i][j], 5600, 7600, 0, 1800);
	SetAxis(hmcp2025av3[i][j], 5600, 7600, 0, 1800);
	SetAxis(hdata_dev[i][j], 5600, 7600, 0, 1800);
	SetAxis(hdata_dec[i][j], 5600, 7600, 0, 1800);
	
      } else if(i==1){
	
	hmcp2024b[i][j]->SetTitle("Induction (Plane 1)");
	
	SetAxis(hmcp2024b[i][j], 7600, 9600, 0, 1800);
	SetAxis(hmcp2025av3[i][j], 7600, 9600, 0, 1800);
	SetAxis(hdata_dev[i][j], 7600, 9600, 0, 1800);
	SetAxis(hdata_dec[i][j], 7600, 9600, 0, 1800);
	
      } else if(i==2){
	
	hmcp2024b[i][j]->SetTitle("Collection Plane");
	
	SetAxis(hmcp2024b[i][j], 4000, 11500, 0, 1800);
	SetAxis(hmcp2025av3[i][j], 4000, 11500, 0, 1800);
	SetAxis(hdata_dev[i][j], 4000, 11500, 0, 1800);
	SetAxis(hdata_dec[i][j], 4000, 11500, 0, 1800);	
      }
      
      
      hmcp2024b[i][j]->Draw("hist");
      hmcp2025av3[i][j]->Draw("hist same");
      hdata_dev[i][j]->Draw("hist same");
      //hdata_dec[i][j]->Draw("hist same");
      
      l->Draw();
      c->Print(Form("plot_dir/plot_channels_per_thetaxz/mediandqdx_%i_%i.pdf", i, j));
      
    }
  }
  

}
