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

void plot_channels() {
  
  TFile *file_mcp2024b = new TFile("output_files/channel_mc2024B_sub123.root", "READ");
  TFile *file_mcp2025av3 = new TFile("output_files/channel_mcp2025av3.root", "READ");
  TFile *file_data_dev = new TFile("output_files/channel_data_dev.root", "READ");
  TFile *file_data_dec = new TFile("output_files/channel_data_runs17742_to_87_sub1.root", "READ");
  
  int nplanes = 3;
  
  TH1F *hmcp2024b[nplanes];
  TH1F *hmcp2025av3[nplanes];
  TH1F *hdata_dev[nplanes];
  TH1F *hdata_dec[nplanes];
  
  for(int i=0; i<nplanes; i++){
    hmcp2024b[i] = (TH1F*)file_mcp2024b->Get(Form("mediandqdx_%i", i));
    hmcp2025av3[i] = (TH1F*)file_mcp2025av3->Get(Form("mediandqdx_%i", i));
    hdata_dev[i] = (TH1F*)file_data_dev->Get(Form("mediandqdx_%i", i));
    hdata_dec[i] = (TH1F*)file_data_dec->Get(Form("mediandqdx_%i", i));
  }

  const char* sample[4] = {"MC 2024B", "MC 2025Av3", "Data Dev", "Data December 2024"};
  
  for(int i=0; i<nplanes; i++){
    
    //TCanvas *c = new TCanvas("c", "c", 800, 400);
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);
    TLegend *l = new TLegend(0.2, 0.17, 0.5, 0.35);
    if(i==2) l = new TLegend(0.35, 0.17, 0.65, 0.35);
    l->SetFillStyle(0);
    
    hmcp2024b[i]->SetLineColor(kBlack);
    hmcp2024b[i]->SetLineWidth(2);
    
    hmcp2025av3[i]->SetLineColor(kGreen-6);
    hmcp2025av3[i]->SetLineWidth(2);
    
    hdata_dev[i]->SetLineColor(kViolet-9);
    hdata_dev[i]->SetLineWidth(2);
    
    hdata_dec[i]->SetLineColor(kOrange+2);
    hdata_dec[i]->SetLineWidth(2);
    
    SetHist(hmcp2024b[i]);
    SetHist(hmcp2025av3[i]);
    SetHist(hdata_dev[i]);
    SetHist(hdata_dec[i]);
    
    /*
    hmcp2024b[i]->Smooth(2);
    hmcp2025av3[i]->Smooth(2);
    hdata_dev[i]->Smooth(2);
    hdata_dec[i]->Smooth(2);
    */

    l->AddEntry(hmcp2024b[i], sample[0], "l");
    l->AddEntry(hmcp2025av3[i], sample[1], "l");
    l->AddEntry(hdata_dev[i], sample[2], "l");
    l->AddEntry(hdata_dec[i], sample[3], "l");
    
    
    if(i==0){
      hmcp2024b[i]->SetTitle("Induction (Plane 0)");
      
      SetAxis(hmcp2024b[i], 5600, 7600, 0, 1800);
      SetAxis(hmcp2025av3[i], 5600, 7600, 0, 1800);
      SetAxis(hdata_dev[i], 5600, 7600, 0, 1800);
      SetAxis(hdata_dec[i], 5600, 7600, 0, 1800);
      
    } else if(i==1){
      //l->AddEntry(hmcp2024b[i], "Induction Plane 2", "l");

      hmcp2024b[i]->SetTitle("Induction (Plane 1)");
      
      SetAxis(hmcp2024b[i], 7600, 9600, 0, 1800);
      SetAxis(hmcp2025av3[i], 7600, 9600, 0, 1800);
      SetAxis(hdata_dev[i], 7600, 9600, 0, 1800);
      SetAxis(hdata_dec[i], 7600, 9600, 0, 1800);
      
    } else if(i==2){
      //l->AddEntry(hmcp2024b[i], "Collection Plane", "l");

      hmcp2024b[i]->SetTitle("Collection Plane");
      
      SetAxis(hmcp2024b[i], 4000, 11500, 0, 1800);
      SetAxis(hmcp2025av3[i], 4000, 11500, 0, 1800);
      SetAxis(hdata_dev[i], 4000, 11500, 0, 1800);
      SetAxis(hdata_dec[i], 4000, 11500, 0, 1800);	
    }

    
    hmcp2024b[i]->Draw("hist");
    hmcp2025av3[i]->Draw("hist same");
    hdata_dev[i]->Draw("hist same");
    hdata_dec[i]->Draw("hist same");
    
    l->Draw();
    c->Print(Form("plot_dir/plot_channels/mediandqdx_%i.pdf", i));
    
  }
  

}
