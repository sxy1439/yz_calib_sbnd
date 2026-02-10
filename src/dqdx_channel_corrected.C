#include "../include/utilities.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>

// to run :
//root [0] .L dqdx_channel_corrected.C
//root [1] dqdx_channel_corrected("doData_runs17742_to_87_ctcMap_sub1", "output_files/channel_data_runs17742_to_87_ctcMap_sub1.root")


void dqdx_channel_corrected(const char* cintyp, const char* inChanOut) {

  TFile* file_ChanOut = new TFile(inChanOut, "READ"); 
  TH2F* CcHist = (TH2F*)file_ChanOut->Get("CcHist");
  
  gSystem->Load("libFileHandling.so");
  if (gSystem->Load("libFileHandling.so") < 0) {
    std::cerr << "Error loading library" << std::endl;
    return;
  }

  auto e_callife = getELifetime(cintyp);
  
  TString output_name = getOutputNamedQdx(cintyp);
  TFile *out_file = new TFile(output_name, "recreate");
  
  // Initialize histograms
  initialize_histograms();
  
  // Open the file containing the tree.
  TChain *tc = new TChain("caloskim/TrackCaloSkim");
  std::ifstream fl = openInputFile(cintyp);
  if (!fl) {
    std::cerr << "Error opening input file.\n";
    return;
  }
  
  TString filename;
  while(fl.peek()!=EOF){
    fl>>filename;
    tc->AddFile(filename);
  }
  
  TTreeReader myReader(tc);

  //  Variables we are going to read
  TTreeReaderValue<int> selected(myReader, "trk.selected");
  TTreeReaderValue<bool> is_mu(myReader, "trk.clear_cosmic_muon");
  TTreeReaderValue<float> startx(myReader, "trk.start.x");
  TTreeReaderValue<float> starty(myReader, "trk.start.y");
  TTreeReaderValue<float> startz(myReader, "trk.start.z");
  TTreeReaderValue<float> endx(myReader, "trk.end.x");
  TTreeReaderValue<float> endy(myReader, "trk.end.y");
  TTreeReaderValue<float> endz(myReader, "trk.end.z");
  TTreeReaderValue<float> dirx(myReader, "trk.dir.x");
  TTreeReaderValue<float> diry(myReader, "trk.dir.y");
  TTreeReaderValue<float> dirz(myReader, "trk.dir.z");
  TTreeReaderValue<float> t0(myReader, "trk.t0");
  float thetaxz, thetayz;

  TTreeReaderArray<float> dqdx(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> tpx(myReader, "trk.hits2.h.sp.x"); // x of track trajectory position  (older was -> trk.hits2.tp.x)
  TTreeReaderArray<float> tpy(myReader, "trk.hits2.h.sp.y"); // y of track trajectory position  
  TTreeReaderArray<float> tpz(myReader, "trk.hits2.h.sp.z"); // z of track trajectory position
  TTreeReaderArray<float> time(myReader, "trk.hits2.h.time"); // in ticks (500 ns), up to 3200 (1.6 ms > 1.25 ms)
  TTreeReaderArray<UShort_t> channel(myReader, "trk.hits2.h.channel");
  

  cout<<"nbinx = "<<nbinx<<" nbiny = "<<nbiny<<" nbinz = "<<nbinz<<endl;
  cout<<"elife[0] = "<<e_callife[0]<<" elife[1] = "<<e_callife[0]<<endl;

  /// FINAL dQdx VALUES (dQ/dx with and without equalization correction)

  int ibinx, ibiny, ibinz, ibinc;
  myReader.Restart();
  while (myReader.Next()) {
    if (!*is_mu) continue; //Pandora clear muon
    if ( !Is_Edge(*startx, *starty, *startz) && !Is_Edge(*endx, *endy, *endz)) continue;//FV
    if(!Is_Cathode_Crossing(*startx, *endx)) continue;
    
    thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi(); 
    if(*dirx<0) thetaxz = -thetaxz;
    thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
    if(*diry<0) thetayz = -thetayz;
    
    if(abs(thetaxz)<115&&abs(thetaxz)>65)continue;//Angle
    if(abs(thetayz)<110&&abs(thetayz)>70)continue;//Angle

    float trk_t0 = *t0;
    
    for (unsigned i = 0; i < dqdx.GetSize(); i++) {
      if(isnan(tpx[i])||isnan(tpy[i])||isnan(tpz[i]))continue;
      if(isnan(dqdx[i]) || isinf(dqdx[i])) continue;

      // masked YZ and X regions    
      if(string(cintyp) == "doData_runs17742_to_87_sub1" || string(cintyp) == "doData_runs17742_to_87_ctcMap_sub1"){
	if(InVeto_region_X(tpx[i])) continue;
	if(tpx[i]<0){
	  if(InVeto_region_eastTPC_C(tpy[i], tpz[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy[i], tpz[i])) continue;
	}
      }

      ibinc = std::round((channel[i]-lowc)/(highc-lowc)*nbinc);
      
      double elife0 = lifetime_correction(time[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      double CF_cg = CcHist->GetBinContent(ibinc+1);  // channel-by-channel gain

      if(tpx[i]<0){
	dqdxHist[0]->Fill(dqdx[i] * elife0);
	dqdxHist[1]->Fill(dqdx[i] * elife0 * CF_cg);
      }
      else{
	dqdxHist[0]->Fill(dqdx[i] * elife1);
	dqdxHist[1]->Fill(dqdx[i] * elife1 * CF_cg);
      }
    }
  }

  
  for(int k=0;k<2;k++)dqdxHist[k]->Write();

  //TGaxis::SetMaxDigits(3);
  TCanvas *c1 = new TCanvas();

  dqdxHist[0]->SetLineColor(kBlack);
  dqdxHist[0]->SetLineWidth(2);
  dqdxHist[0]->SetStats(0);
  
  dqdxHist[1]->SetLineColor(kViolet+2);
  dqdxHist[1]->SetLineWidth(2);
  dqdxHist[1]->SetStats(0);

  /*
  dqdxHist[2]->SetLineColor(kGreen+2);  // kGreen+3 (3365)
  dqdxHist[2]->SetLineWidth(2);
  dqdxHist[2]->SetStats(0);
  */

  dqdxHist[0]->Draw("hist");
  dqdxHist[1]->Draw("hist same");
  //dqdxHist[2]->Draw("hist same");
 
  cout<< "Uncalibrated dQ/dx "<<" Mean : "<< dqdxHist[0]->GetMean()<<" SD : "<<dqdxHist[0]->GetStdDev()<<endl;
  cout<< "Channel Gain CF applied dQ/dx "<<" Mean : "<< dqdxHist[1]->GetMean()<<" SD : "<<dqdxHist[1]->GetStdDev()<<endl;
  //cout<< "YZ CF applied dQ/dx "<<" Mean : "<< dqdxHist[1]->GetMean()<<" SD : "<<dqdxHist[1]->GetStdDev()<<endl;
  //cout<< "X CF applied dQ/dx "<<" Mean : "<< dqdxHist[2]->GetMean()<<" SD : "<<dqdxHist[2]->GetStdDev()<<endl;

 
  TLegend *l1 = new TLegend(0.6, 0.7, 0.85, 0.85);
  l1->SetTextSize(0.030);
  l1->AddEntry(dqdxHist[0], "No Correction", "f");
  l1->AddEntry(dqdxHist[1], "Wire Plane Equalization", "f");
  //l1->AddEntry(dqdxHist[2], "Wire Plane + Drift Equalization", "f");
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);

  l1->Draw();
  c1->SetLeftMargin(0.15);
  c1->SaveAs(Form("plots_dqdxHists/finaldqdx_%s.pdf", cintyp));

  out_file->Close();
  
}
