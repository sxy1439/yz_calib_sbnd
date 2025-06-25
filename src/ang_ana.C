
#include "../include/utilities.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>


// to run : root -l -b -q 'ang_ana.C("doMC2023B_sub1")' > mpv_out_mc.txt


void ang_ana(const char* cintyp) {

  gSystem->Load("libFileHandling.so");
  if (gSystem->Load("libFileHandling.so") < 0) {
    std::cerr << "Error loading library" << std::endl;
    return;
  }

  auto e_callife = getELifetime(cintyp);
 
  TString output_name = getOutputNameAng(cintyp);
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
  TTreeReaderArray<float> time(myReader, "trk.hits2.h.time"); // in ticks (500 ns), up to 3200 (1.6 ms > 1.25 ms) (the comment in TrackCaloSkimmerObj.h does not look right)


  cout<<"nbiny = "<<nbiny<<" nbinz = "<<nbinz<<endl;

  
  /// ANGULAR REQUIREMENT STUDY

  // Loop over all entries of the TTree
  int ntrkpre=0;
  int ntrkpospre=0, ntrknegpre=0;

  while (myReader.Next()) {  
    if (!*is_mu) continue; 
    if ( !Is_Edge(*startx, *starty, *startz) || !Is_Edge(*endx, *endy, *endz)) continue;
    if(!Is_Cathode_Crossing(*startx, *endx)) continue;
    
    thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi();   
    if(*dirx<0) thetaxz = -thetaxz;
    
    thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
    if(*diry<0) thetayz = -thetayz;
    
    thetaxyHist->Fill(thetaxz, thetayz);
    
    if (ntrkpre % 1000 == 0 && ntrkpre > 0) {
      std::cout << "entry => " << ntrkpre / 1000 << "k tracks" << std::endl;
    }
    
    ntrkpre++;
    
    if (*startx >= 0){
      ntrkpospre++;
    } else {
      ntrknegpre++;
    }
    
    float trk_t0 = *t0;
    
    for(unsigned i = 0; i < dqdx.GetSize(); i++){
      if (isnan(tpx[i]) || isnan(tpy[i]) || isnan(tpz[i])) continue;
      if(isnan(dqdx[i]) || isinf(dqdx[i])) continue;
      
      // masked YZ and X regions
      if(string(cintyp) == "doData_runs17742_to_87_sub1" || string(cintyp) == "doData_runs17742_to_87_ctcMap_sub1"){
	if(tpx[i]<0){
	  if(InVeto_region_eastTPC_C(tpy[i], tpz[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy[i], tpz[i])) continue;
	}
      }
      
      if (tpx[i]<0) {
	// for 2D theta plot
	thetaHistNeg->Fill(thetaxz, thetayz, dqdx[i]*lifetime_correction(time[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]));   
	// for the number of dqdx entries per bin to get the average of dqdx
	countHistNeg->Fill(thetaxz, thetayz, 1);
      } else {
	thetaHistPos->Fill(thetaxz, thetayz, dqdx[i]*lifetime_correction(time[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]));        
	countHistPos->Fill(thetaxz, thetayz, 1); 	  
      }
    }
    
  }
  
  for(int ix = 0; ix < thetaHistPos->GetNbinsX(); ix++){
    for(int iy = 0; iy < thetaHistPos->GetNbinsY(); iy++){
      
      float binContentpos = thetaHistPos->GetBinContent(ix+1, iy+1);
      float binEntriespos = countHistPos->GetBinContent(ix+1, iy+1);
      if(binEntriespos != 0){
	thetaHistPos->SetBinContent(ix+1, iy+1, binContentpos / binEntriespos);
      } else {
	thetaHistPos->SetBinContent(ix+1, iy+1, 0.001);
      }
      
      float binContentneg = thetaHistNeg->GetBinContent(ix+1, iy+1);
      float binEntriesneg = countHistNeg->GetBinContent(ix+1, iy+1);
      if(binEntriesneg != 0){
	thetaHistNeg->SetBinContent(ix+1, iy+1, binContentneg / binEntriesneg);
      } else {
	thetaHistNeg->SetBinContent(ix+1, iy+1, 0.001);
      }
      
    }
  }
  
  
  cout<<"Total tracks (before angular cut): "<<ntrkpre<<endl;
  cout<<"# tracks in x>=0 region (before angular cut): "<<ntrkpospre<<endl;
  cout<<"# tracks in x<0 region (before angular cut): "<<ntrknegpre<<endl;
  
  double maxZth = std::max(thetaHistPos->GetMaximum(), thetaHistNeg->GetMaximum());
  
  cout<<"theta maxZth : "<<maxZth<<endl;

  TH2F *thetaHistPos_zfix = (TH2F*)thetaHistPos->Clone("thetaHistPos_zfix");
  TH2F *thetaHistNeg_zfix = (TH2F*)thetaHistPos->Clone("thetaHistNeg_zfix");

  thetaxyHist->GetXaxis()->CenterTitle();
  thetaxyHist->GetYaxis()->CenterTitle();
  thetaxyHist->GetZaxis()->CenterTitle();
  thetaHistPos->GetXaxis()->CenterTitle();
  thetaHistPos->GetYaxis()->CenterTitle();
  thetaHistPos->GetZaxis()->CenterTitle();
  thetaHistNeg->GetXaxis()->CenterTitle();
  thetaHistNeg->GetYaxis()->CenterTitle();
  thetaHistNeg->GetZaxis()->CenterTitle();
  thetaHistPos_zfix->GetXaxis()->CenterTitle();
  thetaHistPos_zfix->GetYaxis()->CenterTitle();
  thetaHistPos_zfix->GetZaxis()->CenterTitle();
  thetaHistNeg_zfix->GetXaxis()->CenterTitle();
  thetaHistNeg_zfix->GetYaxis()->CenterTitle();
  thetaHistNeg_zfix->GetZaxis()->CenterTitle();

  thetaxyHist->SetXTitle("#theta_{xz} [deg]");
  thetaxyHist->SetYTitle("#theta_{yz} [deg]");
  thetaxyHist->SetZTitle("Avg. dQ/dx [ADC/cm]");
  thetaHistPos->SetXTitle("#theta_{xz} [deg]");
  thetaHistPos->SetYTitle("#theta_{yz} [deg]");
  thetaHistPos->SetZTitle("Avg. dQ/dx [ADC/cm]");
  thetaHistNeg->SetXTitle("#theta_{xz} [deg]");
  thetaHistNeg->SetYTitle("#theta_{yz} [deg]");
  thetaHistNeg->SetZTitle("Avg. dQ/dx [ADC/cm]");
  thetaHistPos_zfix->SetXTitle("#theta_{xz} [deg]");
  thetaHistPos_zfix->SetYTitle("#theta_{yz} [deg]");
  thetaHistPos_zfix->SetZTitle("Avg. dQ/dx [ADC/cm]");
  thetaHistNeg_zfix->SetXTitle("#theta_{xz} [deg]");
  thetaHistNeg_zfix->SetYTitle("#theta_{yz} [deg]");
  thetaHistNeg_zfix->SetZTitle("Avg. dQ/dx [ADC/cm]");
  
  
  thetaxyHist->Write();
  thetaHistPos->Write();
  thetaHistNeg->Write();

  thetaHistPos_zfix->GetZaxis()->SetRangeUser(0, maxZth);
  thetaHistNeg_zfix->GetZaxis()->SetRangeUser(0, maxZth);

  thetaHistPos_zfix->Write();
  thetaHistNeg_zfix->Write();
  
  
  out_file->Close();
  
}


