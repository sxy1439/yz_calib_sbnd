
#include "../include/utilities_planes.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>

// to run : root -l -b -q 'ang_ana.C("doMC_sub1")' > mpv_out_mc.txt

// to run :
//root [0] .L ang_ana.C
//root [1] ang_ana("doMC2023B_sub1")


void ang_ana1D(const char* cintyp) {

  gSystem->Load("libFileHandling.so");
  if (gSystem->Load("libFileHandling.so") < 0) {
    std::cerr << "Error loading library" << std::endl;
    return;
  }

  auto e_callife = getELifetime(cintyp);
 
  TString output_name = getOutputNameChannelByAngle(cintyp);
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


  TTreeReaderArray<float> dqdx0(myReader, "trk.hits0.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> tpx0(myReader, "trk.hits0.h.sp.x"); // x of track trajectory position  (older was -> trk.hits2.tp.x)
  TTreeReaderArray<float> tpy0(myReader, "trk.hits0.h.sp.y"); // y of track trajectory position  
  TTreeReaderArray<float> tpz0(myReader, "trk.hits0.h.sp.z"); // z of track trajectory position
  TTreeReaderArray<float> tpdirx0(myReader, "trk.hits0.dir.x");
  TTreeReaderArray<float> tpdiry0(myReader, "trk.hits0.dir.y");
  TTreeReaderArray<float> tpdirz0(myReader, "trk.hits0.dir.z");
  TTreeReaderArray<float> tppitch0(myReader, "trk.hits0.pitch");
  TTreeReaderArray<float> time0(myReader, "trk.hits0.h.time"); // in ticks (500 ns), up to 3200 (1.6 ms > 1.25 ms)
  TTreeReaderArray<UShort_t> channel0(myReader, "trk.hits0.h.channel");
  TTreeReaderArray<float> charge0(myReader, "trk.hits0.h.integral");

  TTreeReaderArray<float> dqdx1(myReader, "trk.hits1.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> tpx1(myReader, "trk.hits1.h.sp.x"); // x of track trajectory position  (older was -> trk.hits2.tp.x)
  TTreeReaderArray<float> tpy1(myReader, "trk.hits1.h.sp.y"); // y of track trajectory position  
  TTreeReaderArray<float> tpz1(myReader, "trk.hits1.h.sp.z"); // z of track trajectory position
  TTreeReaderArray<float> tpdirx1(myReader, "trk.hits1.dir.x");
  TTreeReaderArray<float> tpdiry1(myReader, "trk.hits1.dir.y");
  TTreeReaderArray<float> tpdirz1(myReader, "trk.hits1.dir.z");
  TTreeReaderArray<float> tppitch1(myReader, "trk.hits1.pitch");
  TTreeReaderArray<float> time1(myReader, "trk.hits1.h.time"); // in ticks (500 ns), up to 3200 (1.6 ms > 1.25 ms)
  TTreeReaderArray<UShort_t> channel1(myReader, "trk.hits1.h.channel");
  TTreeReaderArray<float> charge1(myReader, "trk.hits1.h.integral");

  TTreeReaderArray<float> dqdx2(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> tpx2(myReader, "trk.hits2.h.sp.x"); // x of track trajectory position  (older was -> trk.hits2.tp.x)
  TTreeReaderArray<float> tpy2(myReader, "trk.hits2.h.sp.y"); // y of track trajectory position  
  TTreeReaderArray<float> tpz2(myReader, "trk.hits2.h.sp.z"); // z of track trajectory position
  TTreeReaderArray<float> tpdirx2(myReader, "trk.hits2.dir.x");
  TTreeReaderArray<float> tpdiry2(myReader, "trk.hits2.dir.y");
  TTreeReaderArray<float> tpdirz2(myReader, "trk.hits2.dir.z");
  TTreeReaderArray<float> tppitch2(myReader, "trk.hits2.pitch");
  TTreeReaderArray<float> time2(myReader, "trk.hits2.h.time"); // in ticks (500 ns), up to 3200 (1.6 ms > 1.25 ms)
  TTreeReaderArray<UShort_t> channel2(myReader, "trk.hits2.h.channel");
  TTreeReaderArray<float> charge2(myReader, "trk.hits2.h.integral");
  


  cout<<"nbiny = "<<nbiny<<" nbinz = "<<nbinz<<endl;

  
  /// ANGULAR REQUIREMENT STUDY

  // Loop over all entries of the TTree
  int ntrkpre=0;
  int ntrkpospre=0, ntrknegpre=0;

  std::vector<double> dqdx_values_Pos[3][nbinang];  // 3 planes, nbinx θ bins
  std::vector<double> dqdx_values_Neg[3][nbinang];


  myReader.Restart();
  while (myReader.Next()) {

    //if(isnan(*startx) || isnan(*starty) || isnan(*startz)) continue;
    //if(isnan(*endx) || isnan(*endy) || isnan(*endz)) continue;
    //if(isnan(*dirx) || isnan(*diry) || isnan(*dirz)) continue;
    
    if (!*is_mu) continue; 
    if ( !Is_Edge(*startx, *starty, *startz) || !Is_Edge(*endx, *endy, *endz)) continue;
    if(!Is_Cathode_Crossing(*startx, *endx)) continue;
    
    thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi();   
    if(*dirx<0) thetaxz = -thetaxz;
    
    thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
    if(*diry<0) thetayz = -thetayz;
    
    if(abs(thetaxz)<115&&abs(thetaxz)>65)continue;//Angle
    if(abs(thetayz)<110&&abs(thetayz)>70)continue;//Angle
    
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

    // plane 0
    
    for(unsigned i = 0; i < dqdx0.GetSize(); i++){
      if (isnan(tpx0[i]) || isnan(tpy0[i]) || isnan(tpz0[i])) continue;
      if(isnan(dqdx0[i]) || isinf(dqdx0[i])) continue;
      if(isnan(tpdirz0[i]) || isnan(tpdirx0[i])) continue;
      
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx0[i]<0){
	  if(InVeto_region_eastTPC_C(tpy0[i], tpz0[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy0[i], tpz0[i])) continue;
	}
      }

      double elife0 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      // median

      if (tpx0[i] < 0) {
	int bin = thetaHistNeg1D_med[0]->FindBin(thetaxz);
	dqdx_values_Neg[0][bin-1].push_back(dqdx0[i] * elife0);
      } else {
	int bin = thetaHistPos1D_med[0]->FindBin(thetaxz);
	dqdx_values_Pos[0][bin-1].push_back(dqdx0[i] * elife1);
      }
      
    }

    // plane 1
    
    for(unsigned i = 0; i < dqdx1.GetSize(); i++){
      if (isnan(tpx1[i]) || isnan(tpy1[i]) || isnan(tpz1[i])) continue;
      if(isnan(dqdx1[i]) || isinf(dqdx1[i])) continue;
      if(isnan(tpdirz1[i]) || isnan(tpdirx1[i])) continue;
      
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx1[i]<0){
	  if(InVeto_region_eastTPC_C(tpy1[i], tpz1[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy1[i], tpz1[i])) continue;
	}
      }

      double elife0 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      // median

      if (tpx1[i] < 0) {
	int bin = thetaHistNeg1D_med[1]->FindBin(thetaxz);
	dqdx_values_Neg[1][bin-1].push_back(dqdx1[i] * elife0);
      } else {
	int bin = thetaHistPos1D_med[1]->FindBin(thetaxz);
	dqdx_values_Pos[1][bin-1].push_back(dqdx1[i] * elife1);
      }
      
    }

    // plane 2
    
    for(unsigned i = 0; i < dqdx2.GetSize(); i++){
      if (isnan(tpx2[i]) || isnan(tpy2[i]) || isnan(tpz2[i])) continue;
      if(isnan(dqdx2[i]) || isinf(dqdx2[i])) continue;
      if(isnan(tpdirz2[i]) || isnan(tpdirx2[i])) continue;
      
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx2[i]<0){
	  if(InVeto_region_eastTPC_C(tpy2[i], tpz2[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy2[i], tpz2[i])) continue;
	}
      }

      double elife0 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      // median

      if (tpx2[i] < 0) {
	int bin = thetaHistNeg1D_med[2]->FindBin(thetaxz);    // in which thetaxz bin this hit belong
	dqdx_values_Neg[2][bin-1].push_back(dqdx2[i] * elife0);
      } else {
	int bin = thetaHistPos1D_med[2]->FindBin(thetaxz);
	dqdx_values_Pos[2][bin-1].push_back(dqdx2[i] * elife1);
      }
      
    }
      
    
  }

  
  for (int plane = 0; plane < 3; plane++) {
    for (int bin = 1; bin <= thetaHistPos1D_med[plane]->GetNbinsX(); bin++) {
      //if (bin-1 < 0 || bin-1 >= nbinang) continue;
      
      float median_pos = getMedian(dqdx_values_Pos[plane][bin-1]);
      thetaHistPos1D_med[plane]->SetBinContent(bin, median_pos);     
      float median_neg = getMedian(dqdx_values_Neg[plane][bin-1]);
      thetaHistNeg1D_med[plane]->SetBinContent(bin, median_neg);

      if (std::isnan(median_pos) || std::isinf(median_pos)) {
	std::cout << "WARNING: median_pos is NaN or Inf for plane " << plane << ", bin " << bin << std::endl;
      }
      if (std::isnan(median_neg) || std::isinf(median_neg)) {
	std::cout << "WARNING: median_neg is NaN or Inf for plane " << plane << ", bin " << bin << std::endl;
      }

      if (dqdx_values_Pos[plane][bin-1].empty()) {
	std::cout << "Bin " << bin << " for plane " << plane << " is empty (Pos)." << std::endl;
      }
      if (dqdx_values_Neg[plane][bin-1].empty()) {
	std::cout << "Bin " << bin << " for plane " << plane << " is empty (Neg)." << std::endl;
      }


      if (thetaHistPos1D_med[plane]->GetEntries() == 0) {
	thetaHistPos1D_med[plane]->SetBinContent(1, 0.001);  // small dummy value
      }
      if (thetaHistNeg1D_med[plane]->GetEntries() == 0) {
	thetaHistNeg1D_med[plane]->SetBinContent(1, 0.001);
      }

    }
    std::cout << "thetaHistPos1D_med[" << plane << "] nbins: "<< thetaHistPos1D_med[plane]->GetNbinsX() << std::endl;

  }
  

  for(int l=0; l<3; l++){

    thetaHistPos1D_med[l]->GetXaxis()->CenterTitle();
    thetaHistPos1D_med[l]->GetYaxis()->CenterTitle();
    thetaHistNeg1D_med[l]->GetXaxis()->CenterTitle();
    thetaHistNeg1D_med[l]->GetYaxis()->CenterTitle();
     
    thetaHistPos1D_med[l]->SetXTitle("#theta_{xz} [deg]");
    thetaHistPos1D_med[l]->SetYTitle("Median dQ/dx [ADC/cm]");
    thetaHistNeg1D_med[l]->SetXTitle("#theta_{xz} [deg]");
    thetaHistNeg1D_med[l]->SetYTitle("Median dQ/dx [ADC/cm]");

    thetaHistPos1D_med[l]->Write();
    thetaHistNeg1D_med[l]->Write();
    
  }


  cout<<"Total tracks (before angular cut): "<<ntrkpre<<endl;
  cout<<"# tracks in x>=0 region (before angular cut): "<<ntrkpospre<<endl;
  cout<<"# tracks in x<0 region (before angular cut): "<<ntrknegpre<<endl;

  
  
  out_file->Close();
  
}


