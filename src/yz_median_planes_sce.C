
#include "../include/utilities_planes.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>
#include "../include/SCECorr.h"


// need this additional setup in gpvm to run this macro with SCE
// export SBND_DATA_PATH=/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/
// export SBNDDATA_VERSION=v01_28_00


SCECorr *sce_corr_mc = new SCECorr(false);

void yz_median_planes_sce(const char* cintyp) {

  //LoadSCEMaps();
  sce_corr_mc->ReadHistograms();

  gSystem->Load("libFileHandling.so");
  if (gSystem->Load("libFileHandling.so") < 0) {
    std::cerr << "Error loading library" << std::endl;
    return;
  }

  auto e_callife = getELifetime(cintyp);
  
  TString output_name = getOutputNameYZ(cintyp);
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
  //TTreeReaderValue<float> t0(myReader, "trk.t0PFP");
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
  cout<<"elife[0] = "<<e_callife[0]<<" elife[1] = "<<e_callife[1]<<endl;


  /// YZ equalization
  
  int ntrk=0;
  int ntrkpos=0, ntrkneg=0;
  int ibinx, ibiny, ibinz, ibinc;
 
  myReader.Restart();
  while (myReader.Next()) {  
    if (!*is_mu) continue; 
    if ( !Is_Edge(*startx, *starty, *startz) || !Is_Edge(*endx, *endy, *endz)) continue;//FV
    if(!Is_Cathode_Crossing(*startx, *endx)) continue;
    
    thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi();   
    if(*dirx<0) thetaxz = -thetaxz;
    
    thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
    if(*diry<0) thetayz = -thetayz;

    //if(abs(thetaxz)>25 && abs(thetaxz)<35) continue;
    //if(abs(thetaxz)>145 && abs(thetaxz)<155) continue;
    
    if(abs(thetaxz)<115&&abs(thetaxz)>65)continue;//Angle      // TO DO : revist
    if(abs(thetayz)<110&&abs(thetayz)>70)continue;//Angle

    if (ntrk % 1000 == 0 && ntrk > 0) {
      std::cout << "entry => " << ntrk / 1000 << "k tracks" << std::endl;
    }
    
    ntrk++;

    if (*startx >= 0){
      ntrkpos++;
    } else {
      ntrkneg++;
    }

    float trk_t0 = *t0;

    
    // plane 0
    
    for (unsigned i = 0; i < dqdx0.GetSize(); i++) {
      if(isnan(tpx0[i])||isnan(tpy0[i])||isnan(tpz0[i]))continue;
      if(isnan(dqdx0[i]) || isinf(dqdx0[i])) continue;
      
      // sce corrected

      // *** manual shift in z (only till we get this fixed in sbnd geometry)
      // ** remove this line once it's fixed    
      tpz0[i] = tpz0[i] - 4.2;
      
      XYZVector sp_sce_uncorr(tpx0[i], tpy0[i], tpz0[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);

      /*
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx0[i]<0){
	  if(InVeto_region_eastTPC_C(sp_sce_corr.Y(), sp_sce_corr.Z())) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(sp_sce_corr.Y(), sp_sce_corr.Z())) continue;
	}
      }
      */

      //cout<<"spatial postions (Uncorrected) : "<<Form("x(%f), y(%f), z(%f)", tpx0[i], tpy0[i], tpz0[i])<<endl;
      //cout<<"spatial postions (Corrected) : "<<Form("x(%f), y(%f), z(%f)", sp_sce_corr.X(), sp_sce_corr.Y(), sp_sce_corr.Z())<<endl;

      ibinx = floor((sp_sce_corr.X()-lowx)/(highx-lowx)*nbinx);       
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((sp_sce_corr.Y()-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((sp_sce_corr.Z()-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;
      
      double pitch_sce_uncorr = sce_corr_mc->meas_pitch(tpx0[i], tpy0[i], tpz0[i], tpdirx0[i], tpdiry0[i], tpdirz0[i], 0, false);
      double pitch_sce_corr = sce_corr_mc->meas_pitch(tpx0[i], tpy0[i], tpz0[i], tpdirx0[i], tpdiry0[i], tpdirz0[i], 0, true);
      double dqdx_sce_corr = dqdx0[i] * pitch_sce_uncorr/pitch_sce_corr;

      double elife0 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      //cout<<"Plane 0 ==>>>"<<endl;
      //cout<<"Pitch (Uncorrected) : "<<pitch_sce_uncorr<<" Pitch (SCE Corrected) : "<<pitch_sce_corr<<endl;
      //cout<<"dQ/dx (Uncorrected) : "<<dqdx0[i]<<" dQ/dx (SCE Corrected) : "<<dqdx_sce_corr<<endl;
      
      if(tpx0[i]<0){
	dqdxHist[0][0]->Fill(dqdx0[i] * elife0);
	dqdxHist[0][1]->Fill(dqdx_sce_corr * elife0);
	zynhits[0][0]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
	zyHistdqdx[0][0]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife0); 
	nyzHist[0][0][ibiny][ibinz]->Fill(dqdx_sce_corr * elife0);
      }
      else{
	dqdxHist[0][0]->Fill(dqdx0[i] * elife1);
	dqdxHist[0][1]->Fill(dqdx_sce_corr * elife1);
	zynhits[0][1]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
	zyHistdqdx[0][1]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife1); 
	nyzHist[0][1][ibiny][ibinz]->Fill(dqdx_sce_corr * elife1);
      }
      
      
    }

    // plane 1
    
    for (unsigned i = 0; i < dqdx1.GetSize(); i++) {
      if(isnan(tpx1[i])||isnan(tpy1[i])||isnan(tpz1[i]))continue;
      if(isnan(dqdx1[i]) || isinf(dqdx1[i])) continue;

      // sce corrected

      // *** manual shift in z (only till we get this fixed in sbnd geometry)
      // ** remove this line once it's fixed    
      tpz1[i] = tpz1[i] - 4.2;

      XYZVector sp_sce_uncorr(tpx1[i], tpy1[i], tpz1[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);

      /*
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx1[i]<0){
	  if(InVeto_region_eastTPC_C(sp_sce_corr.Y(), sp_sce_corr.Z())) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(sp_sce_corr.Y(), sp_sce_corr.Z())) continue;
	}
      }
      */

      //cout<<"spatial postions (Uncorrected) : "<<Form("x(%f), y(%f), z(%f)", tpx1[i], tpy1[i], tpz1[i])<<endl;
      //cout<<"spatial postions (Corrected) : "<<Form("x(%f), y(%f), z(%f)", sp_sce_corr.X(), sp_sce_corr.Y(), sp_sce_corr.Z())<<endl;
      

      ibinx = floor((sp_sce_corr.X()-lowx)/(highx-lowx)*nbinx);       
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((sp_sce_corr.Y()-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((sp_sce_corr.Z()-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;
      
      double pitch_sce_uncorr = sce_corr_mc->meas_pitch(tpx1[i], tpy1[i], tpz1[i], tpdirx1[i], tpdiry1[i], tpdirz1[i], 1, false);
      double pitch_sce_corr = sce_corr_mc->meas_pitch(tpx1[i], tpy1[i], tpz1[i], tpdirx1[i], tpdiry1[i], tpdirz1[i], 1, true);
      double dqdx_sce_corr = dqdx1[i] * pitch_sce_uncorr/pitch_sce_corr;

      double elife0 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      //cout<<"Plane 1 ==>>>"<<endl;
      //cout<<"Pitch (Uncorrected) : "<<pitch_sce_uncorr<<" Pitch (SCE Corrected) : "<<pitch_sce_corr<<endl;
      //cout<<"dQ/dx (Uncorrected) : "<<dqdx1[i]<<" dQ/dx (SCE Corrected) : "<<dqdx_sce_corr<<endl;
      
      if(tpx1[i]<0){
	dqdxHist[1][0]->Fill(dqdx1[i] * elife0);
	dqdxHist[1][1]->Fill(dqdx_sce_corr * elife0);
	zynhits[1][0]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
	zyHistdqdx[1][0]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife0); 
	nyzHist[1][0][ibiny][ibinz]->Fill(dqdx_sce_corr * elife0);
      }
      else{
	dqdxHist[1][0]->Fill(dqdx1[i] * elife1);
	dqdxHist[1][1]->Fill(dqdx_sce_corr * elife1);
	zynhits[1][1]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
	zyHistdqdx[1][1]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife1); 
	nyzHist[1][1][ibiny][ibinz]->Fill(dqdx_sce_corr * elife1);
      }
      
      
    }

    // plane 2
    
    for (unsigned i = 0; i < dqdx2.GetSize(); i++) {
      if(isnan(tpx2[i])||isnan(tpy2[i])||isnan(tpz2[i]))continue;
      if(isnan(dqdx2[i]) || isinf(dqdx2[i])) continue;

      if(tpx2[i]<0){
	if(tpz2[i]>65 && tpz2[i]<70) continue;   // to avoid floating point leakage
      }

      // sce corrected

      // *** manual shift in z (only till we get this fixed in sbnd geometry)
      // ** remove this line once it's fixed    
      tpz2[i] = tpz2[i] - 4.2;

      XYZVector sp_sce_uncorr(tpx2[i], tpy2[i], tpz2[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);

      /*
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx2[i]<0){
	  if(InVeto_region_eastTPC_C(sp_sce_corr.Y(), sp_sce_corr.Z())) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(sp_sce_corr.Y(), sp_sce_corr.Z())) continue;
	}
      }
      */

      //cout<<"spatial postions (Uncorrected) : "<<Form("x(%f), y(%f), z(%f)", tpx2[i], tpy2[i], tpz2[i])<<endl;
      //cout<<"spatial postions (Corrected) : "<<Form("x(%f), y(%f), z(%f)", sp_sce_corr.X(), sp_sce_corr.Y(), sp_sce_corr.Z())<<endl;

      ibinx = floor((sp_sce_corr.X()-lowx)/(highx-lowx)*nbinx);       
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((sp_sce_corr.Y()-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((sp_sce_corr.Z()-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;
      
      double pitch_sce_uncorr = sce_corr_mc->meas_pitch(tpx2[i], tpy2[i], tpz2[i], tpdirx2[i], tpdiry2[i], tpdirz2[i], 2, false);
      double pitch_sce_corr = sce_corr_mc->meas_pitch(tpx2[i], tpy2[i], tpz2[i], tpdirx2[i], tpdiry2[i], tpdirz2[i], 2, true);
      double dqdx_sce_corr = dqdx2[i] * pitch_sce_uncorr/pitch_sce_corr;

      double elife0 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      //cout<<"Plane 2 ==>>>"<<endl;
      //cout<<"Pitch (Uncorrected) : "<<pitch_sce_uncorr<<" Pitch (SCE Corrected) : "<<pitch_sce_corr<<endl;
      //cout<<"dQ/dx (Uncorrected) : "<<dqdx2[i]<<" dQ/dx (SCE Corrected) : "<<dqdx_sce_corr<<endl;
      
      if(tpx2[i]<0){
	dqdxHist[2][0]->Fill(dqdx2[i] * elife0);
	dqdxHist[2][1]->Fill(dqdx_sce_corr * elife0);
	zynhits[2][0]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
	zyHistdqdx[2][0]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife0); 
	nyzHist[2][0][ibiny][ibinz]->Fill(dqdx_sce_corr * elife0);
      }
      else{
	dqdxHist[2][0]->Fill(dqdx2[i] * elife1);
	dqdxHist[2][1]->Fill(dqdx_sce_corr * elife1);
	zynhits[2][1]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
	zyHistdqdx[2][1]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife1); 
	nyzHist[2][1][ibiny][ibinz]->Fill(dqdx_sce_corr * elife1);
      }
      
      
    }
    
    //cout<<"[DEBUG:] done with the dqdx loop"<<endl;
    
  }
  cout<<"Total tracks: "<<ntrk<<endl;
  cout<<"# tracks in x>=0 region: "<<ntrkpos<<endl;
  cout<<"# tracks in x<0 region: "<<ntrkneg<<endl;
  
  const int nq=1;
  double xq[nq]={0.5}, yq[nq];
  std::vector<double> medianyz_vec[nplanes][2]; // nplanes=3 for three planes U, V, Y
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      for(int i=0;i<nbiny;i++)for(int j=0;j<nbinz;j++){
	  if(nyzHist[l][k][i][j]->Integral()<2)continue;
	  
	  nyzHist[l][k][i][j]->GetQuantiles(nq, yq, xq);
	  
	  if(yq[0]<2)continue;
	  zyHist[l][k]->SetBinContent(j+1, i+1, yq[0]);    // takes the median of dqdx in that bin (to get local dqdx)

	  medianyz_vec[l][k].push_back(yq[0]);     // for median of medians (to get global dqdx)
	  
	}
    }
  }

  cout<<"[DEBUG:] Got median of medians"<<endl;

  /*
  for(int l=0;l<3;l++){
    for(int k=0;k<2;k++){
      for(int i=0;i<nbiny;i++)for(int j=0;j<nbinz;j++){
	  nyzHist[l][k][i][j]->Write();
	}
    }
  }
  */

 

  // for (dQ/dx)_{X}^{global}

  double global_dqdx_medianyz[nplanes][2];
  for(int l=0;l<nplanes;l++){
    for(int k=0; k<2; k++){
      global_dqdx_medianyz[l][k] = getMedian(medianyz_vec[l][k]);
      cout<< "global median : "<<global_dqdx_medianyz[l][k]<<endl;
    }
  }

  cout<<"[DEBUG:] Got global median"<<endl;

  //double maxZ = std::max(zyHist[0]->GetMaximum(), zyHist[1]->GetMaximum());
  
  TH2F *CzyHist[nplanes][2];
  for(int l=0;l<nplanes;l++){
    double maxZ = std::max(zyHist[l][0]->GetMaximum(), zyHist[l][1]->GetMaximum());
    for(int k=0;k<2;k++){
      CzyHist[l][k] = (TH2F*)zyHist[l][k]->Clone(Form("CzyHist_%i_%i",l,k));
      for(int i=0;i<nbiny;i++)for(int j=0;j<nbinz;j++){
	  if(zyHist[l][k]->GetBinContent(j+1, i+1)<1e-2)continue;
	  
	  CzyHist[l][k]->SetBinContent(j+1, i+1, global_dqdx_medianyz[l][k] / zyHist[l][k]->GetBinContent(j+1, i+1));	
	}
      
      CzyHist[l][k]->GetXaxis()->CenterTitle();
      CzyHist[l][k]->GetYaxis()->CenterTitle();
      CzyHist[l][k]->GetZaxis()->CenterTitle();
      
      CzyHist[l][k]->SetXTitle("z [cm]");
      CzyHist[l][k]->SetYTitle("y [cm]");
      CzyHist[l][k]->SetZTitle("YZ correction factor");

      /*
      //zyHist[l][k]->GetZaxis()->SetRangeUser(400, maxZ);
      zyHist[l][k]->GetZaxis()->SetRangeUser(400, 2000);
      if(string(cintyp).find("Data") != string::npos){
	CzyHist[l][k]->GetZaxis()->SetRangeUser(0.8, 1.2);
      }
      else CzyHist[l][k]->GetZaxis()->SetRangeUser(0.9, 1.1);
      */
      
      zyHist[l][k]->Write();     // median/mpv dqdx
      CzyHist[l][k]->Write();    // correction factor
      
      zynhits[l][k]->Write();          // stores number of hits
      zyHistdqdx[l][k]->Write();       // stores dqdx values
    }
    dqdxHist[l][0]->Write();
    dqdxHist[l][1]->Write();
  }
  
  cout<<"[DEBUG:] Done writing histograms"<<endl;

  
  
  out_file->Close();
  
}




// gStyle->SetOptStat(0);
