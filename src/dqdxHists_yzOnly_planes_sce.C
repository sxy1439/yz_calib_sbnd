#include "../include/utilities_planes.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>
#include "../include/SCECorr.h"


// need this additional setup in gpvm to run this macro with SCE
// export SBND_DATA_PATH=/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/
// export SBNDDATA_VERSION=v01_28_00


SCECorr *sce_corr_mc = new SCECorr(false);

// to run :
//root [0] .L dqdxHists.C
//root [1] dqdxHists("doMC2023B_sub1", "output_files_mpv/yz_mc2023B_sub1.root", "output_files_mpv/x_mc2023B_sub1.root")


//void dqdxHists_yzOnly_planes(const char* cintyp, const char* inYZOut, const char* inSCEYZ) {
void dqdxHists_yzOnly_planes_sce(const char* cintyp, const char* inSCEYZ) {

  //LoadSCEMaps();
  sce_corr_mc->ReadHistograms();

  /*
  TH2F *CzyHist[nplanes][2];
  TFile* file_YZOut = new TFile(inYZOut, "READ");
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      CzyHist[l][k] = (TH2F*)file_YZOut->Get(Form("CzyHist_%i_%i",l,k));
    }
  }
  */

  TH2F *CzyHist_sce[nplanes][2];
  TFile* file_SCEYZ = new TFile(inSCEYZ, "READ");
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      CzyHist_sce[l][k] = (TH2F*)file_SCEYZ->Get(Form("CzyHist_%i_%i",l,k));
    }
  }


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


  cout<<"nbinx = "<<nbinx<<" nbiny = "<<nbiny<<" nbinz = "<<nbinz<<endl;
  cout<<"elife[0] = "<<e_callife[0]<<" elife[1] = "<<e_callife[1]<<endl;

  /// FINAL dQdx VALUES (dQ/dx with and without equalization correction)

  int ibinx, ibiny, ibinz;
  int ntrk=0;
  myReader.Restart();
  while (myReader.Next()) {
    if (!*is_mu) continue; //Pandora clear muon
    if ( !Is_Edge(*startx, *starty, *startz) || !Is_Edge(*endx, *endy, *endz)) continue;//FV
    if(!Is_Cathode_Crossing(*startx, *endx)) continue;
    
    thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi(); 
    if(*dirx<0) thetaxz = -thetaxz;
    thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
    if(*diry<0) thetayz = -thetayz;
    
    if(abs(thetaxz)<115&&abs(thetaxz)>65)continue;//Angle
    if(abs(thetayz)<110&&abs(thetayz)>70)continue;//Angle

    if (ntrk % 1000 == 0 && ntrk > 0) {
      std::cout << "entry => " << ntrk / 1000 << "k tracks" << std::endl;
    }
    
    ntrk++;
    

    float trk_t0 = *t0;

    // plane 0
    
    for (unsigned i = 0; i < dqdx0.GetSize(); i++) {
      if(isnan(tpx0[i])||isnan(tpy0[i])||isnan(tpz0[i]))continue;
      if(isnan(dqdx0[i]) || isinf(dqdx0[i])) continue;

      // sce corrected

      tpz0[i] = tpz0[i] - 4.2;
      
      XYZVector sp_sce_uncorr(tpx0[i], tpy0[i], tpz0[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);

      /*
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx0[i]<0){
	  if(InVeto_region_eastTPC_I0(tpy0[i], tpz0[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_I0(tpy0[i], tpz0[i])) continue;
	}
      }
      */

      ibinx = floor((sp_sce_corr.X()-lowx)/(highx-lowx)*nbinx);     
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((sp_sce_corr.Y()-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((sp_sce_corr.Z()-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;

      // electron lifetime
      double elife0 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      // YZ non-uniformity
      double CF_zy0 = CzyHist_sce[0][0]->GetBinContent(ibinz+1, ibiny+1);
      double CF_zy1 = CzyHist_sce[0][1]->GetBinContent(ibinz+1, ibiny+1);

      // sce
      double pitch_sce_uncorr = sce_corr_mc->meas_pitch(tpx0[i], tpy0[i], tpz0[i], tpdirx0[i], tpdiry0[i], tpdirz0[i], 0, false);
      double pitch_sce_corr = sce_corr_mc->meas_pitch(tpx0[i], tpy0[i], tpz0[i], tpdirx0[i], tpdiry0[i], tpdirz0[i], 0, true);
      double dqdx_sce_corr = dqdx0[i] * pitch_sce_uncorr/pitch_sce_corr;


      if(tpx0[i]<0){
	dqdxHist[0][0]->Fill(dqdx0[i] * elife0);
	dqdxHist[0][1]->Fill(dqdx_sce_corr * elife0);
	dqdxHist[0][2]->Fill(dqdx_sce_corr * elife0 * CF_zy0);

	dqdxHist[0][3]->Fill(dqdx0[i] * elife0);    // for ratio
	dqdxHist[0][4]->Fill(dqdx_sce_corr * elife0);
	dqdxHist[0][5]->Fill(dqdx_sce_corr * elife0 * CF_zy0);
      }
      else{

	dqdxHist[0][0]->Fill(dqdx0[i] * elife1);
	dqdxHist[0][1]->Fill(dqdx_sce_corr * elife1);
	dqdxHist[0][2]->Fill(dqdx_sce_corr * elife1 * CF_zy1);

	dqdxHist[0][3]->Fill(dqdx0[i] * elife1);    // for ratio
	dqdxHist[0][4]->Fill(dqdx_sce_corr * elife1);
	dqdxHist[0][5]->Fill(dqdx_sce_corr * elife1 * CF_zy1);
      }      
    }

    // plane 1
    
    for (unsigned i = 0; i < dqdx1.GetSize(); i++) {
      if(isnan(tpx1[i])||isnan(tpy1[i])||isnan(tpz1[i]))continue;
      if(isnan(dqdx1[i]) || isinf(dqdx1[i])) continue;

      // sce corrected

      tpz1[i] = tpz1[i] - 4.2;
      
      XYZVector sp_sce_uncorr(tpx1[i], tpy1[i], tpz1[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);

      /*
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx1[i]<0){
	  if(InVeto_region_eastTPC_C(tpy1[i], tpz1[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy1[i], tpz1[i])) continue;
	}
      }
      */

      ibinx = floor((sp_sce_corr.X()-lowx)/(highx-lowx)*nbinx);     
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((sp_sce_corr.Y()-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((sp_sce_corr.Z()-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;

      // electron lifetime
      double elife0 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      // YZ non-uniformity
      double CF_zy0 = CzyHist_sce[1][0]->GetBinContent(ibinz+1, ibiny+1);
      double CF_zy1 = CzyHist_sce[1][1]->GetBinContent(ibinz+1, ibiny+1);

      // sce
      double pitch_sce_uncorr = sce_corr_mc->meas_pitch(tpx1[i], tpy1[i], tpz1[i], tpdirx1[i], tpdiry1[i], tpdirz1[i], 1, false);
      double pitch_sce_corr = sce_corr_mc->meas_pitch(tpx1[i], tpy1[i], tpz1[i], tpdirx1[i], tpdiry1[i], tpdirz1[i], 1, true);
      double dqdx_sce_corr = dqdx1[i] * pitch_sce_uncorr/pitch_sce_corr;


      if(tpx1[i]<0){
	dqdxHist[1][0]->Fill(dqdx1[i] * elife0);
	dqdxHist[1][1]->Fill(dqdx_sce_corr * elife0);
	dqdxHist[1][2]->Fill(dqdx_sce_corr * elife0 * CF_zy0);

	dqdxHist[1][3]->Fill(dqdx1[i] * elife0);    // for ratio
	dqdxHist[1][4]->Fill(dqdx_sce_corr * elife0);
	dqdxHist[1][5]->Fill(dqdx_sce_corr * elife0 * CF_zy0);
      }
      else{

	dqdxHist[1][0]->Fill(dqdx1[i] * elife1);
	dqdxHist[1][1]->Fill(dqdx_sce_corr * elife1);
	dqdxHist[1][2]->Fill(dqdx_sce_corr * elife1 * CF_zy1);

	dqdxHist[1][3]->Fill(dqdx1[i] * elife1);    // for ratio
	dqdxHist[1][4]->Fill(dqdx_sce_corr * elife1);
	dqdxHist[1][5]->Fill(dqdx_sce_corr * elife1 * CF_zy1);
      }      
    }



    // plane 2
    
    for (unsigned i = 0; i < dqdx2.GetSize(); i++) {
      if(isnan(tpx2[i])||isnan(tpy2[i])||isnan(tpz2[i]))continue;
      if(isnan(dqdx2[i]) || isinf(dqdx2[i])) continue;

      // sce corrected

      tpz2[i] = tpz2[i] - 4.2;
      
      XYZVector sp_sce_uncorr(tpx2[i], tpy2[i], tpz2[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);

      /*
      // masked YZ and X regions
      if(string(cintyp).find("Data") != string::npos){
	if(tpx2[i]<0){
	  if(InVeto_region_eastTPC_C(tpy2[i], tpz2[i])) continue;
	}
	else{
	  if(InVeto_region_westTPC_C(tpy2[i], tpz2[i])) continue;
	}
      }
      */

      ibinx = floor((sp_sce_corr.X()-lowx)/(highx-lowx)*nbinx);     
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((sp_sce_corr.Y()-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((sp_sce_corr.Z()-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;

      // electron lifetime
      double elife0 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      // YZ non-uniformity
      double CF_zy0 = CzyHist_sce[2][0]->GetBinContent(ibinz+1, ibiny+1);
      double CF_zy1 = CzyHist_sce[2][1]->GetBinContent(ibinz+1, ibiny+1);

      // sce
      double pitch_sce_uncorr = sce_corr_mc->meas_pitch(tpx2[i], tpy2[i], tpz2[i], tpdirx2[i], tpdiry2[i], tpdirz2[i], 2, false);
      double pitch_sce_corr = sce_corr_mc->meas_pitch(tpx2[i], tpy2[i], tpz2[i], tpdirx2[i], tpdiry2[i], tpdirz2[i], 2, true);
      double dqdx_sce_corr = dqdx2[i] * pitch_sce_uncorr/pitch_sce_corr;


      if(tpx2[i]<0){
	dqdxHist[2][0]->Fill(dqdx2[i] * 50 * elife0);    // factor of 50 is just to convert the scale to e/cm (dqdx [ADC/cm] = 50 * dqdx [e/cm])
	dqdxHist[2][1]->Fill(dqdx_sce_corr * 50 * elife0);
	dqdxHist[2][2]->Fill(dqdx_sce_corr * 50 * elife0 * CF_zy0);

	dqdxHist[2][3]->Fill(dqdx2[i] * 50 * elife0);    // for ratio
	dqdxHist[2][4]->Fill(dqdx_sce_corr * 50 * elife0);
	dqdxHist[2][5]->Fill(dqdx_sce_corr * 50 * elife0 * CF_zy0);
      }
      else{
	dqdxHist[2][0]->Fill(dqdx2[i] * 50 * elife1);
	dqdxHist[2][1]->Fill(dqdx_sce_corr * 50 * elife1);
	dqdxHist[2][2]->Fill(dqdx_sce_corr * 50 * elife1 * CF_zy1);

	dqdxHist[2][3]->Fill(dqdx2[i] * 50 * elife1);    // for ratio
	dqdxHist[2][4]->Fill(dqdx_sce_corr * 50 * elife1);
	dqdxHist[2][5]->Fill(dqdx_sce_corr * 50 * elife1 * CF_zy1);
      }      
    }
    
  }

  for(int l=0;l<nplanes;l++){
    dqdxHist[l][3]->Divide(dqdxHist[l][0]);
    dqdxHist[l][4]->Divide(dqdxHist[l][0]);
    dqdxHist[l][5]->Divide(dqdxHist[l][0]);
  }
  
  for(int l=0;l<nplanes;l++)for(int k=0;k<6;k++)dqdxHist[l][k]->Write();

  //TGaxis::SetMaxDigits(3);

  for(int l=0;l<nplanes;l++){
    TCanvas *c1 = new TCanvas();
    
    dqdxHist[l][0]->SetLineColor(kBlack);
    dqdxHist[l][0]->SetLineWidth(2);
    dqdxHist[l][0]->SetStats(0);
    
    dqdxHist[l][1]->SetLineColor(kGreen+2);
    dqdxHist[l][1]->SetLineWidth(2);
    dqdxHist[l][1]->SetStats(0);
    
    dqdxHist[l][2]->SetLineColor(kAzure+2);
    dqdxHist[l][2]->SetLineWidth(2);
    dqdxHist[l][2]->SetStats(0);
    
    dqdxHist[l][0]->Draw("hist");
    dqdxHist[l][1]->Draw("hist same");
    dqdxHist[l][2]->Draw("hist same");
    
    
    TLegend *l1 = new TLegend(0.6, 0.7, 0.85, 0.85);
    l1->SetTextSize(0.030);
    l1->AddEntry(dqdxHist[l][0], "No Correction", "l");
    l1->AddEntry(dqdxHist[l][1], "+ SCE correction", "l");
    l1->AddEntry(dqdxHist[l][2], "+ YZ correction", "l");
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    
    l1->Draw();
    c1->SetLeftMargin(0.15);
    c1->SaveAs(Form("plots_dqdxHists_planes/sce/finaldqdx_%i_%s.png", l, cintyp));
    c1->SaveAs(Form("plots_dqdxHists_planes/sce/finaldqdx_%i_%s.pdf", l, cintyp));
  }


  for(int l=0;l<nplanes;l++){
    TCanvas *c2 = new TCanvas();
    
    dqdxHist[l][3]->SetLineColor(kBlack);
    dqdxHist[l][3]->SetLineWidth(2);
    dqdxHist[l][3]->SetStats(0);
    
    dqdxHist[l][4]->SetLineColor(kGreen+2);
    dqdxHist[l][4]->SetLineWidth(2);
    dqdxHist[l][4]->SetStats(0);
    
    dqdxHist[l][5]->SetLineColor(kAzure+2);
    dqdxHist[l][5]->SetLineWidth(2);
    dqdxHist[l][5]->SetStats(0);

    dqdxHist[l][3]->GetYaxis()->SetTitle("Ratio to nominal");

    dqdxHist[l][3]->Draw("hist");
    dqdxHist[l][4]->Draw("hist same");
    dqdxHist[l][5]->Draw("hist same");
    
    TLegend *l2 = new TLegend(0.6, 0.7, 0.85, 0.85);
    l2->SetTextSize(0.030);
    l2->AddEntry(dqdxHist[l][3], "No Correction", "l");
    l2->AddEntry(dqdxHist[l][4], "+ lifetime correction", "l");
    l2->AddEntry(dqdxHist[l][5], "+ SCE + YZ correction", "l");
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    
    l2->Draw();
    c2->SetLeftMargin(0.15);
    c2->SaveAs(Form("plots_dqdxHists_planes/sce/dqdx_ratio_%i_%s.png", l, cintyp));
    c2->SaveAs(Form("plots_dqdxHists_planes/sce/dqdx_ratio_%i_%s.pdf", l, cintyp));
  }

  out_file->Close();
  
}
