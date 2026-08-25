#include "../include/utilities_planes.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>
#include "../include/SCECorr.h"


/*
  This macro is used to make the median dQ/dx YZ map ROOT file in x-slices.
  It is based on yz_median_planes_sce.C, but with additional histograms
  for each x-slice. The number of x-slices and their boundaries can be
  adjusted in the GetXSectionFullRange() function.
*/

SCECorr *sce_corr_mc = new SCECorr(false);

// const int nxsections = 4;

// int GetXSectionFullRange(double x) {
//   if (x >= -200.0 && x < -100.0) return 0;
//   if (x >= -100.0 && x <    0.0) return 1;
//   if (x >=    0.0 && x <  100.0) return 2;
//   if (x >=  100.0 && x <= 200.0) return 3;

//   return -1;
// }


const int nxsections = 20;
const double xslice_low = -200;
const double xslice_high = 200;
const double xslice_width = 20.0;

int GetXSectionFullRange(double x) 
{
  if (x < xslice_low || x > xslice_high) return -1;

  int xs = int((x - xslice_low) / xslice_width);

  // Protect the upper boundary x = 200 cm or whatever you set the upper bound to be
  if (xs == nxsections) xs = nxsections - 1;

  if (xs < 0 || xs >= nxsections) return -1;
  return xs;
}







void yz_median_planes_sce_xslices(const char* cintyp) {
  
  //LoadSCEMaps();
  sce_corr_mc->ReadHistograms();

  gSystem->Load("libFileHandling.so");
  if (gSystem->Load("libFileHandling.so") < 0) {
    std::cerr << "Error loading library" << std::endl;
    return;
  }

  auto e_callife = getELifetime(cintyp);

  TString output_name = getOutputNameYZ(cintyp);
  output_name.ReplaceAll(".root", "_xslices.root");

  TFile *out_file = new TFile(output_name, "recreate");
  
  // Initialize histograms
  initialize_histograms();

  // Extra histograms for x-sliced YZ maps
  TH1F *nyzHist_x[nplanes][2][nxsections][nbiny][nbinz];
  TH2F *zyHist_x[nplanes][2][nxsections];
  TH2F *zyHistdqdx_x[nplanes][2][nxsections];
  TH2F *zynhits_x[nplanes][2][nxsections];

  // New: 1D dQ/dx histograms for each plane/TPC/x-slice.
  // Here k = 0 is East/negative-x TPC and k = 1 is West/positive-x TPC.
  TH1F *dqdxHist_x[nplanes][2][nxsections];

  for (int l = 0; l < nplanes; l++) {
    for (int k = 0; k < 2; k++) {
      for (int xs = 0; xs < nxsections; xs++) {

        zyHist_x[l][k][xs] =
          new TH2F(Form("zy_%i_%i_xsec%i", l, k, xs),
                   Form("zy median dQ/dx xsec%i", xs),
                   nbinz, lowz, highz,
                   nbiny, lowy, highy);

        zyHist_x[l][k][xs]->SetXTitle("z [cm]");
        zyHist_x[l][k][xs]->SetYTitle("y [cm]");
        zyHist_x[l][k][xs]->SetZTitle("Median(dQ/dx [ADC/cm])");

        zyHistdqdx_x[l][k][xs] =
          new TH2F(Form("zydqdx_%i_%i_xsec%i", l, k, xs),
                   Form("zy dQ/dx xsec%i", xs),
                   nbinz, lowz, highz,
                   nbiny, lowy, highy);

        zyHistdqdx_x[l][k][xs]->SetXTitle("z [cm]");
        zyHistdqdx_x[l][k][xs]->SetYTitle("y [cm]");
        zyHistdqdx_x[l][k][xs]->SetZTitle("dQ/dx [ADC/cm]");

        zynhits_x[l][k][xs] =
          new TH2F(Form("zynhits_%i_%i_xsec%i", l, k, xs),
                   Form("zy nhits xsec%i", xs),
                   nbinz, lowz, highz,
                   nbiny, lowy, highy);

        zynhits_x[l][k][xs]->SetXTitle("z [cm]");
        zynhits_x[l][k][xs]->SetYTitle("y [cm]");
        zynhits_x[l][k][xs]->SetZTitle("number of hits");

        dqdxHist_x[l][k][xs] =
          new TH1F(Form("dqdxHist_%i_%i_xsec%i", l, k, xs),
                   Form("dQ/dx plane %i TPC %i xsec%i", l, k, xs),
                   nbinq, lowq, highq);

        dqdxHist_x[l][k][xs]->GetXaxis()->SetTitle("dQ/dx [ADC/cm]");
        dqdxHist_x[l][k][xs]->GetYaxis()->SetTitle("Entries");

        for (int iy = 0; iy < nbiny; iy++) {
          for (int iz = 0; iz < nbinz; iz++) {
            nyzHist_x[l][k][xs][iy][iz] =
              new TH1F(Form("yz%i_%i_xsec%i_%i_%i", l, k, xs, iy, iz),
                       "",
                       nbinq, lowq, highq);
          }
        }
      }
    }
  }



  
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
  //TTreeReaderValue<float> t0(myReader, "trk.t0");
  TTreeReaderValue<float> t0(myReader, "trk.t0PFP");
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


  TTreeReaderValue<int> run(myReader, "trk.meta.run");
  TTreeReaderValue<int> subrun(myReader, "trk.meta.subrun");
  TTreeReaderValue<int> evt(myReader, "trk.meta.evt");


  cout<<"nbiny = "<<nbiny<<" nbinz = "<<nbinz<<endl;
  cout<<"elife[0] = "<<e_callife[0]<<" elife[1] = "<<e_callife[1]<<endl;


  /// YZ equalization
  
  int ntrk=0;
  int ntrkpos=0, ntrkneg=0;
  int ibinx, ibiny, ibinz, ibinc;
  int lastPrintedRun = -1;

  string title="selected-MC-metadata.csv";
  FILE *fptr=fopen(title.c_str(),"a");
 
  myReader.Restart();
  while (myReader.Next()) {  
    if (!*is_mu) continue; 
    if ( !Is_Edge(*startx, *starty, *startz) || !Is_Edge(*endx, *endy, *endz)) continue;//FV
    if(!Is_Cathode_Crossing(*startx, *endx)) continue;

    // pump trip period
    if(string(cintyp).find("Data") != string::npos && *run>=18503 && *run<=18523){
      if(*run != lastPrintedRun){
	std::cout<<"pump trip period run ["<<*run<<"] found. Skipping the run."<<std::endl;
	lastPrintedRun = *run;
      }
      continue;
    }
    
    thetaxz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*dirx,2)))*180/TMath::Pi();   
    if(*dirx<0) thetaxz = -thetaxz;
    
    thetayz = acos(*dirz / sqrt(pow(*dirz,2)+pow(*diry,2)))*180/TMath::Pi();
    if(*diry<0) thetayz = -thetayz;
    
    if(abs(thetaxz)<115&&abs(thetaxz)>65)continue;//Angle      
    if(abs(thetayz)<110&&abs(thetayz)>70)continue;//Angle

    if (ntrk % 1000 == 0 && ntrk > 0) {
      std::cout << "entry => " << ntrk / 1000 << "k tracks" << std::endl;
    }

    
    //std::cout << "Event " << *run << "," << *subrun << " " << *evt  << std::endl;
    string temp=to_string(*run)+",1,"+to_string(*evt);
    //cout<<temp.c_str()<<endl;
    fprintf(fptr,"%s\n",temp.c_str());
    
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
      if (isnan(tpdirx0[i]) || isnan(tpdiry0[i]) || isnan(tpdirz0[i])) continue;
      if(isnan(dqdx0[i]) || isinf(dqdx0[i])) continue;
      
      // sce correction

      // *** manual shift in z (only till we get this fixed in sbnd geometry)
      // ** remove this line once it's fixed
      std::string intype = cintyp;
      if(intype == "doData_dev" || intype == "doMC2024B_sub123" || intype == "doMC2024B_full") tpz0[i] = tpz0[i] - 4.2;

      
      XYZVector sp_sce_uncorr(tpx0[i], tpy0[i], tpz0[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);


      // Nu'26-style fiducial volume cut (DENTFORCE).
      // Putting this BEFORE ibinx/ibiny/ibinz and BEFORE any Fill().
      if (!PassNu26FiducialVolume(sp_sce_corr.X(),
                                  sp_sce_corr.Y(),
                                  sp_sce_corr.Z())) continue;



      int xs = GetXSectionFullRange(sp_sce_corr.X());
      if (xs < 0) continue;


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
	dqdxHist_x[0][0][xs]->Fill(dqdx_sce_corr * elife0);
  zynhits_x[0][0][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
  zyHistdqdx_x[0][0][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife0); 
  nyzHist_x[0][0][xs][ibiny][ibinz]->Fill(dqdx_sce_corr * elife0);
	//nyzHist[0][0][ibiny][ibinz]->Fill(dqdx_sce_corr);
      }
      else{
	dqdxHist[0][0]->Fill(dqdx0[i] * elife1);
	dqdxHist[0][1]->Fill(dqdx_sce_corr * elife1);
	dqdxHist_x[0][1][xs]->Fill(dqdx_sce_corr * elife1);
	zynhits_x[0][1][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
  zyHistdqdx_x[0][1][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife1); 
  nyzHist_x[0][1][xs][ibiny][ibinz]->Fill(dqdx_sce_corr * elife1);


	//nyzHist[0][1][ibiny][ibinz]->Fill(dqdx_sce_corr);
      }
      
      
    }

    // plane 1
    
    for (unsigned i = 0; i < dqdx1.GetSize(); i++) {
      if(isnan(tpx1[i])||isnan(tpy1[i])||isnan(tpz1[i]))continue;
      if (isnan(tpdirx1[i]) || isnan(tpdiry1[i]) || isnan(tpdirz1[i])) continue;
      if(isnan(dqdx1[i]) || isinf(dqdx1[i])) continue;

      // sce correction

      // *** manual shift in z (only till we get this fixed in sbnd geometry)
      // ** remove this line once it's fixed
      std::string intype = cintyp;
      if(intype == "doData_dev" || intype == "doMC2024B_sub123" || intype == "doMC2024B_full") tpz1[i] = tpz1[i] - 4.2;

      XYZVector sp_sce_uncorr(tpx1[i], tpy1[i], tpz1[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);


      // Nu'26-style fiducial volume cut (DENTFORCE).
      if (!PassNu26FiducialVolume(sp_sce_corr.X(),
                                  sp_sce_corr.Y(),
                                  sp_sce_corr.Z())) continue;

      int xs = GetXSectionFullRange(sp_sce_corr.X());
      if (xs < 0) continue;
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
	dqdxHist_x[1][0][xs]->Fill(dqdx_sce_corr * elife0);
  zynhits_x[1][0][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
  zyHistdqdx_x[1][0][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife0); 
  nyzHist_x[1][0][xs][ibiny][ibinz]->Fill(dqdx_sce_corr * elife0);
	//nyzHist[1][0][ibiny][ibinz]->Fill(dqdx_sce_corr);
      }
      else{
	dqdxHist[1][0]->Fill(dqdx1[i] * elife1);
	dqdxHist[1][1]->Fill(dqdx_sce_corr * elife1);
	dqdxHist_x[1][1][xs]->Fill(dqdx_sce_corr * elife1);
  zynhits_x[1][1][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
  zyHistdqdx_x[1][1][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife1); 
  nyzHist_x[1][1][xs][ibiny][ibinz]->Fill(dqdx_sce_corr * elife1);
	//nyzHist[1][1][ibiny][ibinz]->Fill(dqdx_sce_corr);
      }
      


    }

    // plane 2
    
    for (unsigned i = 0; i < dqdx2.GetSize(); i++) {
      if(isnan(tpx2[i])||isnan(tpy2[i])||isnan(tpz2[i]))continue;
      if (isnan(tpdirx2[i]) || isnan(tpdiry2[i]) || isnan(tpdirz2[i])) continue;
      if(isnan(dqdx2[i]) || isinf(dqdx2[i])) continue;

      if(tpx2[i]<0){
	if(tpz2[i]>58 && tpz2[i]<66) continue;   // to avoid floating point leakage
      }

      // sce correction

      // *** manual shift in z (only till we get this fixed in sbnd geometry)
      // ** remove this line once it's fixed
      std::string intype = cintyp;
      if(intype == "doData_dev" || intype == "doMC2024B_sub123" || intype == "doMC2024B_full") tpz2[i] = tpz2[i] - 4.2;


      XYZVector sp_sce_uncorr(tpx2[i], tpy2[i], tpz2[i]);
      XYZVector sp_sce_corr = sce_corr_mc->WireToTrajectoryPosition(sp_sce_uncorr);




      // Nu'26-style fiducial volume cut (DENTFORCE).
      // This filters collection-plane hits before they enter zynhits,
      // zyHistdqdx, and nyzHist, so they do not affect the median.
      if (!PassNu26FiducialVolume(sp_sce_corr.X(),
                                  sp_sce_corr.Y(),
                                  sp_sce_corr.Z())) continue;

      int xs = GetXSectionFullRange(sp_sce_corr.X());
      if (xs < 0) continue;

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
	dqdxHist_x[2][0][xs]->Fill(dqdx_sce_corr * elife0);
  zynhits_x[2][0][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
  zyHistdqdx_x[2][0][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife0); 
  nyzHist_x[2][0][xs][ibiny][ibinz]->Fill(dqdx_sce_corr * elife0);
	//nyzHist[2][0][ibiny][ibinz]->Fill(dqdx_sce_corr);
      }
      else{
	dqdxHist[2][0]->Fill(dqdx2[i] * elife1);
	dqdxHist[2][1]->Fill(dqdx_sce_corr * elife1);
	dqdxHist_x[2][1][xs]->Fill(dqdx_sce_corr * elife1);

  zynhits_x[2][1][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y());
  zyHistdqdx_x[2][1][xs]->Fill(sp_sce_corr.Z(), sp_sce_corr.Y(), dqdx_sce_corr * elife1); 
  nyzHist_x[2][1][xs][ibiny][ibinz]->Fill(dqdx_sce_corr * elife1);
	//nyzHist[2][1][ibiny][ibinz]->Fill(dqdx_sce_corr);
      }
      
      
    }
    
    //cout<<"[DEBUG:] done with the dqdx loop"<<endl;
    
  }
  cout<<"Total tracks: "<<ntrk<<endl;
  cout<<"# tracks in x>=0 region: "<<ntrkpos<<endl;
  cout<<"# tracks in x<0 region: "<<ntrkneg<<endl;
  
  const int nq=1;
  double xq[nq]={0.5}, yq[nq];

  std::vector<double> medianyz_vec_x[nplanes][2][nxsections];

  for (int l = 0; l < nplanes; l++) {
    for (int k = 0; k < 2; k++) {
      for (int xs = 0; xs < nxsections; xs++) {
        for (int i = 0; i < nbiny; i++) {
          for (int j = 0; j < nbinz; j++) {

            if (nyzHist_x[l][k][xs][i][j]->Integral() < 2) continue;

            nyzHist_x[l][k][xs][i][j]->GetQuantiles(nq, yq, xq);

            if (yq[0] < 2) continue;

            zyHist_x[l][k][xs]->SetBinContent(j+1, i+1, yq[0]);
            medianyz_vec_x[l][k][xs].push_back(yq[0]);
          }
        }
      }
    }
  }

  cout<<"[DEBUG:] Got median of medians"<<endl;


  // for (dQ/dx)_{X}^{global}

  double global_dqdx_medianyz_x[nplanes][2][nxsections];

  for (int l = 0; l < nplanes; l++) {
    for (int k = 0; k < 2; k++) {
      for (int xs = 0; xs < nxsections; xs++) {
        global_dqdx_medianyz_x[l][k][xs] = getMedian(medianyz_vec_x[l][k][xs]);

        cout << "global median plane " << l
            << " tpc " << k
            << " xsec " << xs
            << " : " << global_dqdx_medianyz_x[l][k][xs]
            << endl;
      }
    }
  }

  cout<<"[DEBUG:] Got global median"<<endl;



  TH2F *CzyHist_x[nplanes][2][nxsections];

  for (int l = 0; l < nplanes; l++) {
    for (int k = 0; k < 2; k++) {
      for (int xs = 0; xs < nxsections; xs++) {

        CzyHist_x[l][k][xs] =
          (TH2F*)zyHist_x[l][k][xs]->Clone(Form("CzyHist_%i_%i_xsec%i", l, k, xs));

        for (int i = 0; i < nbiny; i++) {
          for (int j = 0; j < nbinz; j++) {

            if (zyHist_x[l][k][xs]->GetBinContent(j+1, i+1) < 2) continue;

            double yzcorr =
              global_dqdx_medianyz_x[l][k][xs] /
              zyHist_x[l][k][xs]->GetBinContent(j+1, i+1);

            CzyHist_x[l][k][xs]->SetBinContent(j+1, i+1, yzcorr);
          }
        }

        CzyHist_x[l][k][xs]->GetXaxis()->CenterTitle();
        CzyHist_x[l][k][xs]->GetYaxis()->CenterTitle();
        CzyHist_x[l][k][xs]->GetZaxis()->CenterTitle();

        CzyHist_x[l][k][xs]->SetXTitle("z [cm]");
        CzyHist_x[l][k][xs]->SetYTitle("y [cm]");
        CzyHist_x[l][k][xs]->SetZTitle("YZ correction factor");

        zyHist_x[l][k][xs]->Write();
        CzyHist_x[l][k][xs]->Write();

        zynhits_x[l][k][xs]->Write();
        zyHistdqdx_x[l][k][xs]->Write();

        // Only write non-empty x-sliced 1D dQ/dx histograms.
        if (dqdxHist_x[l][k][xs]->GetEntries() > 0) {
          dqdxHist_x[l][k][xs]->Write();
        }
      }
    }

    dqdxHist[l][0]->Write();
    dqdxHist[l][1]->Write();
  }
  
  
  cout<<"[DEBUG:] Done writing histograms"<<endl;

  
  
  out_file->Close();
  
}




// gStyle->SetOptStat(0);
