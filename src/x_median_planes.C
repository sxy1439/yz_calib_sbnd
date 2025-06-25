
#include "../include/utilities_planes.h"
#include "file_handling.h"
#include "../include/func_data_veto.h"
#include <TFile.h>

// to run :
//root [0] .L x_median.C
//root [1] x_median("doMC2023B_sub1", "output_files_mpv/yz_mc2023B_sub1.root")


// SCE utils

int num_wires = 1;
double wire_pitch = 0.3; // wire spacing

double angle_to_vertical(double plane, double x){
  
  double wireAngle = 0.;
  if(plane == 0){
    if(x < 0) wireAngle = -60.;
    else wireAngle = 60.;
  }
  else if(plane == 1){
    if(x < 0) wireAngle = 60.;
    else wireAngle = -60.;
  }
  else if(plane == 2){
    wireAngle = 0.;
  }
   
  double angleToVert = wireAngle * TMath::Pi()/180.;
  return angleToVert;
}

double cos_gamma_angle(double plane, double x, double dirY, double dirZ){
  double cos_gamma = std::abs(std::sin(angle_to_vertical(plane, x)) * dirY + std::cos(angle_to_vertical(plane, x)) * dirZ);
  return cos_gamma;
}


TFile* fileHists = new TFile("SCE_hists_E500.root");

TH3F* TrueFwdX;
TH3F* TrueFwdY;
TH3F* TrueFwdZ;

TH3F* TrueBkwdX;
TH3F* TrueBkwdY;
TH3F* TrueBkwdZ;

void LoadSCEMaps()
{
  TrueFwdX = (TH3F*) fileHists->Get("TrueFwd_Displacement_X");
  TrueFwdY = (TH3F*) fileHists->Get("TrueFwd_Displacement_Y");
  TrueFwdZ = (TH3F*) fileHists->Get("TrueFwd_Displacement_Z");

  TrueBkwdX = (TH3F*) fileHists->Get("TrueBkwd_Displacement_X");
  TrueBkwdY = (TH3F*) fileHists->Get("TrueBkwd_Displacement_Y");
  TrueBkwdZ = (TH3F*) fileHists->Get("TrueBkwd_Displacement_Z");
}



double GetFwdOffset(double x_true, double y_true, double z_true, int comp)
{
  if (x_true < -199.999) x_true = -199.999;
  else if (x_true > 199.999) x_true = 199.999;
  
  if (y_true < -199.999) y_true = -199.999;
  else if (y_true > 199.999) y_true = 199.999;
  
  if (z_true < 0.001) z_true = 0.001;
  else if (z_true > 499.999) z_true = 499.999;

  if(x_true < -199.999 || x_true > 199.999 || y_true < -199.999 || y_true > 199.999 || z_true < 0.001 || z_true > 499.999){
    cout<<"x_true : "<<x_true<<" y_true : "<<y_true<<" z_true : "<<z_true<<endl;
  }
  
  double offset = 0.0;
  if(comp == 1){
    offset = TrueFwdX->Interpolate(x_true, y_true, z_true);
  }
  else if(comp == 2){
    offset = TrueFwdY->Interpolate(x_true, y_true, z_true);
  }
  else if(comp == 3){
    offset = TrueFwdZ->Interpolate(x_true, y_true, z_true);
  }
  
  return offset;
}

double GetBkwdOffset(double x_reco, double y_reco, double z_reco, int comp)
{
 
  if (x_reco < -199.999) x_reco = -199.999;
  else if (x_reco > 199.999) x_reco = 199.999;
  
  if (y_reco < -199.999) y_reco = -199.999;
  else if (y_reco > 199.999) y_reco = 199.999;
  
  if (z_reco < 0.001) z_reco = 0.001;
  else if (z_reco > 499.999) z_reco = 499.999;

  if(x_reco < -199.999 || x_reco > 199.999 || y_reco < -199.999 || y_reco > 199.999 || z_reco < 0.001 || z_reco > 499.999){
    cout<<"x_reco : "<<x_reco<<" y_reco : "<<y_reco<<" z_reco : "<<z_reco<<endl;
  }
  
  double offset = 0.0;
  if(comp == 1){
    offset = TrueBkwdX->Interpolate(x_reco,y_reco,z_reco);
  }
  else if(comp == 2){
    offset = TrueBkwdY->Interpolate(x_reco,y_reco,z_reco);
  }
  else if(comp == 3){
    offset = TrueBkwdZ->Interpolate(x_reco,y_reco,z_reco);
  }

  return offset;
}



//std::vector<double> loc_vector(double x, double y, double z){
TVector3 loc_vector(double x, double y, double z){

  TVector3 loc;
  loc = TVector3(x, y, z);

  if(loc(0) < 0.0)
    {
      loc(0) *= -1.0;
      loc(1) *= -1.0;
      loc(2) *= -1.0;
    }

  return loc;
  
}

TVector3 dir_vector(double x, double y, double z, double dirX, double dirY, double dirZ){

  TVector3 loc;
  TVector3 vec;
  
  loc = TVector3(x, y, z);
  vec = TVector3(dirX, dirY, dirZ);
  
  if(loc(0) < 0.0)
    {
      vec(0) *= -1.0;
      vec(1) *= -1.0;
      vec(2) *= -1.0;
    }
  
  TVector3 unit = vec.Unit();
  
  return unit;
}

// get_nominal_spatial_position(loc_vector(x, y, z), dir_vector(x, y, z, dirX, dirY, dirZ))

TVector3 get_nominal_spatial_position(const TVector3& loc, const TVector3& unit)
{
 
  TVector3 nextloc = loc + ((((double) num_wires) * wire_pitch) / unit(2)) * unit;
  TVector3 sizevec_nominal = nextloc - loc;

  return sizevec_nominal;
}

/*
TVector3 get_sce_corrected_spatial_position(const TVector3& loc, const TVector3& unit)
{

  TVector3 nextloc = loc + ((((double) num_wires) * wire_pitch) / unit(2)) * unit;
  
  TVector3 loc_modified = TVector3(loc(0) + GetBkwdOffset(loc(0), loc(1), loc(2), 1),
				   loc(1) + GetBkwdOffset(loc(0), loc(1), loc(2), 2),
				   loc(2) + GetBkwdOffset(loc(0), loc(1), loc(2), 3) );
  
  TVector3 nextloc_modified = TVector3(nextloc(0) + GetBkwdOffset(nextloc(0), nextloc(1), nextloc(2), 1),
				       nextloc(1) + GetBkwdOffset(nextloc(0), nextloc(1), nextloc(2), 2),
				       nextloc(2) + GetBkwdOffset(nextloc(0), nextloc(1), nextloc(2), 3) );
  
  TVector3 sizevec_modified = nextloc_modified - loc_modified;
  
  return sizevec_modified;
}
*/



TVector3 get_sce_corrected_spatial_position(const TVector3& loc, const TVector3& unit, bool isData)
{

  TVector3 nextloc = loc + ((((double) num_wires) * wire_pitch) / unit(2)) * unit;

  TVector3 loc_modified, nextloc_modified;
  
  if(isData==true){
    
    loc_modified = TVector3(loc(0) + GetBkwdOffset(loc(0), loc(1), loc(2), 1),
			    loc(1) + GetBkwdOffset(loc(0), loc(1), loc(2), 2),
			    loc(2) + GetBkwdOffset(loc(0), loc(1), loc(2), 3) );
    
    nextloc_modified = TVector3(nextloc(0) + GetBkwdOffset(nextloc(0), nextloc(1), nextloc(2), 1),
				nextloc(1) + GetBkwdOffset(nextloc(0), nextloc(1), nextloc(2), 2),
				nextloc(2) + GetBkwdOffset(nextloc(0), nextloc(1), nextloc(2), 3) );
  }
  else if(isData==false){
    loc_modified = TVector3(loc(0) + GetFwdOffset(loc(0), loc(1), loc(2), 1),
			    loc(1) + GetFwdOffset(loc(0), loc(1), loc(2), 2),
			    loc(2) + GetFwdOffset(loc(0), loc(1), loc(2), 3) );
    
    nextloc_modified = TVector3(nextloc(0) + GetFwdOffset(nextloc(0), nextloc(1), nextloc(2), 1),
				nextloc(1) + GetFwdOffset(nextloc(0), nextloc(1), nextloc(2), 2),
				nextloc(2) + GetFwdOffset(nextloc(0), nextloc(1), nextloc(2), 3) );
  }
  
  TVector3 sizevec_modified = nextloc_modified - loc_modified;
  
  return sizevec_modified;
}



// SCE utils ends



void x_median_planes(const char* cintyp, const char* inYZOut) {

  LoadSCEMaps();

  TH2F *CzyHist[nplanes][2];
  TFile* file_YZOut = new TFile(inYZOut, "READ");
  for(int l=0;l<nplanes;l++){
    for(int k=0;k<2;k++){
      CzyHist[l][k] = (TH2F*)file_YZOut->Get(Form("CzyHist_%i_%i",l,k));
    }
  }

  
  gSystem->Load("libFileHandling.so");
  if (gSystem->Load("libFileHandling.so") < 0) {
    std::cerr << "Error loading library" << std::endl;
    return;
  }

  auto e_callife = getELifetime(cintyp);
  
  TString output_name = getOutputNameX(cintyp);
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

  /*
  TTreeReaderArray<float> dqdx(myReader, "trk.hits2.dqdx"); // hits on plane 2 (Collection)
  TTreeReaderArray<float> tpx(myReader, "trk.hits2.h.sp.x"); // x of track trajectory position  (older was -> trk.hits2.tp.x)
  TTreeReaderArray<float> tpy(myReader, "trk.hits2.h.sp.y"); // y of track trajectory position  
  TTreeReaderArray<float> tpz(myReader, "trk.hits2.h.sp.z"); // z of track trajectory position
  TTreeReaderArray<float> time(myReader, "trk.hits2.h.time"); // in ticks (500 ns), up to 3200 (1.6 ms > 1.25 ms)
  */

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
  

  cout<<"nbinx = "<<nbinx<<endl;
  cout<<"elife[0] = "<<e_callife[0]<<" elife[1] = "<<e_callife[0]<<endl;

  /// X equalization

  int ntrk=0;
  int ibinx, ibiny, ibinz;

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

    for(unsigned i = 0; i < dqdx0.GetSize(); i++){
      if(isnan(tpx0[i])||isnan(tpy0[i])||isnan(tpz0[i]))continue;
      if(isnan(dqdx0[i]) || isinf(dqdx0[i])) continue;

      ibinx = floor((tpx0[i]-lowx)/(highx-lowx)*nbinx);
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((tpy0[i]-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((tpz0[i]-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;

      TVector3 dir_sce_corrected;   
      if(string(cintyp).find("Data") != string::npos){
	dir_sce_corrected = get_sce_corrected_spatial_position(loc_vector(tpx0[i], tpy0[i], tpz0[i]), dir_vector(tpx0[i], tpy0[i], tpz0[i], tpdirx0[i], tpdiry0[i], tpdirz0[i]), true).Unit();
      }
      else {
	dir_sce_corrected = get_sce_corrected_spatial_position(loc_vector(tpx0[i], tpy0[i], tpz0[i]), dir_vector(tpx0[i], tpy0[i], tpz0[i], tpdirx0[i], tpdiry0[i], tpdirz0[i]), false).Unit();
      }

      double pitch_sce_corrected=0.;
      double cosgamma_sce_corrected = cos_gamma_angle(0, tpx0[i], dir_sce_corrected.Y(), dir_sce_corrected.Z());
      if(cosgamma_sce_corrected) pitch_sce_corrected = wire_pitch / cosgamma_sce_corrected;
      double dqdx_sce_corrected = charge0[i]/pitch_sce_corrected;

      double elife0 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time0[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      
      if(tpx0[i]<0){
	xnhits[0]->Fill(tpx0[i]);
	xHistdqdx[0]->Fill(tpx0[i], dqdx_sce_corrected * elife0);
	nxHist[0][ibinx]->Fill(dqdx_sce_corrected * elife0 * CzyHist[0][0]->GetBinContent(ibinz+1, ibiny+1));
      }
      else{
	xnhits[0]->Fill(tpx0[i]);
	xHistdqdx[0]->Fill(tpx0[i], dqdx_sce_corrected * elife1);
	nxHist[0][ibinx]->Fill(dqdx_sce_corrected * elife1 * CzyHist[0][1]->GetBinContent(ibinz+1, ibiny+1));
      }
      
    }

    // plane 1

    for(unsigned i = 0; i < dqdx1.GetSize(); i++){
      if(isnan(tpx1[i])||isnan(tpy1[i])||isnan(tpz1[i]))continue;
      if(isnan(dqdx1[i]) || isinf(dqdx1[i])) continue;

      ibinx = floor((tpx1[i]-lowx)/(highx-lowx)*nbinx);
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((tpy1[i]-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((tpz1[i]-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;

      TVector3 dir_sce_corrected;   
      if(string(cintyp).find("Data") != string::npos){
	dir_sce_corrected = get_sce_corrected_spatial_position(loc_vector(tpx1[i], tpy1[i], tpz1[i]), dir_vector(tpx1[i], tpy1[i], tpz1[i], tpdirx1[i], tpdiry1[i], tpdirz1[i]), true).Unit();
      }
      else {
	dir_sce_corrected = get_sce_corrected_spatial_position(loc_vector(tpx1[i], tpy1[i], tpz1[i]), dir_vector(tpx1[i], tpy1[i], tpz1[i], tpdirx1[i], tpdiry1[i], tpdirz1[i]), false).Unit();
      }

      double pitch_sce_corrected=0.;
      double cosgamma_sce_corrected = cos_gamma_angle(1, tpx1[i], dir_sce_corrected.Y(), dir_sce_corrected.Z());
      if(cosgamma_sce_corrected) pitch_sce_corrected = wire_pitch / cosgamma_sce_corrected;
      double dqdx_sce_corrected = charge1[i]/pitch_sce_corrected;

      double elife0 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time1[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      
      if(tpx1[i]<0){
	xnhits[1]->Fill(tpx1[i]);
	xHistdqdx[1]->Fill(tpx1[i], dqdx_sce_corrected * elife0);
	nxHist[1][ibinx]->Fill(dqdx_sce_corrected * elife0 * CzyHist[1][0]->GetBinContent(ibinz+1, ibiny+1));
      }
      else{
	xnhits[1]->Fill(tpx1[i]);
	xHistdqdx[1]->Fill(tpx1[i], dqdx_sce_corrected * elife1);
	nxHist[1][ibinx]->Fill(dqdx_sce_corrected * elife1 * CzyHist[1][1]->GetBinContent(ibinz+1, ibiny+1));
      }
      
    }

    // plane 2

    for(unsigned i = 0; i < dqdx2.GetSize(); i++){
      if(isnan(tpx2[i])||isnan(tpy2[i])||isnan(tpz2[i]))continue;
      if(isnan(dqdx2[i]) || isinf(dqdx2[i])) continue;

      ibinx = floor((tpx2[i]-lowx)/(highx-lowx)*nbinx);
      if(ibinx<0||ibinx>=nbinx)continue;
      ibiny = floor((tpy2[i]-lowy)/(highy-lowy)*nbiny);
      if(ibiny<0||ibiny>=nbiny)continue;
      ibinz = floor((tpz2[i]-lowz)/(highz-lowz)*nbinz);
      if(ibinz<0||ibinz>=nbinz)continue;

      TVector3 dir_sce_corrected;   
      if(string(cintyp).find("Data") != string::npos){
	dir_sce_corrected = get_sce_corrected_spatial_position(loc_vector(tpx2[i], tpy2[i], tpz2[i]), dir_vector(tpx2[i], tpy2[i], tpz2[i], tpdirx2[i], tpdiry2[i], tpdirz2[i]), true).Unit();
      }
      else {
	dir_sce_corrected = get_sce_corrected_spatial_position(loc_vector(tpx2[i], tpy2[i], tpz2[i]), dir_vector(tpx2[i], tpy2[i], tpz2[i], tpdirx2[i], tpdiry2[i], tpdirz2[i]), false).Unit();
      }

      double pitch_sce_corrected=0.;
      double cosgamma_sce_corrected = cos_gamma_angle(2, tpx2[i], dir_sce_corrected.Y(), dir_sce_corrected.Z());
      if(cosgamma_sce_corrected) pitch_sce_corrected = wire_pitch / cosgamma_sce_corrected;
      double dqdx_sce_corrected = charge2[i]/pitch_sce_corrected;

      double elife0 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[0]);
      double elife1 = lifetime_correction(time2[i], trk_t0, ticksToMs, preTriggerWindow, nsToMs, e_callife[1]);

      
      if(tpx2[i]<0){
	xnhits[2]->Fill(tpx2[i]);
	xHistdqdx[2]->Fill(tpx2[i], dqdx_sce_corrected * elife0);
	nxHist[2][ibinx]->Fill(dqdx_sce_corrected * elife0 * CzyHist[2][0]->GetBinContent(ibinz+1, ibiny+1));
      }
      else{
	xnhits[2]->Fill(tpx2[i]);
	xHistdqdx[2]->Fill(tpx2[i], dqdx_sce_corrected * elife1);
	nxHist[2][ibinx]->Fill(dqdx_sce_corrected * elife1 * CzyHist[2][1]->GetBinContent(ibinz+1, ibiny+1));
      }
      
    }
  }

  const int nq=1;
  double xq[nq]={0.5}, yq[nq];
  std::vector<double> medianx_vec[nplanes];
  for(int l=0;l<nplanes;l++){
    for(int i=0;i<nbinx;i++){
      if(nxHist[l][i]->Integral()<2)continue;
      nxHist[l][i]->GetQuantiles(nq, yq, xq);
      xHist[l]->SetBinContent(i+1, yq[0]);       /// (xHist->SetBinContent(i+1, yq[0]) = (dQ/dx)_{X}^{local})
      
      medianx_vec[l].push_back(yq[0]); 
    }
  }

  /*
  for(int i=0;i<nbinx;i++){
    nxHist[i]->Write();
  }
  */

  // for (dQ/dx)_{X}^{global}
  double global_dqdx_medianx[nplanes];
  for(int l=0;l<nplanes;l++){
    global_dqdx_medianx[l] = getMedian(medianx_vec[l]);
    cout<< "global median : "<<global_dqdx_medianx[l]<<endl;
  }
  
  
  // for correction factor
  TH2F *CxHist[nplanes];
  for(int l=0;l<nplanes;l++){
    CxHist[l] = (TH2F*)xHist[l]->Clone(Form("CxHist_%i", l));
    for(int i=0;i<nbinx;i++){
      if(xHist[l]->GetBinContent(i+1)<1e-2)continue;
      CxHist[l]->SetBinContent(i+1, global_dqdx_medianx[l] / xHist[l]->GetBinContent(i+1));
    }

    CxHist[l]->GetXaxis()->CenterTitle();
    CxHist[l]->GetYaxis()->CenterTitle();
    
    CxHist[l]->SetXTitle("x [cm]");
    CxHist[l]->SetYTitle("X correction factor");
    
    xnhits[l]->Write();      // stores number of hits
    xHist[l]->Write();       // median dQdx
    xHistdqdx[l]->Write();   // dQdx values (not the median)
    CxHist[l]->Write();       // correction factor
  }
  


  
  out_file->Close();
  
}
