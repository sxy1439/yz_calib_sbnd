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



// SCE utils ends
