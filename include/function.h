bool Is_Edge(float x, float y, float z){   // FV
  bool res = false;
  if(abs(x)>190) res = true;
  //if(abs(y)>160) res = true;
  if(abs(y)>190) res = true;
  if(z<10 || z>490) res = true;
  return res;
}

bool Is_EdgeYZ(float y, float z){  
  bool res = false;
  if(abs(y)>195) res = true;
  if(z<5 || z>495) res = true;
  return res;
}

bool Is_Cathode_Crossing(float startx, float endx){  //start and end points on opposite sides of the cathode
  if ((startx < 0 && endx > 0) || (startx > 0 && endx < 0)) return true;     
  return false;
}

float getMedian(std::vector<double>& values){
  if (!values.empty()){
    size_t size = values.size();
    std::sort(values.begin(), values.end());
    if (size % 2 == 0) {
      return (values[size / 2 - 1] + values[size / 2]) / 2;
    } else {
      return values[size / 2];
    }
  }
  return 0.;
}

float getMean(std::vector<double>& values){
  size_t size = values.size();
  if (size == 0) return 0.;
  
  double sum=0.;
  for(int i=0; i<size; i++){
    sum += values[i];
  }
  return sum/size;
}

float getFWHM(TH1F* hist){
  if (!hist) {
    cout << "Histogram not found!" << endl;
    return -1;
  }
  
  double maxY = hist->GetMaximum();  
  double half_max = maxY / 2.;       
  
  int bin1 = -1, bin2 = -1, bin_peak = -1;
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    double y = hist->GetBinContent(i);
    
    if (y >= half_max && bin1 == -1) { 
      bin1 = i; 
    }
    if (y >= half_max) { // checks when bin1 will be filled
      bin2 = i;  
    }
  }
  
  if (bin1 == -1 || bin2 == -1) {
    cout << "couldn't determine FWHM for " << hist->GetName() << endl;
    return -1;
  }

  bin_peak = hist->GetMaximumBin();
  
  double x1 = hist->GetBinCenter(bin1);
  double x2 = hist->GetBinCenter(bin2);

  cout<< "max hits (y-axis) : "<< maxY << "; half max hits (y-axis) : "<< half_max <<endl;

  cout<< " bin at max dqdx : "<<bin_peak << "; dQ/dx at peak : "<< hist->GetBinCenter(bin_peak) <<endl;
  cout<< " bin at fwhm start : "<< bin1 << "; corresponding dQ/dx : "<< hist->GetBinCenter(bin1) <<endl;
  cout<< " bin at fwhm end : "<< bin2 << "; corresponding dQ/dx : "<< hist->GetBinCenter(bin2) <<endl;
  
  cout<< "hits (y-axis) just before the fwhm start : "<< hist->GetBinContent(bin1-1) <<endl;
  cout<< "hits (y-axis) just after the fwhm end : "<< hist->GetBinContent(bin2+1) <<endl;
  
  cout<< "location 1 : "<<x1<<endl;
  cout<< "location 2 : "<<x2<<endl;
  cout<< "width : "<< x2 - x1 <<endl;

  
  return x2 - x1;  // FWHM = x2 - x1
}

float getFWHMresolution(TH1F* hist){
  if (!hist) {
    cout << "Histogram not found!" << endl;
    return -1;
  }
  
  double maxY = hist->GetMaximum();  
  double half_max = maxY / 2.;       
  
  int bin1 = -1, bin2 = -1, bin_peak = -1;
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    double y = hist->GetBinContent(i);
 
    if (y >= half_max && bin1 == -1) { 
      bin1 = i; 
    }
    if (y >= half_max) { 
      bin2 = i;  
    }
  }
  
  if (bin1 == -1 || bin2 == -1) {
    cout << "couldn't determine FWHM for " << hist->GetName() << endl;
    return -1;
  }

  bin_peak = hist->GetMaximumBin();
  
  double x1 = hist->GetBinCenter(bin1);
  double x2 = hist->GetBinCenter(bin2);

  cout<< " bin at max dqdx : "<<bin_peak << "; dQ/dx at peak : "<< hist->GetBinCenter(bin_peak) <<endl;
  cout<< " bin at fwhm start : "<< bin1 << "; corresponding dQ/dx : "<< hist->GetBinCenter(bin1) <<endl;
  cout<< " bin at fwhm end : "<< bin2 << "; corresponding dQ/dx : "<< hist->GetBinCenter(bin2) <<endl;
  
  cout<< "hits (y-axis) just before the fwhm start : "<< hist->GetBinContent(bin1-1) <<endl;
  cout<< "hits (y-axis) just after the fwhm end : "<< hist->GetBinContent(bin2+1) <<endl;
  
  cout<< "location 1 : "<<x1<<endl;
  cout<< "location 2 : "<<x2<<endl;
  cout<< "width : "<< x2 - x1 <<endl;

  cout<< "max hits (y-axis) : "<< maxY << "; half max hits (y-axis) : "<< half_max <<endl;
  cout<< "width : "<< x2 - x1 <<endl;

  double resolution = -1.;
  double dqdx_max = hist->GetBinCenter(bin_peak);
  
  if(maxY>0) resolution = ((x2 - x1)/dqdx_max)*100;

  
  return resolution;  
}

/*
// gaussian fit

float getFWHMfit(TH1F* hist){
  if (!hist) {
    cout << "Histogram not found!" << endl;
    return -1;
  }

  TF1 *fgaus = new TF1("fgaus","gaus",950, 1100);

  hist->Fit("fgaus", "R");

  double sd = fgaus->GetParameter(2);
  double fwhm = 2.355*sd;

  cout<<"SD from the fit : "<<sd<<endl;
  cout<<"FWHM from the fit : "<<fwhm<<endl;
 
  
  return fwhm;  
}
*/

// LG fit

float getFWHMfit(TH1F* hist) {
  if (!hist) {
    cout << "Histogram not found!" << endl;
    return -1;
  }
  
  TF1 *fitf = new TF1("fitf", langaufun, 950, 1250, 4);
  fitf->SetParameter(0, 40);      
  fitf->SetParameter(1, 1000);    
  fitf->SetParameter(2, 1.5e5);  
  fitf->SetParameter(3, 50);      
  
  fitf->SetRange(950, 1250);
  
  hist->Fit(fitf, "R"); 
  
  
  double mpv = fitf->GetParameter(1);
  double maxY = fitf->Eval(mpv);
  double halfMax = maxY/2.;
  
  double xLow = fitf->GetX(halfMax, 950, mpv);  
  double xHigh = fitf->GetX(halfMax, mpv, 1250); 
  
  double fwhm = xHigh - xLow;
  
  cout << "mpv from the fit: " << mpv << endl;
  cout << "FWHM from the fit: " << fwhm << endl;
  
  return fwhm;
}

float getFWHMfitres(TH1F* hist) {
  if (!hist) {
    cout << "Histogram not found!" << endl;
    return -1;
  }
  
  TF1 *fitf = new TF1("fitf", langaufun, 950, 1250, 4);
  fitf->SetParameter(0, 40);      
  fitf->SetParameter(1, 1000);    
  fitf->SetParameter(2, 1.5e5);  
  fitf->SetParameter(3, 50);      
  
  fitf->SetRange(950, 1250);
  
  hist->Fit(fitf, "R"); 
  
  
  double mpv = fitf->GetParameter(1);
  double maxY = fitf->Eval(mpv);
  double halfMax = maxY/2.;
  
  double xLow = fitf->GetX(halfMax, 950, mpv);  
  double xHigh = fitf->GetX(halfMax, mpv, 1250); 
  
  double fwhm = xHigh - xLow;

  double bin_peak = hist->GetMaximumBin();
  
  cout << "mpv from the fit: " << mpv << endl;
  cout << "FWHM from the fit: " << fwhm << endl;

  double resolution = -1.;
  double dqdx_max = hist->GetBinCenter(bin_peak); // need to change this
  
  if(maxY>0) resolution = (fwhm/dqdx_max)*100;
  
  return resolution;
}




/// formula to get the drift time
// float drift_time = (time[i] * ticksToMs) - preTriggerWindow - (t0 * nsToMs);

float drift_time(float ti, float t0, float ticksToMs, float preTriggerWindow, float nsToMs){
  float tdrift = (ti * ticksToMs) - preTriggerWindow - (t0 * nsToMs);
  return tdrift;
}


float lifetime_correction(float ti, float t0, float ticksToMs, float preTriggerWindow, float nsToMs, float tau){
  float correction_factor = 1.;
  float tdrift = (ti * ticksToMs) - preTriggerWindow - (t0 * nsToMs);

  correction_factor = exp(tdrift/2/tau);

  return correction_factor;
}
