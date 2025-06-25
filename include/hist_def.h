#ifndef HISTOGRAM_DEFINITIONS_H
#define HISTOGRAM_DEFINITIONS_H

#include <TH1F.h>
#include <TH2F.h>


const int nbinz=100, nbinx=40, nbiny=80, nbinq=200, nbint=100, nbinc=8000, nbinrr=120;
const float lowz=0, highz=500, lowx=-200, highx=200, lowy=-200, highy=200, lowq=200, highq=2000, lowt=0, hight=1.2, lowc=4000, highc=12000;
const float lowrr=0, highrr=600;

// lowc=4000, highc=12000, nbinc=80;

TH2F *thetaxyHist = new TH2F("thetaxy","thetaxy", 36, -180., 180., 36, -180., 180.);
 
TH2F *thetaHistPos = new TH2F("thetaHistPos", "x>0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);
TH2F *thetaHistNeg = new TH2F("thetaHistNeg", "x<0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);

TH2F* countHistPos = new TH2F("countHistPos", "Count of entries for Pos;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);
TH2F* countHistNeg = new TH2F("countHistNeg", "Count of entries for Neg;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);

TH1F *nyzHist[2][nbiny][nbinz];
TH2F *zyHist[2];      // this takes median values of dqdx
TH2F *zyHistdqdx[2];     // this takes all the dqdx values (not considering only the median values)
TH2F *zynhits[2];     // this takes the number of hits
TH2F *rmsHist[2];     // to get rms value of dqdx

TH1F *nxHist[nbinx];
TH1F *ncHist[nbinc]; 
TH1F *xHist = new TH1F("xmediandqdx","xmediandqdx", nbinx, lowx, highx);   // this takes median values of dqdx
TH1F *cHist = new TH1F("channel_avgdqdx","channel_avgdqdx", nbinc, lowc, highc); // for avg dqdx vs channel number
TH1F *xnhits = new TH1F("xnhits","xnhits", nbinx, lowx, highx);
TH1F *xHistdqdx = new TH1F("xHistdqdx","xHistdqdx", nbinx, lowx, highx);   // this takes all the dqdx values (not considering only the median values) (Using : Fill(x, weight))
TH1F *ndriftdqdx[nbint];
TH1F *driftdqdx = new TH1F("driftdqdx","driftdqdx", nbint, lowt, hight);


TH1F *dqdxHist[6];

void initialize_histograms(){
  
  thetaxyHist->GetXaxis()->CenterTitle();
  thetaxyHist->GetYaxis()->CenterTitle();
  thetaxyHist->SetXTitle("#theta_{xz} [deg]");
  thetaxyHist->SetYTitle("#theta_{yz} [deg]");

  
  for(int k=0;k<2;k++)for(int i=0;i<nbiny;i++)for(int j=0;j<nbinz;j++){
	nyzHist[k][i][j]=new TH1F(Form("yz%i_%i_%i",k,i,j),"",nbinq,lowq,highq);
      }


  for(int i=0;i<2;i++){
    zyHist[i]= new TH2F(Form("zy_%i",i),"zy (Median of dqdx)", nbinz, lowz, highz, nbiny, lowy, highy);
    zyHist[i]->SetXTitle("z [cm]");
    zyHist[i]->SetYTitle("y [cm]");
    zyHist[i]->SetZTitle("Median(dQ/dx [ADC/cm])");

    zyHistdqdx[i]= new TH2F(Form("zydqdx_%i",i),"zy (dqdx)", nbinz, lowz, highz, nbiny, lowy, highy);
    zyHistdqdx[i]->SetXTitle("z [cm]");
    zyHistdqdx[i]->SetYTitle("y [cm]");
    zyHistdqdx[i]->SetZTitle("dQ/dx [ADC/cm]");
    
    zynhits[i]= new TH2F(Form("zynhits_%i",i),"zy (nhits)", nbinz, lowz, highz, nbiny, lowy, highy);
    zynhits[i]->SetXTitle("z [cm]");
    zynhits[i]->SetYTitle("y [cm]");
    zynhits[i]->SetZTitle("number of hits");

    rmsHist[i] = new TH2F(Form("rmsHist_%i", i), "zy (rms of dq/dx)", nbinz, lowz, highz, nbiny, lowy, highy);
    rmsHist[i]->SetXTitle("z [cm]");
    rmsHist[i]->SetYTitle("y [cm]");
    rmsHist[i]->SetZTitle("RMS of dQ/dx [ADC/cm]");
  }

  for(int i=0;i<nbinx;i++){
    nxHist[i]=new TH1F(Form("x%i",i),"",nbinq,lowq,highq);
  }

  for(int i=0;i<nbinc;i++){
    ncHist[i]=new TH1F(Form("c%i",i),"",nbinq,lowq,highq);
  }

  for(int i=0;i<nbint;i++){
    ndriftdqdx[i]=new TH1F(Form("drift_t%i",i),"",nbinq,lowq,highq);
  }

  xHist->SetXTitle("x [cm]");
  xHist->SetYTitle("Median(dQ/dx [ADC/cm])");
  xHist->Sumw2();

  for(int i=0;i<6;i++){
    dqdxHist[i]=new TH1F(Form("dqdxHist_%i",i),"dqdx",nbinq,lowq,highq);

    dqdxHist[i]->SetXTitle("dQ/dx [ADC/cm]");
    dqdxHist[i]->SetYTitle("number of hits");

    dqdxHist[i]->GetXaxis()->CenterTitle();
    dqdxHist[i]->GetYaxis()->CenterTitle();
  }

}

#endif // HISTOGRAM_DEFINITIONS_H
