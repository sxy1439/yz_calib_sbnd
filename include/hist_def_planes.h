#ifndef HISTOGRAM_DEFINITIONS_H
#define HISTOGRAM_DEFINITIONS_H

#include <TH1F.h>
#include <TH2F.h>

const int nplanes = 3;
const int nbinz=100, nbinx=40, nbiny=80, nbinq=200, nbint=100, nbinc=8000, nbinang=100;
const float lowz=0, highz=500, lowx=-200, highx=200, lowy=-200, highy=200, lowq=200, highq=2000, lowt=0, hight=1.2, lowc=4000, highc=12000, lowang=-180, highang=180;

//const int nbinz=100, nbinx=40, nbiny=80, nbinq=200, nbint=100, nbinc=8000, nbinang=100;
//const float lowz=0, highz=500, lowx=-200, highx=200, lowy=-200, highy=200, lowq=200*50, highq=2000*50, lowt=0, hight=1.2, lowc=4000, highc=12000, lowang=-180, highang=180;

//const int nbinz=600, nbinx=40, nbiny=80, nbinq=200, nbint=100, nbinc=8000, nbinang=100;
//const float lowz=-50, highz=550, lowx=-200, highx=200, lowy=-200, highy=200, lowq=200, highq=2000, lowt=0, hight=1.2, lowc=4000, highc=12000, lowang=-180, highang=180;



// lowc=4000, highc=12000, nbinc=80;

TH2F *thetaxyHist = new TH2F("thetaxy","thetaxy", 36, -180., 180., 36, -180., 180.);
 
TH2F *thetaHistPos = new TH2F("thetaHistPos", "x>0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);
TH2F *thetaHistNeg = new TH2F("thetaHistNeg", "x<0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);

TH2F* countHistPos = new TH2F("countHistPos", "Count of entries for Pos;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);
TH2F* countHistNeg = new TH2F("countHistNeg", "Count of entries for Neg;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180., 100, -180., 180.);

TH1F *thetaHistPos1D[nplanes];
TH1F *thetaHistNeg1D[nplanes];
TH1F* countHistPos1D[nplanes];
TH1F* countHistNeg1D[nplanes];
TH1F *thetaHistPos1D_med[nplanes];
TH1F *thetaHistNeg1D_med[nplanes];


TH1F *nyzHist[nplanes][2][nbiny][nbinz];
TH2F *zyHist[nplanes][2];      // this takes median values of dqdx
TH2F *zyHistdqdx[nplanes][2];     // this takes all the dqdx values (not considering only the median values)
TH2F *zynhits[nplanes][2];     // this takes the number of hits
TH2F *rmsHist[nplanes][2];     // to get rms value of dqdx

TH1F *nxHist[nplanes][nbinx];
TH1F *xHist[nplanes];
TH1F *xnhits[nplanes];
TH1F *xHistdqdx[nplanes];

TH1F *h_trkstartx_beforeAngCut = new TH1F("trkstartx_beforeAngCut"," ; trk start x [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkstarty_beforeAngCut = new TH1F("trkstarty_beforeAngCut"," ; trk start y [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkstartz_beforeAngCut = new TH1F("trkstartz_beforeAngCut"," ; trk start z [cm]; number of tracks", nbinz, lowz, highz);
TH1F *h_trkendx_beforeAngCut = new TH1F("trkendx_beforeAngCut"," ; trk end x [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkendy_beforeAngCut = new TH1F("trkendy_beforeAngCut"," ; trk end y [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkendz_beforeAngCut = new TH1F("trkendz_beforeAngCut"," ; trk end z [cm]; number of tracks", nbinz, lowz, highz);

TH1F *h_trkstartx_afterAngCut = new TH1F("trkstartx_afterAngCut"," ; trk start x [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkstarty_afterAngCut = new TH1F("trkstarty_afterAngCut"," ; trk start y [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkstartz_afterAngCut = new TH1F("trkstartz_afterAngCut"," ; trk start z [cm]; number of tracks", nbinz, lowz, highz);
TH1F *h_trkendx_afterAngCut = new TH1F("trkendx_afterAngCut"," ; trk end x [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkendy_afterAngCut = new TH1F("trkendy_afterAngCut"," ; trk end y [cm]; number of tracks", 300, -300, 300);
TH1F *h_trkendz_afterAngCut = new TH1F("trkendz_afterAngCut"," ; trk end z [cm]; number of tracks", nbinz, lowz, highz);

TH1F *h_dqdx_beforeAngCut = new TH1F("dqdx_beforeAngCut"," ; dQ/dx [ADC/cm]; number of hits", 400, 0, highq);
TH1F *h_dqdx_afterAngCut = new TH1F("dqdx_afterAngCut"," ; dQ/dx [ADC/cm]; number of hits", 400, 0, highq);




TH1F *h_hit_thetaxz[nplanes][2];
TH1F *h_hit_mediandqdx_thetaxz[nplanes][2];

TH1F *ncHist[nplanes][nbinc];
TH1F *nangHist[nplanes][2][nbinang];
TH1F *cHist[nplanes];
//TH1F *cHist = new TH1F("avgdqdx","avgdqdx", nbinc, lowc, highc);

TH1F *ndriftdqdx[nbint];
TH1F *driftdqdx = new TH1F("driftdqdx","driftdqdx", nbint, lowt, hight);


TH1F *dqdxHist[nplanes][8];

void initialize_histograms(){
  
  thetaxyHist->GetXaxis()->CenterTitle();
  thetaxyHist->GetYaxis()->CenterTitle();
  thetaxyHist->SetXTitle("#theta_{xz} [deg]");
  thetaxyHist->SetYTitle("#theta_{yz} [deg]");

  
  for(int l=0;l<nplanes;l++)for(int k=0;k<2;k++)for(int i=0;i<nbiny;i++)for(int j=0;j<nbinz;j++){
	  nyzHist[l][k][i][j]=new TH1F(Form("yz%i_%i_%i_%i",l,k,i,j),"",nbinq,lowq,highq);
	}


  for(int i=0;i<nplanes;i++)for(int j=0;j<2;j++){
      zyHist[i][j]= new TH2F(Form("zy_%i_%i",i,j),"zy (Median of dqdx)", nbinz, lowz, highz, nbiny, lowy, highy);
      zyHist[i][j]->SetXTitle("z [cm]");
      zyHist[i][j]->SetYTitle("y [cm]");
      zyHist[i][j]->SetZTitle("Median(dQ/dx [ADC/cm])");
      
      zyHistdqdx[i][j]= new TH2F(Form("zydqdx_%i_%i",i,j),"zy (dqdx)", nbinz, lowz, highz, nbiny, lowy, highy);
      zyHistdqdx[i][j]->SetXTitle("z [cm]");
      zyHistdqdx[i][j]->SetYTitle("y [cm]");
      zyHistdqdx[i][j]->SetZTitle("dQ/dx [ADC/cm]");
      
      zynhits[i][j]= new TH2F(Form("zynhits_%i_%i",i,j),"zy (nhits)", nbinz, lowz, highz, nbiny, lowy, highy);
      zynhits[i][j]->SetXTitle("z [cm]");
      zynhits[i][j]->SetYTitle("y [cm]");
      zynhits[i][j]->SetZTitle("number of hits");
      
      rmsHist[i][j] = new TH2F(Form("rmsHist_%i_%i",i,j), "zy (rms of dq/dx)", nbinz, lowz, highz, nbiny, lowy, highy);
      rmsHist[i][j]->SetXTitle("z [cm]");
      rmsHist[i][j]->SetYTitle("y [cm]");
      rmsHist[i][j]->SetZTitle("RMS of dQ/dx [ADC/cm]");
    }


  for(int i=0;i<nplanes;i++){
    thetaHistPos1D[i] = new TH1F(Form("thetaHistPos_%i", i), "avgdqdx x>0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180.);
    thetaHistNeg1D[i] = new TH1F(Form("thetaHistNeg_%i", i), "avgdqdx x<0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180.);
    countHistPos1D[i] = new TH1F(Form("countHistPos_%i", i), "Count of entries for Pos;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180.);
    countHistNeg1D[i] = new TH1F(Form("countHistNeg_%i", i), "Count of entries for Neg;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180.);

    thetaHistPos1D_med[i] = new TH1F(Form("thetaHistPos_med_%i", i), "mediandqdx x>0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180.);
    thetaHistNeg1D_med[i] = new TH1F(Form("thetaHistNeg_med_%i", i), "mediandqdx x<0;#theta_{xz} [deg];#theta_{yz} [deg]", 100, -180., 180.);
  }

  
  for(int i=0;i<nplanes;i++){
    for(int j=0;j<nbinx;j++){
      nxHist[i][j]=new TH1F(Form("x_%i_%i",i,j),"",nbinq,lowq,highq);
    }
  }

  for(int i=0;i<nplanes;i++){
    for(int j=0;j<nbinc;j++){
      ncHist[i][j]=new TH1F(Form("c_%i_%i",i,j),"",nbinq,lowq,highq);
    }
  }

  for(int i=0;i<nplanes;i++)for(int k=0;k<2;k++){
      for(int j=0;j<nbinang;j++){
	nangHist[i][k][j]=new TH1F(Form("ang_%i_%i_%i",i,k,j),"",nbinang,lowang,highang);
      }
    }

  for(int i=0;i<nplanes;i++){
    cHist[i] = new TH1F(Form("mediandqdx_%i", i),"avgdqdx", nbinc, lowc, highc);
  }
  
  for(int i=0;i<nplanes;i++){
    xHist[i] = new TH1F(Form("xmediandqdx_%i", i),"xmediandqdx", nbinx, lowx, highx);
    xnhits[i] = new TH1F(Form("xnhits_%i", i),"xnhits", nbinx, lowx, highx);
    xHistdqdx[i] = new TH1F(Form("xHistdqdx_%i", i),"xHistdqdx", nbinx, lowx, highx);

    xHist[i]->SetXTitle("x [cm]");
    xHist[i]->SetYTitle("Median(dQ/dx [ADC/cm])");
  }


  for(int i=0;i<nplanes;i++)for(int j=0;j<2;j++){
      h_hit_thetaxz[i][j] = new TH1F(Form("hhit_thetaxz_%i_%i", i,j),"dqdx wrt theta_xz", 100, -180., 180.);
    }

  for(int i=0;i<nplanes;i++)for(int k=0;k<2;k++){
      h_hit_mediandqdx_thetaxz[i][k] = new TH1F(Form("hhit_mediandqdx_thetaxz_%i_%i", i,k),"mediandqdx wrt theta_xz", 100, -180., 180.);
    }

  
  for(int i=0;i<nbint;i++){
    ndriftdqdx[i]=new TH1F(Form("drift_t%i",i),"",nbinq,lowq,highq);
  }


  for(int i=0;i<nplanes;i++){
    for(int j=0;j<8;j++){
      dqdxHist[i][j]=new TH1F(Form("dqdxHist_%i_%i",i,j),"dqdx",nbinq,lowq,highq);
      
      dqdxHist[i][j]->SetXTitle("dQ/dx [ADC/cm]");
      dqdxHist[i][j]->SetYTitle("number of hits");
      
      dqdxHist[i][j]->GetXaxis()->CenterTitle();
      dqdxHist[i][j]->GetYaxis()->CenterTitle();
    }
  }

}

#endif // HISTOGRAM_DEFINITIONS_H
