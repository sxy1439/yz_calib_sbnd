#include "root_includes.h"
#include "landaugausfit.h"
#include "function.h"
#include "hist_def.h"


using std::cout;
using std::endl;
using std::ifstream;
using std::string;

//const float e_callife[2]={10.e3, 10.e3};// in us

const float ticksToMs = 0.0005; // 500 ns per tick
const float nsToMs = 1e-6; // 1 ns to 0.000001 ms
const float preTriggerWindow = 0.2; /* my pre-trigger window in ms */;

const float lowfq=800, highfq=1400;



bool Is_Edge(float x, float y, float z);
bool Is_Cathode_Crossing(float startx, float endx);
float getMedian(std::vector<double>& values);
float getMean(std::vector<double>& values);
float drift_time(float ti, float t0, float ticksToMs, float preTriggerWindow, float nsToMs);
float lifetime_correction(float ti, float t0, float ticksToMs, float preTriggerWindow, float nsToMs, float tau);

double langaufun(double *x, double *par);
TF1 *langaufit(TH1F *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF);
