const int nBins=6;
double ptBins[nBins+1] = {30.,35.,40.,45.,50.,60.,100.};
const int ntriggers=3;
TString triggerHLTPP[ntriggers]={"HLT_DmesonPPTrackingGlobal_Dpt15_v1","HLT_DmesonPPTrackingGlobal_Dpt30_v1","HLT_DmesonPPTrackingGlobal_Dpt50_v1"};
int triggerassignmentPP[nBins]= {0,0,1,1,1,2};
TString triggerHLTPbPb[ntriggers]={"HLT_HIDmesonHITrackingGlobal_Dpt20_v1","HLT_HIDmesonHITrackingGlobal_Dpt40_v1","HLT_HIDmesonHITrackingGlobal_Dpt60_v1"};
int triggerassignmentPbPb[nBins]= {0,0,0,1,1,2};

//nothing from here on needs to be changed
const int BIN_NUM= 1196;
const int HMIN=1;
const int HMAX=300;

const double binsize=((double)(HMAX)-(double)(HMIN))/(double)(BIN_NUM);
Double_t BRchain=0.0388;
