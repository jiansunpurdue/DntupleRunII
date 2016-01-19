/*
const int nBins=6;
double ptBins[nBins+1] = {30.,35.,40.,45.,50.,60.,100.};
const int ntriggers=3;
TString triggerHLTPP[ntriggers]={"HLT_DmesonPPTrackingGlobal_Dpt15_v1","HLT_DmesonPPTrackingGlobal_Dpt30_v1","HLT_DmesonPPTrackingGlobal_Dpt50_v1"};
int triggerassignmentPP[nBins]= {0,0,1,1,1,2};
TString triggerHLTPbPb[ntriggers]={"HLT_HIDmesonHITrackingGlobal_Dpt20_v1","HLT_HIDmesonHITrackingGlobal_Dpt40_v1","HLT_HIDmesonHITrackingGlobal_Dpt60_v1"};
int triggerassignmentPbPb[nBins]= {0,0,0,1,1,2};
*/
/*
const int nBins=4;
double ptBins[nBins+1] = {15.,20.,30.,40.,60.};
*/

const int nBins=1;
double ptBins[nBins+1] = {60.,100.};
const int ntriggers=3;
TString triggerHLTPP[ntriggers]={"HLT_DmesonPPTrackingGlobal_Dpt15_v1","HLT_DmesonPPTrackingGlobal_Dpt30_v1","HLT_DmesonPPTrackingGlobal_Dpt50_v1"};
int triggerassignmentPP[nBins]= {2};
TString triggerHLTPbPb[ntriggers]={"HLT_HIDmesonHITrackingGlobal_Dpt20_v1","HLT_HIDmesonHITrackingGlobal_Dpt40_v1","HLT_HIDmesonHITrackingGlobal_Dpt60_v1"};
int triggerassignmentPbPb[nBins]= {2};

/*
//ALICE binning
const int nBins=9;
double ptBins[nBins+1] = {1.,2.,3.,4.,5.,6.,7.,8.,12.,16.};
const int ntriggers=1;
TString triggerHLTPP[ntriggers]={"HLT_DmesonPPTrackingGlobal_Dpt15_v1"};
int triggerassignmentPP[nBins]= {0,0,0,0,0,0,0,0,0};
TString triggerHLTPbPb[ntriggers]={"HLT_HIDmesonHITrackingGlobal_Dpt20_v1"};
int triggerassignmentPbPb[nBins]= {0,0,0,0,0,0,0,0,0};
*/


/*
//ATLAS binning for pp 7 TeV for D*
const int nBins=9;
double ptBins[nBins+1] = {3.5,5.0,6.5,8.0,12.0,20.0,30.,40.,60.,100.};
*/
//nothing from here on needs to be changed
const int BIN_NUM= 1196;
const int HMIN=1;
const int HMAX=300;

const double binsize=((double)(HMAX)-(double)(HMIN))/(double)(BIN_NUM);
Double_t BRchain=0.0388;
