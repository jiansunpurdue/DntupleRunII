const int nBins=6 ;
double ptBins[nBins+1] = {20.,25., 30.,35.,40.,50.,60.};
int triggerassignment[nBins] = {0,0,1,1,1,2};

const int ntriggerspp=3;
TString triggerHLTpp[ntriggerspp]={"HLT_DmesonPPTrackingGlobal_Dpt15_v1","HLT_DmesonPPTrackingGlobal_Dpt30_v1","HLT_DmesonPPTrackingGlobal_Dpt50_v1"};
TString triggerL1pp[ntriggerspp]={"L1_SingleJet24_BptxAND","L1_SingleJet40_BptxAND","L1_SingleJet48_BptxAND"};

//nothing from here on needs to be changed
const int BIN_NUM= 792;
const int HMIN=2;
const int HMAX=200;
const double binsize=((double)(HMAX)-(double)(HMIN))/(double)(BIN_NUM);
Double_t BRchain=0.0388;
