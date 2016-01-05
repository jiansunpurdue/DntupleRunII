#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"
#include "TNtuple.h"
#include "TH1F.h"
using namespace std;


void triggercombination(TString ispp="PP",TString inputdata="/data/HeavyFlavourRun2/DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p/merged_ntuple.root",TString output="output.root"){

  TFile* inf = new TFile(inputdata.Data());
  TTree* nt = (TTree*) inf->Get("ntDkpi");
  TTree* HltTree = (TTree*) inf->Get("ntHlt");
  HltTree->AddFriend(nt);
  nt->AddFriend(HltTree);

  TH1D* hPrescalesPtBins = new TH1D("hPrescalesPtBins","",nBins,0,nBins);
  
  TString triggerHLT[ntriggers];
   int triggerassignment[nBins];
  
  if(ispp=="PP"){
     for (int index=0; index<ntriggers;index++) triggerHLT[index]=triggerHLTPP[index];
     for (int index=0; index<nBins;index++) triggerassignment[index]=triggerassignmentPP[index];
  }
  if(ispp=="PbPb"){
     for (int index=0; index<ntriggers;index++) triggerHLT[index]=triggerHLTPbPb[index];
     for (int index=0; index<nBins;index++) triggerassignment[index]=triggerassignmentPbPb[index];
  }

  int triggervariable[ntriggers];                    
  double ntriggerscounters[ntriggers];       
  double prescale[ntriggers];
  double nflag[ntriggers];               
  double ncounters[ntriggers];
  double ncountersANDunprescaled[ntriggers];

  for (int index=0; index<ntriggers;index++){  
    triggervariable[index]=0;
    ntriggerscounters[index]=0.;
    prescale[index]=0.;
    nflag[index]=0.;
    ncounters[index]=0.;
    ncountersANDunprescaled[index]=0.;
  }

  for (int index=0; index<ntriggers;index++) nt->SetBranchAddress(triggerHLT[index].Data(),&triggervariable[index]);
  //int nevents_total = nt->GetEntries();    
  int nevents_total = 200000;
  for(int entry=0; entry<nevents_total; entry++){
    if((entry%10000)==0) printf("Loading event #%d of %d.\n",entry,nevents_total);
    nt->GetEntry(entry);
    
    for (int index=0; index<ntriggers;index++) {
      ncountersANDunprescaled[index]=ncountersANDunprescaled[index]+(triggervariable[index]&&(triggervariable[2]==1));
      ncounters[index]=ncounters[index]+triggervariable[index];
    }
  }//end of for over events
    for (int index=0; index<ntriggers;index++){
     prescale[index]=ncountersANDunprescaled[index]/ncounters[2];
     cout<<"------index------"<<endl;
     cout<<"ncounters="<<ncounters[index]<<endl; 
     cout<<"ncountersANDunprescaled="<<ncountersANDunprescaled[index]<<endl; 
     cout<<"prescale="<<prescale[index]<<endl; 
    }
    
    for (int index=0; index<nBins;index++){
     hPrescalesPtBins->SetBinContent(index+1,prescale[triggerassignment[index]]);
    }
    hPrescalesPtBins->Draw();
    TFile*foutput=new TFile(output.Data(),"recreate");
    foutput->cd();
    hPrescalesPtBins->Write();
}


int main(int argc, char *argv[])
{
  if((argc != 4))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }

  if(argc == 4)
    triggercombination(argv[1],argv[2],argv[3]);
  return 0;
}


