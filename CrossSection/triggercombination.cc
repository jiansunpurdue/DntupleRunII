#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"
#include "TNtuple.h"
#include "TH1F.h"

void triggercombination(TString inputdata="/data/HeavyFlavourRun2/DfinderData_pp_20151218_dPt0tkPt1_D0Dstar3p5p/merged_ntuple.root",TString output="output.root"){

  TFile* inf = new TFile(inputdata.Data());
  TTree* nt = (TTree*) inf->Get("ntDkpi");
  TTree* HltTree = (TTree*) inf->Get("ntHlt");
  HltTree->AddFriend(nt);
  nt->AddFriend(HltTree);

  TH1D* hPrescalesPtBins = new TH1D("hPrescalesPtBins","",nBins,0,nBins);

  int triggervariable[ntriggerspp];                    
  double ntriggerscounters[ntriggerspp];       
  double prescale[ntriggerspp];
  double nflag[ntriggerspp];               
  double ncounters[ntriggerspp];
  double ncountersANDunprescaled[ntriggerspp];

  for (int index=0; index<ntriggerspp;index++){  
    triggervariable[index]=0;
    ntriggerscounters[index]=0.;
    prescale[index]=0.;
    nflag[index]=0.;
    ncounters[index]=0.;
    ncountersANDunprescaled[index]=0.;
  }

  for (int index=0; index<ntriggerspp;index++) nt->SetBranchAddress(triggerHLTpp[index].Data(),&triggervariable[index]);
  
  //int nevents_total = nt->GetEntries();    
  int nevents_total = 100000;
  for(int entry=0; entry<nevents_total; entry++){
    if((entry%10000)==0) printf("Loading event #%d of %d.\n",entry,nevents_total);
    nt->GetEntry(entry);
    
    for (int index=0; index<ntriggerspp;index++) {
      ncountersANDunprescaled[index]=ncountersANDunprescaled[index]+(triggervariable[index]&&triggervariable[2]);
      ncounters[index]=ncounters[index]+triggervariable[index];
    }
  }//end of for over events
    for (int index=0; index<ntriggerspp;index++){
     prescale[index]=ncountersANDunprescaled[index]/ncounters[2];
     cout<<"ncounters for trigger seed="<<prescale[index]<<endl; 
    }
    
    for (int index=0; index<nBins;index++){
     hPrescalesPtBins->SetBinContent(index+1,prescale[triggerassignment[index]]);
    }

    
    
    TFile*foutput=new TFile(output.Data(),"recreate");
    foutput->cd();
    hPrescalesPtBins->Write();
    delete foutput;
}


int main(int argc, char *argv[])
{
  if((argc != 3))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }

  if(argc == 3)
    triggercombination(argv[1],argv[2]);
  return 0;
}


