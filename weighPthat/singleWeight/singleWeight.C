using namespace std;
#include "uti.h"

int singleWeight(TString ifname = "/afs/cern.ch/work/w/wangj/public/RunII/weighPthat/test.root",
                 Float_t weight = 1.)
{
  cout<<endl;
  cout<<" -- Processing file"<<endl;
  cout<<"    "<<ifname<<endl;
  cout<<" -- Weight is"<<endl;
  cout<<"    "<<weight<<endl;
  cout<<" -- Opening unweighed sample"<<endl;
  TFile* inf = TFile::Open(ifname,"update");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  Int_t nentries = ntHi->GetEntries();
  cout<<" -- Building weight branch"<<endl;
  Float_t pthatweight;
  TBranch* newBr_pthatweight = ntHi->Branch("pthatweight", &pthatweight, "pthatweight/F");
  cout<<" -- Filling weight branch"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      ntHi->GetEntry(i);
      if(i%10000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      pthatweight = weight;
      newBr_pthatweight->Fill();
    }
  ntHi->Write("", TObject::kOverwrite);
  cout<<" -- End"<<endl;
  return 1;
}

int main(int argc, char *argv[])
{
  float fweight;
  if(argc==3)
    {
      fweight = atof(argv[2]);
      singleWeight(argv[1], fweight);
      return 1;
    }
  else
    {
      std::cout<<"Invalid parameter number"<<std::endl;
      return 0;
    }
}
