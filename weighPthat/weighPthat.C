using namespace std;
#include "weighPthat.h"

int weighPthat(TString ifname = "/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dstar_D0kpi/ntD_EvtBase_20160112_Dfinder_20151229_pp_Pythia8D0kpi_noweight.root",
               TString ofname = "/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dstar_D0kpi/ntD_EvtBase_20160112_Dfinder_20151229_pp_Pythia8D0kpi_withweight.root")
{
  Bool_t isInsidebin(Float_t xpthat, Float_t xmaxgenpt, Int_t i);
  cout<<endl;
  cout<<" -- Checking if input and output files are same"<<endl;
  if(ifname==ofname)
    {
      cout<<"    Error: Input file will be overwritten."<<endl;
      return 0;
    }
  cout<<" -- Opening unweighed sample"<<endl;
  TFile* inf = TFile::Open(ifname);
  TTree* ntGen = (TTree*)inf->Get("ntGen");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  Int_t Gsize; ntGen->SetBranchAddress("Gsize",&Gsize);
  Float_t Gpt[MAX_GEN]; ntGen->SetBranchAddress("Gpt",Gpt);
  Int_t GisSignal[MAX_GEN]; ntGen->SetBranchAddress("GisSignal",GisSignal);
  Float_t pthat; ntHi->SetBranchAddress("pthat",&pthat);

  Float_t weight[nBins],nweight[nBins];
  for(Int_t j=0;j<nBins;j++)
    {
      weight[j]=0;
      nweight[j]=0;
    }
  cout<<" -- Checking event number"<<endl;
  if(ntGen->GetEntries()!=ntHi->GetEntries())
    {
      cout<<"    Error: Gen tree and Hi tree have different event number."<<endl;
      return 0;
    }
  Int_t nentries = ntGen->GetEntries();
  cout<<" -- Calculating weights"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      ntGen->GetEntry(i);
      ntHi->GetEntry(i);
      if(i%100000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      Float_t maxpt=0;
      for(Int_t k=0;k<Gsize;k++)
        {
          if((GisSignal[k]==1||GisSignal[k]==2)&&Gpt[k]>maxpt) maxpt=Gpt[k];
        }
      for(Int_t j=0;j<nBins;j++)
        {
          if(isInsidebin(pthat,maxpt,j)) nweight[j]++;
        }
    }
  cout<<" -- Weight results"<<endl;
  for(Int_t j=0;j<nBins;j++)
    {
      if(nweight[j]==0)
        {
          cout<<"    Error: Weight fails."<<endl;
          return 0;
        }
      weight[j] = (crosssec[j]-crosssec[j+1])/nweight[j];
      cout<<"    Pthat"<<setiosflags(ios::left)<<setw(3)<<pthatBin[j]<<": "<<weight[j]<<endl;
    }

  cout<<" -- Building weight branch"<<endl;
  TFile* otf = TFile::Open(ofname,"update");
  TTree* ntHinew = (TTree*)otf->Get("ntHi");
  Float_t pthatweight,maxDgenpt;
  TBranch* newBr_pthatweight = ntHinew->Branch("pthatweight", &pthatweight, "pthatweight/F");
  TBranch* newBr_maxDgenpt = ntHinew->Branch("maxDgenpt", &maxDgenpt, "maxDgenpt/F");
  cout<<" -- Filling weight branch"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      ntGen->GetEntry(i);
      ntHi->GetEntry(i);
      if(i%100000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      pthatweight=0;
      Float_t maxpt=0;
      for(Int_t k=0;k<Gsize;k++)
        {
          if((GisSignal[k]==1||GisSignal[k]==2)&&Gpt[k]>maxpt) maxpt=Gpt[k];
        }
      for(Int_t j=0;j<nBins;j++)
        {
          maxDgenpt = maxpt;
          if(isInsidebin(pthat,maxpt,j)) pthatweight = weight[j];
        }
      newBr_pthatweight->Fill();
      newBr_maxDgenpt->Fill();
    }
  ntHinew->Write("", TObject::kOverwrite);

  cout<<" -- End"<<endl;
  cout<<endl;
  return 1;
}

Bool_t isInsidebin(Float_t xpthat, Float_t xmaxgenpt, Int_t i)
{
  if(i>=nBins)
    {
      cout<<"    Error: invalid input"<<endl;
      return false;
    }
  if(i<(nBins-1)&&xpthat>pthatBin[i]&&xpthat<pthatBin[i+1]&&xmaxgenpt>pthatBin[i]&&xmaxgenpt<pthatBin[i+1]) return true;
  else if(i==(nBins-1)&&xpthat>pthatBin[i]&&xmaxgenpt>pthatBin[i]) return true;
  else return false;
}

int main(int argc, char *argv[])
{
  if(argc==3)
    {
      weighPthat(argv[1], argv[2]);
      return 1;
    }
  else
    {
      std::cout<<"Invalid parameter number"<<std::endl;
      return 0;
    }
}
