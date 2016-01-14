#include "uti.h"
#include "parameters.h"
Double_t minhisto=	1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;


void countsDzeromc(){

  TString inputmc="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi/ntD_EvtBase_20160112_Dfinder_20151229_pp_Pythia8_prompt_D0pt30p0_Pthat30_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20151212_dPt1tkPt1_D0Ds.root";
  TString cut="Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)&&(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5";
  TString selmcgen="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))";
  TString outputfile="mytest.root";

  TFile* infMC = new TFile(inputmc.Data());
  TTree* nt = (TTree*)infMC->Get("ntDkpi");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  
  TFile*fitfile = new TFile("testingWeight.root");
  TF1*myfit=(TF1*)fitfile->Get("myfit");

  nt->AddFriend(ntGen);

  float PVx;
  int Dsize;
  int Gsize;
  float Dpt[100000];
  float Gpt[100000];
  float Gy[100000];
  Int_t GisSignal[100000];
  float Dmass[100000];
  float Dtrk1Pt[100000];
  float Dtrk2Pt[100000];
  float Dchi2cl[100000];
  float Dalpha[100000];
  bool Dtrk1highPurity[100000];
  bool Dtrk2highPurity[100000];
  float Dy[100000];
  float Dgen[100000];
  float DsvpvDistance[100000];
  float DsvpvDisErr[100000];

  nt->SetBranchAddress("PVx",&PVx);
  nt->SetBranchAddress("Dsize",&Dsize);
  nt->SetBranchAddress("Gsize",&Gsize);
  nt->SetBranchAddress("Dpt",Dpt);  
  nt->SetBranchAddress("Dgen",Dgen);  
  nt->SetBranchAddress("Gpt",Gpt);  
  nt->SetBranchAddress("Gy",Gy);  
  nt->SetBranchAddress("Dmass",Dmass);  
  nt->SetBranchAddress("Dtrk1Pt",Dtrk1Pt);  
  nt->SetBranchAddress("Dtrk2Pt",Dtrk2Pt);  
  nt->SetBranchAddress("Dchi2cl",Dchi2cl);  
  nt->SetBranchAddress("Dalpha",Dalpha);  
  nt->SetBranchAddress("Dtrk1highPurity",Dtrk1highPurity);  
  nt->SetBranchAddress("Dtrk2highPurity",Dtrk2highPurity);  
  nt->SetBranchAddress("Dy",Dy);  
  nt->SetBranchAddress("DsvpvDistance",DsvpvDistance);  
  nt->SetBranchAddress("DsvpvDisErr",DsvpvDisErr);  
  nt->SetBranchAddress("GisSignal",GisSignal);  
  

  TH1D* hPtMCWeight = new TH1D("hPtMCWeight","",nBins,ptBins);
  TH1D* hPtGenWeight = new TH1D("hPtGenWeight","",nBins,ptBins);
  TH1D* hPtMC = new TH1D("hPtMC","",nBins,ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",nBins,ptBins);

  double par0=myfit->GetParameter(0);
  double par1=myfit->GetParameter(1);
  double par2=myfit->GetParameter(2);
  double par3=myfit->GetParameter(3);
  
  for(int i = 0; i<nt->GetEntries(); i++){
    nt->GetEntry(i);
    for (int j=0; j<Gsize; j++){     
      if((GisSignal[j]==1 || GisSignal[j]==2)&&(Gy[j]>-1&&Gy[j]<1.)){ 
        double weight=par0+par1*Gpt[j]+par2*Gpt[j]*Gpt[j]+par3*Gpt[j]*Gpt[j]*Gpt[j];
        cout<<"Gpt[j]="<<Gpt[j]<<", weight="<<weight<<endl;
        hPtGenWeight->Fill(Gpt[j],weight);
        hPtGen->Fill(Gpt[j]);
      }
    }
    for (int j=0; j<Dsize; j++){      
      if(Dtrk1highPurity[j]==1&& Dtrk2highPurity[j]==1 && Dy[j]>-1. && Dy[j]<1. && DsvpvDistance[j]/DsvpvDisErr[j]>3.5 && Dalpha[j]<0.12 && Dtrk1Pt[j]>1.5 && Dtrk2Pt[j]>1.5 && Dgen[j]==23333) {
        double weight=par0+par1*Dpt[j]+par2*Dpt[j]*Dpt[j]+par3*Dpt[j]*Dpt[j]*Dpt[j];
        hPtMCWeight->Fill(Dpt[j],weight);     
        hPtMC->Fill(Dpt[j]);
      }
    }
  }//end of global loop
    divideBinWidth(hPtGen);
    divideBinWidth(hPtMC);
    divideBinWidth(hPtGenWeight);
    divideBinWidth(hPtMCWeight);

    TH1D* heffWeight=(TH1D*)hPtMCWeight->Clone("heffWeight");
    heffWeight->Divide(hPtGenWeight);
    TH1D* heff=(TH1D*)hPtMC->Clone("heff");
    heff->Divide(hPtGen);

    TCanvas*canvas=new TCanvas("canvas","canvas",1000,500);
    canvas->Divide(3,2);
    canvas->cd(1);
    hPtGen->SetTitle(";generated pt no weight");
    hPtGen->Draw("p");
    canvas->cd(2);
    hPtMC->SetTitle(";reco pt no weight");
    hPtMC->Draw("p");
    canvas->cd(3);
    heff->SetTitle(";eff no weight");
    heff->SetMinimum(0);
    heff->SetMaximum(1);
    heff->Draw("p");
    canvas->cd(4);
    hPtGenWeight->SetTitle(";generated pt with weight");
    hPtGenWeight->Draw("p");
    canvas->cd(5);
    hPtMCWeight->SetTitle(";reco pt with weight");
    hPtMCWeight->Draw("p");
    canvas->cd(6);
    heffWeight->SetMinimum(0);
    heffWeight->SetMaximum(1);
    heffWeight->SetTitle(";eff with weight");
    heffWeight->Draw("p");

    TFile*fout=new TFile("test.root","recreate");

}