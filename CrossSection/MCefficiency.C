#include "uti.h"
#include "parameters.h"

Double_t setparam0=100.;
Double_t setparam1=1.865;
Double_t setparam2=0.03;
Double_t setparam10=0.005;
Double_t setparam8=0.1;
Double_t setparam9=0.1;
Double_t fixparam1=1.865;
Double_t minhisto=	1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

TString weight = "1";
TString selmc;

void MCefficiency(TString inputmc, TString selmcgen,TString selmcgenacceptance, TString cut_recoonly, TString cut,TString label)
{
  selmc = Form("%s",cut.Data());

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);


  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("ntDkpi");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  ntMC->AddFriend(ntGen);
  

  TTree* nthi = (TTree*)infMC->Get("ntHi");
  ntGen->AddFriend(nthi);
  nthi->AddFriend(ntMC);
  ntMC->AddFriend(nthi);
  
  TH1D* hPtMC = new TH1D("hPtMC","",nBins,ptBins);
  TH1D* hPtMCrecoonly = new TH1D("hPtMCrecoonly","",nBins,ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",nBins,ptBins);
  TH1D* hPtGenAcc = new TH1D("hPtGenAcc","",nBins,ptBins);
  TH1D* hpthat = new TH1D("hpthat","",100,0,100);
  TH1D* hpthatweight = new TH1D("hpthatweight","",100,0,100);

  ntMC->Project("hPtMC","Dpt",TCut(weight)*(TCut(selmc.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMC);
  ntMC->Project("hPtMCrecoonly","Dpt",TCut(weight)*(TCut(cut_recoonly.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMCrecoonly);
  ntGen->Project("hPtGen","Gpt",TCut(weight)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGen);
  ntGen->Project("hPtGenAcc","Gpt",TCut(weight)*(TCut(selmcgenacceptance.Data())));
  divideBinWidth(hPtGenAcc);

  hPtMC->Sumw2();
  TH1D* hEff = (TH1D*)hPtMC->Clone("hEff");
  hEff->Divide(hPtGen);

  TH1D* hEffReco = (TH1D*)hPtMCrecoonly->Clone("hEffReco");
  hEffReco->Sumw2();
  hEffReco->Divide(hPtGen);

  TH1D* hEffAcc = (TH1D*)hPtGenAcc->Clone("hEffAcc");
  hEffAcc->Sumw2();
  hEffAcc->Divide(hEffAcc,hPtGen,1,1,"b");
  
  TH1D* hEffSelection = (TH1D*)hPtMC->Clone("hEffSelection");
  hEffSelection->Sumw2();
  hEffSelection->Divide(hEffSelection,hPtMCrecoonly,1,1,"b");
  
  
  TCanvas*canvasEff=new TCanvas("canvasEff","canvasEff",1000.,1000);
  canvasEff->Divide(2,2);
  canvasEff->cd(1);
  
  hEffAcc->SetXTitle("Gen D p_{T}");
  hEffAcc->SetYTitle("#alpha");
  hEffAcc->SetMinimum(0);
  hEffAcc->SetMaximum(1.5);
  hEffAcc->GetYaxis()->SetTitleOffset(1.2);
  hEffAcc->Draw();
  canvasEff->cd(2);
  hEffReco->SetMinimum(0);
  hEffReco->SetMaximum(1.5);
  hEffReco->SetXTitle("D p_{T}");
  hEffReco->SetYTitle("#alpha x #epsilon_{reco}");
  hEffReco->GetYaxis()->SetTitleOffset(1.2);
  hEffReco->Draw();
  canvasEff->cd(3);
  hEffSelection->SetMinimum(0);
  hEffSelection->SetMaximum(1.5);
  hEffSelection->SetXTitle("D p_{T}");
  hEffSelection->SetYTitle("#epsilon_{sel}");
  hEffSelection->GetYaxis()->SetTitleOffset(1.2);
  hEffSelection->Draw();  
  canvasEff->cd(4);
  hEff->SetMinimum(0);
  hEff->SetMaximum(1.5);
  hEff->Sumw2();
  hEff->SetXTitle("D p_{T}");
  hEff->SetYTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hEff->GetYaxis()->SetTitleOffset(1.2);
  hEff->Draw();
  canvasEff->SaveAs(Form("canvasEff_study%s.pdf",Form(label.Data())));

}

int main(int argc, char *argv[])
{
  if((argc != 7))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 7)
    MCefficiency(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
  return 0;
}

