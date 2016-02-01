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
  
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,10.,120.,10.,0,1.5);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  hemptyEff->GetYaxis()->SetTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyEff->GetXaxis()->SetTitleOffset(0.9);
  hemptyEff->GetYaxis()->SetTitleOffset(0.95);
  hemptyEff->GetXaxis()->SetTitleSize(0.05);
  hemptyEff->GetYaxis()->SetTitleSize(0.05);
  hemptyEff->GetXaxis()->SetTitleFont(42);
  hemptyEff->GetYaxis()->SetTitleFont(42);
  hemptyEff->GetXaxis()->SetLabelFont(42);
  hemptyEff->GetYaxis()->SetLabelFont(42);
  hemptyEff->GetXaxis()->SetLabelSize(0.035);
  hemptyEff->GetYaxis()->SetLabelSize(0.035);  
  hemptyEff->SetMaximum(2);
  hemptyEff->SetMinimum(0.);
  hemptyEff->Draw();

  TH2F* hemptyEffAcc=(TH2F*)hemptyEff->Clone("hemptyEffAcc");
  TH2F* hemptyEffReco=(TH2F*)hemptyEff->Clone("hemptyEffReco");
  TH2F* hemptyEffSelection=(TH2F*)hemptyEff->Clone("hemptyEffSelection");
  
  TCanvas*canvasEff=new TCanvas("canvasEff","canvasEff",1100.,1000);
  canvasEff->Divide(2,2);
  canvasEff->cd(1);
  
  hemptyEffAcc->SetYTitle("#alpha");
  hemptyEffAcc->Draw();
  hEffAcc->Draw("same");
  canvasEff->cd(2);
  hEffReco->GetYaxis()->SetTitleOffset(1.2);
  hemptyEffReco->SetYTitle("#alpha x #epsilon_{reco} ");
  hemptyEffReco->Draw();
  hEffReco->Draw("same");
  canvasEff->cd(3);
  hemptyEffSelection->SetYTitle(" #epsilon_{sel}");
  hemptyEffSelection->Draw();  
  hEffSelection->Draw("same");  
  canvasEff->cd(4);
  hemptyEff->Draw();
  hEff->Draw("same");
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

