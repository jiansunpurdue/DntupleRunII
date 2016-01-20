#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"


void NuclearModification(TString inputPP="hPtSpectrumDzeroPP.root", TString inputPbPb="hPtSpectrumDzeroPbPb.root", TString inputprescalesPP="prescalePP.root", TString inputprescalesPbPb="prescalePbPb.root")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);
  
  TFile* filePP = new TFile(inputPP.Data());
  TH1F* hEffPP = (TH1F*)filePP->Get("hEff");
  TH1F* hSigmaPPStat = (TH1F*)filePP->Get("hPtSigma");
  TFile*fprescalesPP=new TFile(inputprescalesPP.Data());
  TFile*fprescalesPbPb=new TFile(inputprescalesPbPb.Data());


  TFile* filePbPb = new TFile(inputPbPb.Data());
  TH1F* hEffPbPb = (TH1F*)filePbPb->Get("hEff");
  TH1F* hRAA = (TH1F*)filePbPb->Get("hPtSigma");
  TH1F* hSigmaPbPbStat = (TH1F*)hRAA->Clone("hSigmaPbPbStat");
  hRAA->SetName("hRAA");
  hSigmaPbPbStat->SetName("hSigmaPbPbStat");
  hSigmaPbPbStat->Scale(1./(208.*208.));
  hRAA->Divide(hSigmaPPStat);
  hRAA->Scale(1./(208.*208.));
  std::cout<<hRAA->GetBinContent(1)<<std::endl;

  
  TH1F*hPrescalesPtBinsPP=(TH1F*)fprescalesPP->Get("hPrescalesPtBins");
  TH1F*hPrescalesPtBinsPbPb=(TH1F*)fprescalesPbPb->Get("hPrescalesPtBins");
  for (int i=0;i<nBins;i++) {
    hSigmaPPStat->SetBinContent(i+1,hSigmaPPStat->GetBinContent(i+1)/hPrescalesPtBinsPP->GetBinContent(i+1));
    hRAA->SetBinContent(i+1,hRAA->GetBinContent(i+1)/hPrescalesPtBinsPbPb->GetBinContent(i+1));
  }

  
  TCanvas* cRAA = new TCanvas("cRAA","",1000,500);
  cRAA->Divide(2,1);
  cRAA->cd(1);
  TH2F* hemptySigma=new TH2F("hemptySigma","",50,10.,110.,10.,0.11,10000.0);
  hemptySigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptySigma->GetYaxis()->CenterTitle();
  hemptySigma->GetYaxis()->SetTitle("R_{AA}");
  hemptySigma->GetXaxis()->SetTitleOffset(1.0);
  hemptySigma->GetYaxis()->SetTitleOffset(1.0);
  hemptySigma->GetXaxis()->SetTitleSize(0.04);
  hemptySigma->GetYaxis()->SetTitleSize(0.04);
  hemptySigma->GetXaxis()->SetTitleFont(42);
  hemptySigma->GetYaxis()->SetTitleFont(42);
  hemptySigma->GetXaxis()->SetLabelFont(42);
  hemptySigma->GetYaxis()->SetLabelFont(42);
  hemptySigma->GetXaxis()->SetLabelSize(0.03);
  hemptySigma->GetYaxis()->SetLabelSize(0.03);  
  hemptySigma->Draw();
  hSigmaPPStat->Draw("psame");
  cRAA->cd(2);
  TH2F* hemptyRAA=new TH2F("hemptyRAA","",50,10.,110.,10.,0.,20.0);
  hemptyRAA->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyRAA->GetYaxis()->CenterTitle();
  hemptyRAA->GetYaxis()->SetTitle("R_{AA}");
  hemptyRAA->GetXaxis()->SetTitleOffset(1.0);
  hemptyRAA->GetYaxis()->SetTitleOffset(1.0);
  hemptyRAA->GetXaxis()->SetTitleSize(0.04);
  hemptyRAA->GetYaxis()->SetTitleSize(0.04);
  hemptyRAA->GetXaxis()->SetTitleFont(42);
  hemptyRAA->GetYaxis()->SetTitleFont(42);
  hemptyRAA->GetXaxis()->SetLabelFont(42);
  hemptyRAA->GetYaxis()->SetLabelFont(42);
  hemptyRAA->GetXaxis()->SetLabelSize(0.03);
  hemptyRAA->GetYaxis()->SetLabelSize(0.03);  
  hemptyRAA->Draw();
  hRAA->Draw("psame");
  cRAA->SaveAs("cRAA.pdf");

  TCanvas* cEffPbPb = new TCanvas("cEffPbPb","",600,500);
  cEffPbPb->cd();
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,10.,110.,10.,0,1.);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  hemptyEff->GetYaxis()->SetTitle("efficiency");
  hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyEff->GetXaxis()->SetTitleOffset(1.);
  hemptyEff->GetYaxis()->SetTitleOffset(1.1);
  hemptyEff->GetXaxis()->SetTitleSize(0.045);
  hemptyEff->GetYaxis()->SetTitleSize(0.045);
  hemptyEff->GetXaxis()->SetTitleFont(42);
  hemptyEff->GetYaxis()->SetTitleFont(42);
  hemptyEff->GetXaxis()->SetLabelFont(42);
  hemptyEff->GetYaxis()->SetLabelFont(42);
  hemptyEff->GetXaxis()->SetLabelSize(0.04);
  hemptyEff->GetYaxis()->SetLabelSize(0.04);  
  hemptyEff->SetMaximum(2);
  hemptyEff->SetMinimum(0.);
  hemptyEff->Draw();
  hEffPbPb->Draw("same");
  cEffPbPb->SaveAs("efficiencyPbPb.pdf");  

}


int main(int argc, char *argv[])
{
  if((argc != 5))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc ==5)
    NuclearModification(argv[1], argv[2], argv[3], argv[4]);
  return 0;
}


