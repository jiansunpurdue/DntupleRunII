#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"


void NuclearModificationFactor(TString inputPP="CrossSectionFONLLPP.root", TString inputPbPb="CrossSectionFONLLPbPb.root")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);
  
  TFile *fpp=new TFile(inputPP.Data());
  TFile *fPbPb=new TFile(inputPbPb.Data());
  
  TH1D*hSigmaPPStat=(TH1D*)fpp->Get("hPtSigma");
  hSigmaPPStat->SetName("hSigmaPPStat");
  TH1D*hNuclearModification=(TH1D*)fPbPb->Get("hPtSigma");
  hNuclearModification->SetName("hNuclearModification");
  //hNuclearModification->SetBinContent(1,0);
  //hNuclearModification->SetBinError(1,0);

  hNuclearModification->Scale(1/208./208.);
  hNuclearModification->Divide(hSigmaPPStat);
  
  TCanvas*canvasRAA=new TCanvas("canvasRAA","canvasRAA",500,500);
  canvasRAA->cd();
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,0.,120.,10.,0,1.5);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  hemptyEff->GetYaxis()->SetTitle("Nuclear Modification");
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
  hNuclearModification->Draw("same");
  canvasRAA->SaveAs("canvasRAA.pdf");
}


int main(int argc, char *argv[])
{
  if((argc != 3))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc ==3)
    NuclearModificationFactor(argv[1], argv[2]);
  return 0;
}


