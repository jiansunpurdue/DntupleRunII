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

  hNuclearModification->Scale(1/208./208.);
  hNuclearModification->Divide(hSigmaPPStat);
  
  TCanvas*canvasRAA=new TCanvas("canvasRAA","canvasRAA",500,500);
  canvasRAA->cd();
  hNuclearModification->SetMinimum(0.);
  hNuclearModification->SetMaximum(1.5);
  hNuclearModification->GetYaxis()->SetTitle("Nuclear Modification");
  hNuclearModification->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hNuclearModification->Draw();
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


