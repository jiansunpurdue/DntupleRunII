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
  
  TCanvas*canvasRAA=new TCanvas("canvasRAA","canvasRAA",550,500);
  canvasRAA->cd();
  canvasRAA->SetLogx();
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,10.,120.,10.,0,1.5);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  hemptyEff->GetYaxis()->SetTitle("D^{0} R_{AA}");
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
  TLatex * tlatexeff=new TLatex(0.1612903,0.8525793,"CMS Preliminary     PbPb #sqrt{s}= 5.02 TeV");
  tlatexeff->SetNDC();
  tlatexeff->SetTextColor(1);
  tlatexeff->SetTextFont(42);
  tlatexeff->SetTextSize(0.04);
  tlatexeff->Draw();
  TLatex * tlatexeff2=new TLatex(0.1612903,0.7925793,"Centrality 0-100%");
  tlatexeff2->SetNDC();
  tlatexeff2->SetTextColor(1);
  tlatexeff2->SetTextFont(42);
  tlatexeff2->SetTextSize(0.04);
  tlatexeff2->Draw();
  TLatex * tlatexeff3=new TLatex(0.1612903,0.7325793,"D^{0} #rightarrow K#pi");
  tlatexeff3->SetNDC();
  tlatexeff3->SetTextColor(1);
  tlatexeff3->SetTextFont(42);
  tlatexeff3->SetTextSize(0.04);
  tlatexeff3->Draw();
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


