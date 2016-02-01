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


  double apt[nBins];
  //bin half width
  double aptl[nBins];
  //number of every rebined bin
  double bin_num[nBins];


  for (int ibin=0; ibin<nBins; ibin++){
    apt[ibin]=(ptBins[ibin+1]+ptBins[ibin])/2.;
    aptl[ibin] = (ptBins[ibin+1]-ptBins[ibin])/2;
    bin_num[ibin]=aptl[ibin]/binsize*2;
  }

  Double_t xr[nBins], yr[nBins], xrlow[nBins], yrlow[nBins],xrhigh[nBins],yrhigh[nBins];
  for(int i=0;i<nBins;i++)
  {
    yr[i] = hNuclearModification->GetBinContent(i+1);
    yrlow[i] = hNuclearModification->GetBinContent(i+1)*0.2;
    yrhigh[i] =hNuclearModification->GetBinContent(i+1)*0.2;
  }


  TGraphAsymmErrors* gNuclearModification = new TGraphAsymmErrors(nBins,apt,yr,aptl,aptl,yrlow,yrhigh);
  gNuclearModification->SetName("gNuclearModification");
  gNuclearModification->SetMarkerStyle(20);
  gNuclearModification->SetMarkerSize(0.8);

  TCanvas*canvasRAA=new TCanvas("canvasRAA","canvasRAA",550,500);
  canvasRAA->cd();
  canvasRAA->SetLogx();
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,10.,120.,10.,0,1.5);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  hemptyEff->GetYaxis()->SetTitle("D^{0} R_{AA}");
  hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyEff->GetXaxis()->SetTitleOffset(0.9);
  hemptyEff->GetYaxis()->SetTitleOffset(1.);
  hemptyEff->GetXaxis()->SetTitleSize(0.045);
  hemptyEff->GetYaxis()->SetTitleSize(0.045);
  hemptyEff->GetXaxis()->SetTitleFont(42);
  hemptyEff->GetYaxis()->SetTitleFont(42);
  hemptyEff->GetXaxis()->SetLabelFont(42);
  hemptyEff->GetYaxis()->SetLabelFont(42);
  hemptyEff->GetXaxis()->SetLabelSize(0.035);
  hemptyEff->GetYaxis()->SetLabelSize(0.035);  
  hemptyEff->SetMaximum(2);
  hemptyEff->SetMinimum(0.);
  hemptyEff->Draw();
  gNuclearModification->SetFillColor(4);
  gNuclearModification->SetFillStyle(3001); 
  gNuclearModification->SetLineWidth(1);
  gNuclearModification->SetLineColor(4);
  gNuclearModification->Draw("5same");
  hNuclearModification->SetLineWidth(3);
  hNuclearModification->Draw("same");

  TLatex * tlatexeff=new TLatex(0.1612903,0.8525793,"CMS Preliminary     PbPb #sqrt{s}= 5.02 TeV");
  tlatexeff->SetNDC();
  tlatexeff->SetTextColor(1);
  tlatexeff->SetTextFont(42);
  tlatexeff->SetTextSize(0.038);
  tlatexeff->Draw();
  TLatex * tlatexeff2=new TLatex(0.1612903,0.7925793,"Centrality 0-100%");
  tlatexeff2->SetNDC();
  tlatexeff2->SetTextColor(1);
  tlatexeff2->SetTextFont(42);
  tlatexeff2->SetTextSize(0.038);
  tlatexeff2->Draw();
  TLatex * tlatexeff3=new TLatex(0.1612903,0.7325793,"L^{pp}_{int} = 26.31pb^{-1}, L^{PbPb}_{int} = 404.4 #mub^{-1}");
  tlatexeff3->SetNDC();
  tlatexeff3->SetTextColor(1);
  tlatexeff3->SetTextFont(42);
  tlatexeff3->SetTextSize(0.038);
  tlatexeff3->Draw();

  TLegend *legendSigma=new TLegend(0.5100806,0.5168644,0.8084677,0.6605932,"");
  legendSigma->SetBorderSize(0);
  legendSigma->SetLineColor(0);
  legendSigma->SetFillColor(0);
  legendSigma->SetFillStyle(1001);
  legendSigma->SetTextFont(42);
  legendSigma->SetTextSize(0.045);

  TLegendEntry *ent_SigmaPP=legendSigma->AddEntry(hNuclearModification,"R_{AA} stat. unc.","pf");
  ent_SigmaPP->SetTextFont(42);
  ent_SigmaPP->SetLineColor(1);
  ent_SigmaPP->SetMarkerColor(1);
  ent_SigmaPP->SetTextSize(0.03);


  TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(gNuclearModification,"R_{AA} syst. (arbitrary 20%)","f");
  ent_Sigmapp->SetTextFont(42);
  ent_Sigmapp->SetLineColor(5);
  ent_Sigmapp->SetMarkerColor(1);
  ent_Sigmapp->SetTextSize(0.03);

  legendSigma->Draw("same");

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


