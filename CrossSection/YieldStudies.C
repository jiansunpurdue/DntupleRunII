#include "uti.h"
#include "TLegendEntry.h"


void YieldStudies(TString inputPP="hPtSpectrumDzeroPP.root",TString inputPPMC="hPtSpectrumDzeroPPMCClosure.root")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);


  TFile* filesPP=new TFile(inputPP.Data());
  TFile* filesPPMC=new TFile(inputPPMC.Data());
  TH1D* hMeanPP=(TH1D*)filesPP->Get("hMean");
  TH1D* hMeanPPMC=(TH1D*)filesPPMC->Get("hMean");
  TH1D* hSigmaGaus1PP=(TH1D*)filesPP->Get("hSigmaGaus1");
  TH1D* hSigmaGaus1PPMC=(TH1D*)filesPPMC->Get("hSigmaGaus1");
  TH1D* hSigmaGaus2PP=(TH1D*)filesPP->Get("hSigmaGaus2");
  TH1D* hSigmaGaus2PPMC=(TH1D*)filesPPMC->Get("hSigmaGaus2");
  TH1D* hRelMagnGaus1Gaus2PP=(TH1D*)filesPP->Get("hRelMagnGaus1Gaus2");
  TH1D* hRelMagnGaus1Gaus2PPMC=(TH1D*)filesPPMC->Get("hRelMagnGaus1Gaus2");
  
  TH1D* hClosure=(TH1D*)filesPPMC->Get("hPtCor");
  TH1D* hPtGen=(TH1D*)filesPPMC->Get("hPtGen");
  hClosure->Sumw2();
  hClosure->Divide(hPtGen);


  TH2F* hemptySigma1=new TH2F("hemptySigma1","",50,0.,100.,10.,0.,0.1);  
  hemptySigma1->GetXaxis()->CenterTitle();
  hemptySigma1->GetYaxis()->CenterTitle();
  hemptySigma1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptySigma1->GetYaxis()->SetTitle("#sigma_{1} (GeV/c^{2})");
  hemptySigma1->GetXaxis()->SetTitleOffset(1.);
  hemptySigma1->GetYaxis()->SetTitleOffset(1.4);
  hemptySigma1->GetXaxis()->SetTitleSize(0.045);
  hemptySigma1->GetYaxis()->SetTitleSize(0.05);
  hemptySigma1->GetXaxis()->SetTitleFont(42);
  hemptySigma1->GetYaxis()->SetTitleFont(42);
  hemptySigma1->GetXaxis()->SetLabelFont(42);
  hemptySigma1->GetYaxis()->SetLabelFont(42);
  hemptySigma1->GetXaxis()->SetLabelSize(0.04);
  hemptySigma1->GetYaxis()->SetLabelSize(0.04);  
 
  TH2F* hemptyRelMag=new TH2F("hemptyRelMag","",50,0.,100.,10.,0.,0.3);  
  hemptyRelMag->GetXaxis()->CenterTitle();
  hemptyRelMag->GetYaxis()->CenterTitle();
  hemptyRelMag->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyRelMag->GetYaxis()->SetTitle("norm(#sigma_1)/norm(#sigma_1)+norm(#sigma_2)");
  hemptyRelMag->GetXaxis()->SetTitleOffset(1.);
  hemptyRelMag->GetYaxis()->SetTitleOffset(1.4);
  hemptyRelMag->GetXaxis()->SetTitleSize(0.045);
  hemptyRelMag->GetYaxis()->SetTitleSize(0.05);
  hemptyRelMag->GetXaxis()->SetTitleFont(42);
  hemptyRelMag->GetYaxis()->SetTitleFont(42);
  hemptyRelMag->GetXaxis()->SetLabelFont(42);
  hemptyRelMag->GetYaxis()->SetLabelFont(42);
  hemptyRelMag->GetXaxis()->SetLabelSize(0.04);
  hemptyRelMag->GetYaxis()->SetLabelSize(0.04);  

  TH2F* hemptySigma2=new TH2F("hemptySigma2","",50,0.,100.,10.,0.,0.03);  
  hemptySigma2->GetXaxis()->CenterTitle();
  hemptySigma2->GetYaxis()->CenterTitle();
  hemptySigma2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptySigma2->GetYaxis()->SetTitle("#sigma_{2} (GeV/c^{2})");
  hemptySigma2->GetXaxis()->SetTitleOffset(1.);
  hemptySigma2->GetYaxis()->SetTitleOffset(1.4);
  hemptySigma2->GetXaxis()->SetTitleSize(0.045);
  hemptySigma2->GetYaxis()->SetTitleSize(0.05);
  hemptySigma2->GetXaxis()->SetTitleFont(42);
  hemptySigma2->GetYaxis()->SetTitleFont(42);
  hemptySigma2->GetXaxis()->SetLabelFont(42);
  hemptySigma2->GetYaxis()->SetLabelFont(42);
  hemptySigma2->GetXaxis()->SetLabelSize(0.04);
  hemptySigma2->GetYaxis()->SetLabelSize(0.04);  
  
  TH2F* hemptyMean=new TH2F("hemptyMean","",50,0.,110.,10.,1.855,1.88);  
  hemptyMean->GetXaxis()->CenterTitle();
  hemptyMean->GetYaxis()->CenterTitle();
  hemptyMean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyMean->GetYaxis()->SetTitle("Mean (GeV/c^{2})");
  hemptyMean->GetXaxis()->SetTitleOffset(1.);
  hemptyMean->GetYaxis()->SetTitleOffset(1.4);
  hemptyMean->GetXaxis()->SetTitleSize(0.045);
  hemptyMean->GetYaxis()->SetTitleSize(0.05);
  hemptyMean->GetXaxis()->SetTitleFont(42);
  hemptyMean->GetYaxis()->SetTitleFont(42);
  hemptyMean->GetXaxis()->SetLabelFont(42);
  hemptyMean->GetYaxis()->SetLabelFont(42);
  hemptyMean->GetXaxis()->SetLabelSize(0.04);
  hemptyMean->GetYaxis()->SetLabelSize(0.04);  

  TCanvas* cYields = new TCanvas("cYields","",1000,1000);
  cYields->Divide(2,2);
  cYields->cd(1);
  hemptyMean->Draw();
  hMeanPP->SetLineWidth(2);
  hMeanPP->Draw("same");
  hMeanPPMC->SetLineColor(2);
  hMeanPPMC->SetLineWidth(2);
  hMeanPPMC->Draw("same");
  TLegend *legendMean=new TLegend(0.4958166,0.7558707,0.7949297,0.9299148,"");
  legendMean->SetBorderSize(0);
  legendMean->SetLineColor(0);
  legendMean->SetFillColor(0);
  legendMean->SetFillStyle(1001);
  legendMean->SetTextFont(42);
  legendMean->SetTextSize(0.045);
  TLegendEntry *ent_MeanPP=legendMean->AddEntry(hMeanPP,"Data","pf");
  ent_MeanPP->SetTextFont(42);
  ent_MeanPP->SetLineColor(1);
  TLegendEntry *ent_MeanPPMC=legendMean->AddEntry(hMeanPPMC,"MC","pf");
  ent_MeanPPMC->SetTextFont(42);
  ent_MeanPPMC->SetLineColor(2);
  legendMean->Draw("same");
  TLine* l = new TLine(5,1.865,100,1.865);
  l->SetLineWidth(2);
  l->SetLineStyle(2);
  l->Draw();
  cYields->cd(2);
  hemptySigma1->Draw();
  hSigmaGaus1PP->SetLineWidth(2);
  hSigmaGaus1PP->Draw("same");
  hSigmaGaus1PPMC->SetLineColor(2);
  hSigmaGaus1PPMC->SetLineWidth(2);
  hSigmaGaus1PPMC->Draw("same");
  TLegend *legendSigma1=new TLegend(0.4958166,0.7558707,0.7949297,0.9299148,"");
  legendSigma1->SetBorderSize(0);
  legendSigma1->SetLineColor(0);
  legendSigma1->SetFillColor(0);
  legendSigma1->SetFillStyle(1001);
  legendSigma1->SetTextFont(42);
  legendSigma1->SetTextSize(0.045);
  TLegendEntry *ent_Sigma1PP=legendSigma1->AddEntry(hSigmaGaus1PP,"Data","pf");
  ent_Sigma1PP->SetTextFont(42);
  ent_Sigma1PP->SetLineColor(1);
  TLegendEntry *ent_Sigma1PPMC=legendSigma1->AddEntry(hSigmaGaus1PPMC,"MC","pf");
  ent_Sigma1PPMC->SetTextFont(42);
  ent_Sigma1PPMC->SetLineColor(2);
  legendSigma1->Draw("same");

  cYields->cd(3);
  hemptySigma2->Draw();
  hSigmaGaus2PP->SetLineWidth(2);
  hSigmaGaus2PP->Draw("same");
  hSigmaGaus2PPMC->SetLineColor(2);
  hSigmaGaus2PPMC->SetLineWidth(2);
  hSigmaGaus2PPMC->Draw("same");
  TLegend *legendSigma2=new TLegend(0.4958166,0.7558707,0.7949297,0.9299148,"");
  legendSigma2->SetBorderSize(0);
  legendSigma2->SetLineColor(0);
  legendSigma2->SetFillColor(0);
  legendSigma2->SetFillStyle(1001);
  legendSigma2->SetTextFont(42);
  legendSigma2->SetTextSize(0.045);
  TLegendEntry *ent_Sigma2PP=legendSigma2->AddEntry(hSigmaGaus2PP,"Data","pf");
  ent_Sigma2PP->SetTextFont(42);
  ent_Sigma2PP->SetLineColor(1);
  TLegendEntry *ent_Sigma2PPMC=legendSigma2->AddEntry(hSigmaGaus2PPMC,"MC","pf");
  ent_Sigma2PPMC->SetTextFont(42);
  ent_Sigma2PPMC->SetLineColor(2);
  legendSigma2->Draw("same");
  cYields->cd(4);
  hemptyRelMag->Draw();
  hRelMagnGaus1Gaus2PP->SetLineWidth(2);
  hRelMagnGaus1Gaus2PP->Draw("same");
  
  TH2F* hemptyClosure=new TH2F("hemptyClosure","",50,0.,100.,10.,0.7,1.3);  
  hemptyClosure->GetXaxis()->CenterTitle();
  hemptyClosure->GetYaxis()->CenterTitle();
  hemptyClosure->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyClosure->GetYaxis()->SetTitle("Corr. yield/Gen");
  hemptyClosure->GetXaxis()->SetTitleOffset(1.);
  hemptyClosure->GetYaxis()->SetTitleOffset(1.4);
  hemptyClosure->GetXaxis()->SetTitleSize(0.045);
  hemptyClosure->GetYaxis()->SetTitleSize(0.05);
  hemptyClosure->GetXaxis()->SetTitleFont(42);
  hemptyClosure->GetYaxis()->SetTitleFont(42);
  hemptyClosure->GetXaxis()->SetLabelFont(42);
  hemptyClosure->GetYaxis()->SetLabelFont(42);
  hemptyClosure->GetXaxis()->SetLabelSize(0.04);
  hemptyClosure->GetYaxis()->SetLabelSize(0.04);  
  cYields->SaveAs("cYieldStudies.pdf");
  
  TCanvas* cClosure = new TCanvas("cClosure","",500,500);
  cClosure->cd(1);
  hemptyClosure->Draw();
  hClosure->SetLineWidth(2);
  hClosure->Draw("same");
  TLine* lunity = new TLine(5,1,100,1);
  lunity->SetLineWidth(2);
  lunity->SetLineStyle(2);
  lunity->Draw("same");
  cClosure->SaveAs("cClosure.pdf");

  }


int main(int argc, char *argv[])
{
  if((argc != 3))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 3)
    YieldStudies(argv[1], argv[2]);
  return 0;
}


