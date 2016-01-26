#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"


void CrossSectionRatio(TString inputFONLL="output_pp_d0meson_5TeV_y1.root", TString inputPP="hPtSpectrumDzeroPP.root", TString inputprescalesPP="prescalePP.root",int usePrescaleCorr=1)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);

  TFile* filePPReference = new TFile(inputFONLL.Data());  
  TGraphAsymmErrors* gaeBplusReference = (TGraphAsymmErrors*)filePPReference->Get("gaeSigmaDzero");
  
  TFile* filePP = new TFile(inputPP.Data());
  TH1F* hEffPP = (TH1F*)filePP->Get("hEff");
  TH1F* hSigmaPPStat = (TH1F*)filePP->Get("hPtSigma");
  TFile*fprescalesPP=new TFile(inputprescalesPP.Data());
  
  TH1F*hPrescalesPtBinsPP=(TH1F*)fprescalesPP->Get("hPrescalesPtBins");
  if (usePrescaleCorr==1){
    for (int i=0;i<nBins;i++) {
      hSigmaPPStat->SetBinContent(i+1,hSigmaPPStat->GetBinContent(i+1)/hPrescalesPtBinsPP->GetBinContent(i+1));
    }
  }
  Double_t xr[nBins], yr[nBins], xrlow[nBins], yrlow[nBins],xrhigh[nBins],yrhigh[nBins];
  Double_t yf[nBins], yflow[nBins],yfhigh[nBins];
  for(int i=0;i<nBins;i++)
    {
      gaeBplusReference->GetPoint(i,xr[i],yf[i]);
      xrlow[i] = gaeBplusReference->GetErrorXlow(i);
      xrhigh[i] = gaeBplusReference->GetErrorXhigh(i);
      yr[i] = hSigmaPPStat->GetBinContent(i+1)/yf[i];
      yrlow[i] = hSigmaPPStat->GetBinError(i+1)/yf[i];
      yrhigh[i] = hSigmaPPStat->GetBinError(i+1)/yf[i];
      yflow[i] = gaeBplusReference->GetErrorYlow(i)/yf[i];
      yfhigh[i] = gaeBplusReference->GetErrorYhigh(i)/yf[i];
      yf[i] = yf[i]/yf[i];
    }
    

  TGraphAsymmErrors* gaeRatio = new TGraphAsymmErrors(nBins,xr,yr,xrlow,xrhigh,yrlow,yrhigh);
  gaeRatio->SetName("gaeRatio");
  gaeRatio->SetMarkerStyle(20);
  gaeRatio->SetMarkerSize(0.8);
  TGraphAsymmErrors* gaeRatioFo = new TGraphAsymmErrors(nBins,xr,yf,xrlow,xrhigh,yflow,yfhigh);
  gaeRatioFo->SetName("gaeRatioFo");
  gaeRatioFo->SetLineWidth(2);
  gaeRatioFo->SetLineColor(2);
  gaeRatioFo->SetFillColor(2);
  gaeRatioFo->SetFillStyle(3002);
  TCanvas* cSigma = new TCanvas("cSigma","",600,750);
  cSigma->SetFrameBorderMode(0);
  cSigma->SetFrameBorderMode(0);
  cSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
  cSigma->SetFillColor(0);
  cSigma->SetBorderMode(0);
  cSigma->SetBorderSize(2);
  cSigma->SetLeftMargin(0.1451613);
  cSigma->SetRightMargin(0.05443548);
  cSigma->SetTopMargin(0.08474576);
  cSigma->SetBottomMargin(0.1165254);
  cSigma->SetFrameBorderMode(0);
  cSigma->SetFrameBorderMode(0);
  cSigma->cd();
  TPad* pSigma = new TPad("pSigma","",0,0.25,1,1);
  pSigma->SetFillColor(0);
  pSigma->SetBorderMode(0);
  pSigma->SetBorderSize(2);
  pSigma->SetLeftMargin(0.1451613);
  pSigma->SetRightMargin(0.05443548);
  pSigma->SetTopMargin(0.08474576);
  pSigma->SetBottomMargin(0);
  pSigma->SetLogy();
  pSigma->Draw();
  pSigma->cd();

  TH2F* hemptySigma=new TH2F("hemptySigma","",50,0.,110.,10.,1.1e-2,1.e6);  
  hemptySigma->GetXaxis()->CenterTitle();
  hemptySigma->GetYaxis()->CenterTitle();
  hemptySigma->GetYaxis()->SetTitle("d#sigma / dp_{T}( pb GeV^{-1}c)");
  hemptySigma->GetXaxis()->SetTitleOffset(1.);
  hemptySigma->GetYaxis()->SetTitleOffset(1.3);
  hemptySigma->GetXaxis()->SetTitleSize(0.045);
  hemptySigma->GetYaxis()->SetTitleSize(0.045);
  hemptySigma->GetXaxis()->SetTitleFont(42);
  hemptySigma->GetYaxis()->SetTitleFont(42);
  hemptySigma->GetXaxis()->SetLabelFont(42);
  hemptySigma->GetYaxis()->SetLabelFont(42);
  hemptySigma->GetXaxis()->SetLabelSize(0.04);
  hemptySigma->GetYaxis()->SetLabelSize(0.04);  
  hemptySigma->SetMaximum(2);
  hemptySigma->SetMinimum(0.);
  hemptySigma->Draw();
  gaeBplusReference->Draw("5same");
  hSigmaPPStat->Draw("epsame");  
  
  TLegend *legendSigma=new TLegend(0.5100806,0.5868644,0.8084677,0.7605932,"");
  legendSigma->SetBorderSize(0);
  legendSigma->SetLineColor(0);
  legendSigma->SetFillColor(0);
  legendSigma->SetFillStyle(1001);
  legendSigma->SetTextFont(42);
  legendSigma->SetTextSize(0.045);
  
  TLegendEntry *ent_SigmaPP=legendSigma->AddEntry(hSigmaPPStat,"pp","pf");
  ent_SigmaPP->SetTextFont(42);
  ent_SigmaPP->SetLineColor(1);
  ent_SigmaPP->SetMarkerColor(1);
  
  TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(gaeBplusReference,"FONLL pp ref.","f");
  ent_Sigmapp->SetTextFont(42);
  ent_Sigmapp->SetLineColor(5);
  ent_Sigmapp->SetMarkerColor(1);
  legendSigma->Draw("same");
    
  TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS Preliminary     pp #sqrt{s}= 5.02 TeV");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(42);
  tlatex1->SetTextSize(0.045);
  tlatex1->Draw();
    
  TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 9.97 pb^{-1}");
  tlatexlumi->SetNDC();
  tlatexlumi->SetTextColor(1);
  tlatexlumi->SetTextFont(42);
  tlatexlumi->SetTextSize(0.045);

  cSigma->cd();
  TPad* pRatio = new TPad("pRatio","",0,0,1,0.25);
  pRatio->SetLeftMargin(0.1451613);
  pRatio->SetRightMargin(0.05443548);
  pRatio->SetTopMargin(0);
  pRatio->SetBottomMargin(0.25);
  pRatio->Draw();
  pRatio->cd();

  TH2F* hemptyRatio=new TH2F("hemptyRatio","",50,0.,110.,10.,0.,3.0);
  hemptyRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyRatio->GetYaxis()->CenterTitle();
  hemptyRatio->GetYaxis()->SetTitle("Data / FONLL");
  hemptyRatio->GetXaxis()->SetTitleOffset(1.);
  hemptyRatio->GetYaxis()->SetTitleOffset(0.5);
  hemptyRatio->GetXaxis()->SetTitleSize(0.12);
  hemptyRatio->GetYaxis()->SetTitleSize(0.12);
  hemptyRatio->GetXaxis()->SetTitleFont(42);
  hemptyRatio->GetYaxis()->SetTitleFont(42);
  hemptyRatio->GetXaxis()->SetLabelFont(42);
  hemptyRatio->GetYaxis()->SetLabelFont(42);
  hemptyRatio->GetXaxis()->SetLabelSize(0.1);
  hemptyRatio->GetYaxis()->SetLabelSize(0.1);  
  hemptyRatio->Draw();

  TLine* l = new TLine(10,1,210,1);
  l->SetLineWidth(1);
  l->SetLineStyle(2);
  gaeRatioFo->Draw("5same");
  gaeRatio->Draw("epsame");
  l->Draw("same");
  cSigma->SaveAs("canvasSigmaDzeroRatio.pdf");
  
  
  TCanvas* cEffPP = new TCanvas("cEffPP","",600,500);
  
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,0.,110.,10.,0,1.);  
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
  
  cEffPP->cd();
  hemptyEff->Draw();
  hEffPP->Draw("same");
  cEffPP->SaveAs("efficiencyPP.pdf");

}


int main(int argc, char *argv[])
{
  if((argc != 5))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 5)
    CrossSectionRatio(argv[1], argv[2], argv[3],atoi(argv[4]));
  return 0;
}


