#include "uti.h"

void CrossSectionRatioPbPb(TString particle="Dzero", bool isluminormalised=true)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);
  
  const int nbins=1;

  TFile* filePPReference = new TFile("../../FONLL/Results/output_pp_d0meson5_5TeV_y1_PbPbBinning.root");  
  TGraphAsymmErrors* gaeBplusReference = (TGraphAsymmErrors*)filePPReference->Get(Form("gaeSigma%s",particle.Data()));

  TFile* filepPb = new TFile(Form("ResultsD0_PbPb/PtSigma%s.root",particle.Data()));
  TH1F* hSigmapPbStat = (TH1F*)filepPb->Get("hPtSigma");
  
  Double_t xr[nbins], yr[nbins], xrlow[nbins], yrlow[nbins],xrhigh[nbins],yrhigh[nbins];
  Double_t yf[nbins], yflow[nbins],yfhigh[nbins];
  for(int i=0;i<nbins;i++)
    {
      gaeBplusReference->GetPoint(i,xr[i],yf[i]);
      xrlow[i] = gaeBplusReference->GetErrorXlow(i);
      xrhigh[i] = gaeBplusReference->GetErrorXhigh(i);
      yr[i] = hSigmapPbStat->GetBinContent(i+1)/yf[i];
      yrlow[i] = hSigmapPbStat->GetBinError(i+1)/yf[i];
      yrhigh[i] = hSigmapPbStat->GetBinError(i+1)/yf[i];
      yflow[i] = gaeBplusReference->GetErrorYlow(i)/yf[i];
      yfhigh[i] = gaeBplusReference->GetErrorYhigh(i)/yf[i];
      yf[i] = yf[i]/yf[i];
    }
  
  TGraphAsymmErrors* gaeRatio = new TGraphAsymmErrors(nbins,xr,yr,xrlow,xrhigh,yrlow,yrhigh);
  gaeRatio->SetName("gaeRatio");
  gaeRatio->SetMarkerStyle(20);
  gaeRatio->SetMarkerSize(0.8);
  TGraphAsymmErrors* gaeRatioFo = new TGraphAsymmErrors(nbins,xr,yf,xrlow,xrhigh,yflow,yfhigh);
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

  TH2F* hemptySigma=new TH2F("hemptySigma","",50,60.,100,10.,1.1,1e8);  
  hemptySigma->GetXaxis()->CenterTitle();
  hemptySigma->GetYaxis()->CenterTitle();
  hemptySigma->GetYaxis()->SetTitle("d#sigma / dp_{T}( pb GeV^{-1}c)");
  //hemptySigma->GetYaxis()->SetTitle("arbitrary units");
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
  hSigmapPbStat->Draw("epsame");  
  
  TLegend *legendSigma=new TLegend(0.5100806,0.5868644,0.8084677,0.7605932,"");
  legendSigma->SetBorderSize(0);
  legendSigma->SetLineColor(0);
  legendSigma->SetFillColor(0);
  legendSigma->SetFillStyle(1001);
  legendSigma->SetTextFont(42);
  legendSigma->SetTextSize(0.045);
  
  TLegendEntry *ent_SigmapPb=legendSigma->AddEntry(hSigmapPbStat,"pp","pf");
  ent_SigmapPb->SetTextFont(42);
  ent_SigmapPb->SetLineColor(1);
  ent_SigmapPb->SetMarkerColor(1);
  
  TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(gaeBplusReference,"FONLL pp ref.","f");
  ent_Sigmapp->SetTextFont(42);
  ent_Sigmapp->SetLineColor(5);
  ent_Sigmapp->SetMarkerColor(1);
  legendSigma->Draw("same");
    
  TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS Preliminary     PbPb #sqrt{s}= 5.02 TeV");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(42);
  tlatex1->SetTextSize(0.045);
  tlatex1->Draw();
    
  TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 9.5 #mu b^{-1}");
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

  TH2F* hemptyRatio=new TH2F("hemptyRatio","",50,60.,100,10.,0,1.0);
  hemptyRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyRatio->GetYaxis()->CenterTitle();
  hemptyRatio->GetYaxis()->SetTitle("FONLL R_{AA}");
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
  //gaeRatioFo->Draw("5same");
  gaeRatio->Draw("epsame");
  //l->Draw("same");
  cSigma->SaveAs("canvasSigmaDzeroRatio.pdf");
}
