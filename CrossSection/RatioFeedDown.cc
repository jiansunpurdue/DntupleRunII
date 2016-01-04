#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>

#include "parameters.h"

using namespace std;


void RatioFeedDown(TString inputDpromptFONLL="output_pp_d0meson5_5TeV_y1.root", TString inputnonpromptDFONLL="output_pp_Btod0meson5_5TeV_y1.root", TString outputinclusiveD="output_inclusiveDd0meson5_5TeV_y1.root"){

  TFile*fileprompt=new TFile(inputDpromptFONLL.Data());
  TFile*fileFeedDown=new TFile(inputnonpromptDFONLL.Data());
  
  TGraphAsymmErrors *fprompt = (TGraphAsymmErrors*)fileprompt->Get("gaeSigmaDzero");
  TGraphAsymmErrors *ffeed = (TGraphAsymmErrors*)fileFeedDown->Get("gaeSigmaDzero");

  //bin middle
  double apt[nBins];
  //bin half width
  double aptl[nBins];
  //number of every rebined bin
  double bin_num[nBins];
  
  //feed/prompt ratio
  double asigmaRatio[nBins];
  double aerrorRatiol[nBins];
  double aerrorRatioh[nBins];
  //feed+prompt cross section
  double asigmaInclusive[nBins];
  double aerrorInclusivel[nBins];
  double aerrorInclusiveh[nBins];
  //dummy 1 ratio
  double asigmaOne[nBins];
  //relative errors for prompt, feed down and inclusive
  double aerrorRatioRelPromptl[nBins];
  double aerrorRatioRelPrompth[nBins];
  double aerrorRatioRelFeedl[nBins];
  double aerrorRatioRelFeedh[nBins];
  double aerrorRatioRelInclusivel[nBins];
  double aerrorRatioRelInclusiveh[nBins];
  
  for (int ibin=0; ibin<nBins; ibin++){
    apt[ibin]=(ptBins[ibin+1]+ptBins[ibin])/2.;
    aptl[ibin] = (ptBins[ibin+1]-ptBins[ibin])/2;
    bin_num[ibin]=aptl[ibin]/binsize*2;
    
    asigmaRatio[ibin]=0.;
    aerrorRatiol[ibin]=0.;
    aerrorRatioh[ibin]=0.;
    asigmaInclusive[ibin]=0.;
    aerrorInclusivel[ibin]=0.;
    aerrorInclusiveh[ibin]=0.;
    asigmaOne[ibin]=1.;
    aerrorRatioRelPromptl[ibin]=0.;
    aerrorRatioRelPrompth[ibin]=0.;
    aerrorRatioRelFeedl[ibin]=0.;
    aerrorRatioRelFeedh[ibin]=0.;
    aerrorRatioRelInclusivel[ibin]=0.;
    aerrorRatioRelInclusiveh[ibin]=0.;
  }

  double temp=0.;
  double yprompt=0.;
  double yfeed=0.;

  double yprompt_errup=0.;
  double yprompt_errlow=0.;
  double yfeed_errup=0.;
  double yfeed_errlow=0.;
    
  for(int i=0; i<nBins; i++){

    //prompt and feed down and total cross section  
    fprompt->GetPoint(i,temp,yprompt);
    ffeed->GetPoint(i,temp,yfeed);

    yprompt_errup=fprompt->GetErrorYhigh(i);
    yprompt_errlow=fprompt->GetErrorYlow(i);
    yfeed_errup=ffeed->GetErrorYhigh(i);
    yfeed_errlow=ffeed->GetErrorYlow(i);
    
    asigmaRatio[i]=yfeed/yprompt;    
    aerrorRatiol[i]=0.;    
    aerrorRatioh[i]=0.;    
    
    asigmaInclusive[i]=yfeed+yprompt;    
    aerrorInclusivel[i]=TMath::Sqrt(yprompt_errlow*yprompt_errlow+yfeed_errlow*yfeed_errlow);
    aerrorInclusiveh[i]=TMath::Sqrt(yprompt_errup*yprompt_errup+yfeed_errup*yfeed_errup);

    //relative errors  

    aerrorRatioRelPromptl[i]=yprompt_errlow/yprompt;    
    aerrorRatioRelPrompth[i]=yprompt_errup/yprompt;    
            
    aerrorRatioRelFeedl[i]=yfeed_errlow/yfeed;    
    aerrorRatioRelFeedh[i]=yfeed_errup/yfeed;    

    aerrorRatioRelInclusivel[i]=aerrorInclusivel[i]/asigmaInclusive[i];    
    aerrorRatioRelInclusiveh[i]=aerrorInclusiveh[i]/asigmaInclusive[i];    

  }
  
  TGraphAsymmErrors* gaeInclusive = new TGraphAsymmErrors(nBins, apt, asigmaInclusive, aptl, aptl, aerrorInclusivel, aerrorInclusiveh);  
  TGraphAsymmErrors* gaeRatio = new TGraphAsymmErrors(nBins, apt, asigmaRatio, aptl, aptl, aerrorRatiol, aerrorRatiol);  
  TGraphAsymmErrors* gaeErrPrompt = new TGraphAsymmErrors(nBins, apt, asigmaOne, aptl, aptl, aerrorRatioRelPromptl, aerrorRatioRelPrompth);  
  TGraphAsymmErrors* gaeErrFeed = new TGraphAsymmErrors(nBins, apt, asigmaOne, aptl, aptl, aerrorRatioRelFeedl, aerrorRatioRelFeedh);  
  TGraphAsymmErrors* gaeErrInclusive = new TGraphAsymmErrors(nBins, apt, asigmaOne, aptl, aptl, aerrorRatioRelInclusivel, aerrorRatioRelInclusiveh);  
  
  TCanvas* cFeedOverPrompt = new TCanvas("cFeedOverPrompt","cFeedOverPrompt",600,500);
  TH2F* hemptyratio=new TH2F("hemptyratio","",10,0,100.,10.,0.,0.5);  
  hemptyratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyratio->GetYaxis()->SetTitle("Feed down/prompt ratio)");
  hemptyratio->GetXaxis()->SetTitleOffset(1.);
  hemptyratio->GetYaxis()->SetTitleOffset(1.2);
  hemptyratio->GetXaxis()->SetTitleSize(0.045);
  hemptyratio->GetYaxis()->SetTitleSize(0.045);
  hemptyratio->GetXaxis()->SetTitleFont(42);
  hemptyratio->GetYaxis()->SetTitleFont(42);
  hemptyratio->GetXaxis()->SetLabelFont(42);
  hemptyratio->GetYaxis()->SetLabelFont(42);
  hemptyratio->GetXaxis()->SetLabelSize(0.04);
  hemptyratio->GetYaxis()->SetLabelSize(0.04);  
  hemptyratio->Draw();
  gaeRatio->SetLineWidth(3);
  gaeRatio->Draw("psame");
  
  TLatex * tlatexratio=new TLatex(0.22,0.85,"Feed-down / prompt, |y|<1, no uncertainty");
  tlatexratio->SetNDC();
  tlatexratio->SetTextColor(1);
  tlatexratio->SetTextFont(42);
  tlatexratio->SetTextSize(0.04);
  tlatexratio->Draw();
  cFeedOverPrompt->SaveAs("cFeedOverPrompt.pdf");

  TCanvas* cInclusiveCross = new TCanvas("cInclusiveCross","cInclusiveCross",600,500);
  cInclusiveCross->SetLogy();
  TH2F* hemptyinclusive=new TH2F("hemptyinclusive","",10,0,100.,10.,1.,100);  
  hemptyinclusive->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyinclusive->GetYaxis()->SetTitle("d#sigma(D)/dp_{T}(pb GeV-1c)");
  hemptyinclusive->GetXaxis()->SetTitleOffset(1.);
  hemptyinclusive->GetYaxis()->SetTitleOffset(1.1);
  hemptyinclusive->GetXaxis()->SetTitleSize(0.045);
  hemptyinclusive->GetYaxis()->SetTitleSize(0.045);
  hemptyinclusive->GetXaxis()->SetTitleFont(42);
  hemptyinclusive->GetYaxis()->SetTitleFont(42);
  hemptyinclusive->GetXaxis()->SetLabelFont(42);
  hemptyinclusive->GetYaxis()->SetLabelFont(42);
  hemptyinclusive->GetXaxis()->SetLabelSize(0.04);
  hemptyinclusive->GetYaxis()->SetLabelSize(0.04);  
  hemptyinclusive->Draw();
  gaeInclusive->SetLineWidth(3);
  gaeInclusive->Draw("psame");
  
  TLatex * tlatexinclusive=new TLatex(0.22,0.85,"Total prompt+feed down production, |y|<1");
  tlatexinclusive->SetNDC();
  tlatexinclusive->SetTextColor(1);
  tlatexinclusive->SetTextFont(42);
  tlatexinclusive->SetTextSize(0.04);
  tlatexinclusive->Draw();
  cInclusiveCross->SaveAs("cInclusiveCross.pdf");


  TCanvas* crelativeuncertainty = new TCanvas("crelativeuncertainty","crelativeuncertainty",600,500);
  TH2F* hemptyunc=new TH2F("hemptyunc","",10,0,100.,10.,0.7,1.3);  
  hemptyunc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyunc->GetYaxis()->SetTitle("Relative uncertainty");
  hemptyunc->GetXaxis()->SetTitleOffset(1.);
  hemptyunc->GetYaxis()->SetTitleOffset(.9);
  hemptyunc->GetXaxis()->SetTitleSize(0.045);
  hemptyunc->GetYaxis()->SetTitleSize(0.045);
  hemptyunc->GetXaxis()->SetTitleFont(42);
  hemptyunc->GetYaxis()->SetTitleFont(42);
  hemptyunc->GetXaxis()->SetLabelFont(42);
  hemptyunc->GetYaxis()->SetLabelFont(42);
  hemptyunc->GetXaxis()->SetLabelSize(0.04);
  hemptyunc->GetYaxis()->SetLabelSize(0.04);  
  hemptyunc->Draw();
  gaeErrInclusive->SetLineWidth(7);
  gaeErrInclusive->SetLineColor(1);
  gaeErrInclusive->SetMarkerColor(1);
  gaeErrInclusive->Draw("psame");
  gaeErrPrompt->SetLineWidth(3);
  gaeErrPrompt->SetLineColor(2);
  gaeErrPrompt->SetMarkerColor(2);
  gaeErrPrompt->Draw("psame");
  gaeErrFeed->SetLineWidth(3);
  gaeErrFeed->SetLineColor(4);
  gaeErrFeed->SetMarkerColor(4);
  gaeErrFeed->Draw("psame");
  
  TLegend* legend= new TLegend(0.1862416,0.7987288,0.8741611,0.9216102);
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  TLegendEntry* entprompt=legend->AddEntry(gaeErrPrompt,"Total relative error on prompt cross","PL");
  TLegendEntry* entfeed=legend->AddEntry(gaeErrFeed,"Total relative error on feed-down cross","PL");
  TLegendEntry* entinclusive=legend->AddEntry(gaeErrInclusive,"Total relative error on inclusive cross","PL");
  legend->Draw();
  crelativeuncertainty->SaveAs("crelativeuncertainty.pdf");
  
  TFile*foutput=new TFile(outputinclusiveD.Data(),"recreate");
  foutput->cd();
  gaeInclusive->SetName("gaeSigmaDzero");
  gaeInclusive->Write();
  
}

int main(int argc, char *argv[])
{
  if((argc != 4))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 4)
    RatioFeedDown(argv[1], argv[2], argv[3]);
  return 0;
}

