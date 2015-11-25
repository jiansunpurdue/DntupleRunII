#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include <cmath>
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include "TStyle.h"
#include "TFile.h"
#include "TROOT.h"

//#define BIN_NUM 220 //pPb_pt:220,pPb_y:40,pp_pt:222,pp_y:45
#define BIN_NUM 800 //pt0~200

//#define REBIN 3     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
#define REBIN 1     //pPb_pt:6,pPb_y:4,pp_pt:8,pp_y:4
//#define REBINp 4    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5
#define REBINp 2    //pPb_pt:7,pPb_y:5,pp_pt:9,pp_y:5

//#define HMIN 5      //pPb_pt:5,pPb_y:-2,pp_pt:9,pp_y:0
#define HMIN 0      //
//#define HMAX 60     //pPb_pt:55,pPb_y:2,pp_pt:120,pp_y:2.25
#define HMAX 200     //

void Bplusdsigmadpt_1ptbins(bool ispPb, double & diffXsec, double ptmin = 5, double ptmax = 100, string inf = "../fonllData/B_pp_pt_rap24_5p5_pt0to200.dat")
{
  TString infile,outfile;
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  if(ispPb) {
    infile="../FONLLInputs/fo_pPb_pt_rap24.dat";
    outfile="../outputBplus.root";
  }
  else {
    //infile="../FONLLInputs/fo_pp_pt_rap24.dat";
    //infile="../fonllData/B_pp_pt_rap24_5p5_pt0to200.dat";
    //infile="../fonllData/D_pp_pt_rap24_5p5_pt0to200.dat";
    infile = inf.c_str();
    //outfile="outputBplus_pp_1ptbins.root";
  }
  
  ifstream getdata(infile.Data());

  if(!getdata.is_open())
    {
      cout<<"Opening the file fails"<<endl;
    }

//  float central[BIN_NUM];
//  float min_all[BIN_NUM],max_all[BIN_NUM],min_sc[BIN_NUM],max_sc[BIN_NUM],min_mass[BIN_NUM],max_mass[BIN_NUM],min_pdf[BIN_NUM],max_pdf[BIN_NUM];
  float central[BIN_NUM+1];
  float min_all[BIN_NUM+1],max_all[BIN_NUM+1],min_sc[BIN_NUM+1],max_sc[BIN_NUM+1],min_mass[BIN_NUM+1],max_mass[BIN_NUM+1],min_pdf[BIN_NUM+1],max_pdf[BIN_NUM+1];
  int i;
  float tem;
//  for(i=0;i<BIN_NUM;i++)
  for(i=0;i<BIN_NUM+1;i++)
    {
      getdata>>tem;
      getdata>>central[i];
	  //cout<<central[i]<<endl;
      getdata>>min_all[i];
      getdata>>max_all[i];
      getdata>>min_sc[i];
      getdata>>max_sc[i];
      getdata>>min_mass[i];
      getdata>>max_mass[i];
//      getdata>>min_pdf[i];
//      getdata>>max_pdf[i];
	  //printf("%f\n",central[i]);
    }
  
  TH1F* hpt = new TH1F("hpt","",BIN_NUM,HMIN,HMAX);
  TH1F* hminall = new TH1F("hminall","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxall = new TH1F("hmaxall","",BIN_NUM,HMIN,HMAX);
  TH1F* hminsc = new TH1F("hminsc","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxsc = new TH1F("hmaxsc","",BIN_NUM,HMIN,HMAX);
  TH1F* hminmass = new TH1F("hminmass","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxmass = new TH1F("hmaxmass","",BIN_NUM,HMIN,HMAX);
  TH1F* hminpdf = new TH1F("hminpdf","",BIN_NUM,HMIN,HMAX);
  TH1F* hmaxpdf = new TH1F("hmaxpdf","",BIN_NUM,HMIN,HMAX);

  for(i=0;i<BIN_NUM;i++)
    {
//      hpt->SetBinContent(i+1,central[i]);
//      hminall->SetBinContent(i+1,min_all[i]);
//      hmaxall->SetBinContent(i+1,max_all[i]);
//      hminsc->SetBinContent(i+1,min_sc[i]);
//      hmaxsc->SetBinContent(i+1,max_sc[i]);
//      hminmass->SetBinContent(i+1,min_mass[i]);
//      hmaxmass->SetBinContent(i+1,max_mass[i]);
//      hminpdf->SetBinContent(i+1,min_pdf[i]);
//      hmaxpdf->SetBinContent(i+1,max_pdf[i]);
      hpt->SetBinContent(i+1,(central[i]+central[i+1])/2);
      hminall->SetBinContent(i+1,(min_all[i]+min_all[i+1])/2);
      hmaxall->SetBinContent(i+1,(max_all[i]+max_all[i+1])/2);
      hminsc->SetBinContent(i+1,(min_sc[i]+min_sc[i+1])/2);
      hmaxsc->SetBinContent(i+1,(max_sc[i]+max_sc[i+1])/2);
      hminmass->SetBinContent(i+1,(min_mass[i]+min_mass[i+1])/2);
      hmaxmass->SetBinContent(i+1,(max_mass[i]+max_mass[i+1])/2);
//      hminpdf->SetBinContent(i+1,(min_pdf[i]+min_pdf[i+1])/2);
//      hmaxpdf->SetBinContent(i+1,(max_pdf[i]+max_pdf[i+1])/2);
    }
  //Rebin Edge
//  double rebin[REBINp] = {10,15,20,60};
//  double rebin[REBINp] = {10,60};
//  double rebin[REBINp] = {8,60};
//  double rebin[REBINp] = {5,100};
//  double rebin[REBINp] = {0,100};
	double rebin[REBINp] = {ptmin, ptmax};

  TH1F* hpt_rebin = (TH1F*)hpt->Rebin(REBIN,"hpt_rebin",rebin);
  TH1F* hminall_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminall_rebin",rebin);
  TH1F* hmaxall_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxall_rebin",rebin);
  TH1F* hminsc_rebin = (TH1F*)hminsc->Rebin(REBIN,"hminsc_rebin",rebin);
  TH1F* hmaxsc_rebin = (TH1F*)hmaxsc->Rebin(REBIN,"hmaxsc_rebin",rebin);
  TH1F* hminmass_rebin = (TH1F*)hminmass->Rebin(REBIN,"hminmass_rebin",rebin);
  TH1F* hmaxmass_rebin = (TH1F*)hmaxmass->Rebin(REBIN,"hmaxmass_rebin",rebin);
  TH1F* hminpdf_rebin = (TH1F*)hminpdf->Rebin(REBIN,"hminpdf_rebin",rebin);
  TH1F* hmaxpdf_rebin = (TH1F*)hmaxpdf->Rebin(REBIN,"hmaxpdf_rebin",rebin);

  //bin middle
//  double apt[REBIN] = {12.5,17.5,40.};//pPb_pt
//  double apt[REBIN] = {35.};//pt10~60
//  double apt[REBIN] = {34.};//pt8~60
//  double apt[REBIN] = {52.5};//pt5~100
//  double apt[REBIN] = {50.};//pt0~100
  double apt[REBIN] = {(ptmin+ptmax)/2};

  //bin half width
//  double aptl[REBIN] = {2.5,2.5,20.};//pPb_pt
//  double aptl[REBIN] = {25.};//pt10~60
//  double aptl[REBIN] = {26.};//pt8~60
//  double aptl[REBIN] = {47.5};//pt5~100
//  double aptl[REBIN] = {50.};//pt0~100
  double aptl[REBIN] = {(ptmax-ptmin)/2};
  double asigma[REBIN],aminall[REBIN],amaxall[REBIN],aminsc[REBIN],amaxsc[REBIN],aminmass[REBIN],amaxmass[REBIN],aminpdf[REBIN],amaxpdf[REBIN],aerrorl[REBIN],aerrorh[REBIN];

  //number of every rebined bin
//  double bin_num[REBIN] = {20,20,160};//pPb_pt
//  double bin_num[REBIN] = {200};//pt10~60
//  double bin_num[REBIN] = {208};//pt8~60
//  double bin_num[REBIN] = {380};//pt5~100
//  double bin_num[REBIN] = {400};//pt0~100
  double bin_num[REBIN] = {(ptmax-ptmin)*4};
  
  int j;
  double norm=1.;
  
  for(j=0;j<REBIN;j++)
    {

      tem = hpt_rebin->GetBinContent(j+1);
      asigma[j] = tem*norm/bin_num[j];

      tem = hminall_rebin->GetBinContent(j+1);
      aminall[j] = tem*norm/bin_num[j];

      tem = hmaxsc_rebin->GetBinContent(j+1);
      amaxall[j] = tem*norm/bin_num[j];

      tem = hminsc_rebin->GetBinContent(j+1);
      aminsc[j] = tem*norm/bin_num[j];

      tem = hmaxsc_rebin->GetBinContent(j+1);
      amaxsc[j] = tem*norm/bin_num[j];

      tem = hminmass_rebin->GetBinContent(j+1);
      aminmass[j] = tem*norm/bin_num[j];

      tem = hmaxmass_rebin->GetBinContent(j+1);
      amaxmass[j] = tem*norm/bin_num[j];

      tem = hminpdf_rebin->GetBinContent(j+1);
      aminpdf[j] = tem*norm/bin_num[j];

      tem = hmaxpdf_rebin->GetBinContent(j+1);
      amaxpdf[j] = tem*norm/bin_num[j];

      aerrorl[j] = asigma[j]-aminall[j];//all,sc,mass,pdf
      aerrorh[j] = amaxall[j]-asigma[j];//all,sc,mass,pdf
    }
  //cout<<"------- pPb_5.02------"<<endl;
  cout<<endl;
 
  TGraphAsymmErrors* gae = new TGraphAsymmErrors(REBIN, apt, asigma, aptl, aptl, aerrorl, aerrorh);
  gae->SetTitle("p_{T}(GeV/c);d#sigma (B admix) /dp_{T}(pb c/GeV)");
  gae->SetFillColor(2);
  gae->SetFillStyle(3001);

  TGraphAsymmErrors* gaeSigmaDecay=(TGraphAsymmErrors*)gae->Clone();
  gaeSigmaDecay->SetName("gaeSigmaDecay");
  //double BRchain=6.09604e-5;
  double BRchain=1.;
  //double Fraction=0.401;
  double Fraction=1;
  //double norm=208.;
  norm=1.;
  double BRFraction=BRchain*Fraction;
  
  for (int k=0;k<gaeSigmaDecay->GetN();k++){
    gaeSigmaDecay->GetY()[k] *= BRFraction*norm;
	printf("diff xsec pt : %.3f ~ %.3f= %.3f pb\n", -gaeSigmaDecay->GetErrorXlow(k)+gaeSigmaDecay->GetX()[k], gaeSigmaDecay->GetErrorXhigh(k)+gaeSigmaDecay->GetX()[k], gaeSigmaDecay->GetY()[k]);
	printf("total xsec pt : %.3f pb\n", gaeSigmaDecay->GetY()[k]*(gaeSigmaDecay->GetErrorXhigh(k) + gaeSigmaDecay->GetErrorXlow(k)));
	diffXsec = gaeSigmaDecay->GetY()[k];
    gaeSigmaDecay->SetPointEYhigh(k,gaeSigmaDecay->GetErrorYhigh(k)*BRFraction*norm);
    gaeSigmaDecay->SetPointEYlow(k,gaeSigmaDecay->GetErrorYlow(k)*BRFraction*norm); 
  }
  
  /*
  TCanvas* cr = new TCanvas("cr","cr",600,500);
  cr->SetLogy();
  //TH2F* hempty=new TH2F("hempty","",10,5,60.,10.,10,100000000);  
  TH2F* hempty=new TH2F("hempty","",1,ptmin,ptmax,10.,10,100000000);  
  hempty->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hempty->GetYaxis()->SetTitle("d#sigma(B admix)/dp_{T}(pb/GeV)");
  hempty->GetXaxis()->SetTitleOffset(1.);
  hempty->GetYaxis()->SetTitleOffset(.9);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);  
  hempty->Draw();
  hminall->SetLineColor(2);
  hmaxall->SetLineColor(2);
  hpt->SetLineColor(2);
  hminall->Draw("same");
  hmaxall->Draw("same");
  hpt->Draw("same");
  gae->SetLineWidth(3);
  gae->Draw("psame");
  
  //TLatex * tlatex=new TLatex(0.18,0.85,"pp collisions at 5.02 from FONLL, |y_{LAB}|<2.4");
  TLatex * tlatex=new TLatex(0.18,0.85,"pp collisions from FONLL, |y_{LAB}|<2.4");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();
  tlatex=new TLatex(0.18,0.80,"Total syst uncertainties shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.04);
  tlatex->Draw();
  if(ispPb)cr->SaveAs("Plots/cBmesonPredFONLLBplusBinning.eps");
  else cr->SaveAs("Plots/cBmesonPredFONLLBplusBinning_pp.eps");

  gaeSigmaDecay->SetFillColor(2);
  gaeSigmaDecay->SetFillStyle(3001); 
  gaeSigmaDecay->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) #times A (GeV^{-1}c)");
   
  hempty=new TH2F("hempty","",10,0,70.,10.,1000.,5000000000);  
  hempty->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hempty->GetXaxis()->SetTitleOffset(1.);
  hempty->GetYaxis()->SetTitleOffset(.9);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.035);  
  hempty->GetYaxis()->SetTitle("d#sigma/dp_{T}(B^{+}) #times A (pb c/GeV)");

  TCanvas*canvas=new TCanvas("canvas","canvas",600,500);
  canvas->SetLogy();
  hempty->Draw();
  gaeSigmaDecay->Draw("psame");
  
  //TLatex * tlatex=new TLatex(0.2,0.85,"B^{+}=40.1%, |y|<1.93, BR unc not shown");
  tlatex=new TLatex(0.2,0.85,"B^{+},|y_{LAB}|<2.4, BR unc not shown");
  tlatex->SetNDC();
  tlatex->SetTextColor(1);
  tlatex->SetTextFont(42);
  tlatex->SetTextSize(0.05);
  tlatex->Draw();
  
  gae->SetName("gaeBplus");
  gaeSigmaDecay->SetName("gaeSigmaDecayBplus");
  if(ispPb){
    canvas->SaveAs("Plots/canvasBplus.eps");
    canvas->SaveAs("Plots/canvasBplus.pdf");
  }
  else{
    canvas->SaveAs("Plots/canvasBplus_pp_1ptbins.eps");
    canvas->SaveAs("Plots/canvasBplus_pp_1ptbins.pdf");
  }
  TFile*foutput=new TFile(outfile.Data(),"recreate");
  foutput->cd();
  gae->Write();
  gaeSigmaDecay->Write();
  */
  
}
