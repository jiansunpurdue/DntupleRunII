#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom2.h>
#include <TExec.h>
#include "parameters.h"


using namespace std;

void getMaximum(TH1D *hResult,TH1D *h)
{
   for (int i=0;i<=hResult->GetNbinsX()+1;i++){
      if (h->GetBinContent(i)>hResult->GetBinContent(i)) hResult->SetBinContent(i,h->GetBinContent(i));
   }
}

void getMinimum(TH1D *hResult,TH1D *h)
{
   for (int i=0;i<=hResult->GetNbinsX()+1;i++){
      if (h->GetBinContent(i)<hResult->GetBinContent(i)) hResult->SetBinContent(i,h->GetBinContent(i));
   }
}


void setHist(TH1D *h,double val)
{
   for (int i=0;i<=h->GetNbinsX();i++){
      h->SetBinContent(i,val);
      h->SetBinError(i,0);
   }
}

void normalize(TH2D* h2D)
{
   TH1D *htmp = h2D->ProjectionX("htmp");
   TH1D *htmp2 = h2D->ProjectionY("htmp2");
   
   for (int x=1;x<=h2D->GetNbinsX()+1;x++){
      if (htmp->GetBinContent(x)==0) continue;
      double ratio = 1/htmp->GetBinContent(x);
       
      for (int y=1;y<=h2D->GetNbinsY()+1;y++){
         double val = h2D->GetBinContent(x,y)*ratio;
         double valError = h2D->GetBinError(x,y)*ratio;
	 h2D->SetBinContent(x,y,val);
	 h2D->SetBinError(x,y,valError);
      }   
   }
   delete htmp;
   delete htmp2;
   
}


void reweighthisto(TH1D* hBPt, TH2D* h2D, TH1D *hRAA=0, bool weightByJpsiRAA=0, bool weightByBRAA=0)
{
   TH1D *htmp = h2D->ProjectionX("htmp");
   TH1D *htmp2 = h2D->ProjectionY("htmp2");
   
   for (int x=1;x<=h2D->GetNbinsX()+1;x++){
      if (htmp->GetBinContent(x)==0) continue;
      double ratio = hBPt->GetBinContent(x)/htmp->GetBinContent(x);
      if (weightByBRAA) ratio *= hRAA->GetBinContent(hRAA->FindBin(htmp->GetBinLowEdge(x)));   
       
      for (int y=1;y<=h2D->GetNbinsY()+1;y++){
         double ratio2 = ratio;
         if (weightByJpsiRAA) ratio2 *= hRAA->GetBinContent(hRAA->FindBin(htmp2->GetBinLowEdge(y)));   
         double val = h2D->GetBinContent(x,y)*ratio2;
         double valError = h2D->GetBinError(x,y)*ratio2;
	 h2D->SetBinContent(x,y,val);
	 h2D->SetBinError(x,y,valError);
      }   
   }
   delete htmp;
   
}
void plotFeedDown(TString inputDfile="output_pp_d0meson_5TeV_y1.root",TString inputBfile="output_pp_Bmeson_5TeV_y1.root", TString ntuplePythia="/data/HeavyFlavourRun2/BtoDPythia/treefile_ptall_11january2016.root", int ntest=1, int centL=0,int centH=100)
{
   // B cross-section
   TFile *inf = new TFile(inputBfile.Data());
//   TFile *inf = new TFile("outputBplus_D_pp_rap24.root");
//    TFile *inf = new TFile("outputBplus_pp.root");
   TH1D *hBPtMax = (TH1D*)inf->Get("hmaxall");
   TH1D *hBPtMin = (TH1D*)inf->Get("hminall");
   TH1D *hBPt = (TH1D*)inf->Get("hpt");
   hBPt->SetName("hBPt");
   hBPtMax->SetName("hBPtMax");
   hBPtMin->SetName("hBPtMin");

   TH1D *hBMaxRatio = (TH1D*)hBPt->Clone("hBMaxRatio");
   hBMaxRatio->Divide(hBPtMax);

   TH1D *hBMinRatio = (TH1D*)hBPt->Clone("hBMinRatio");
   hBMinRatio->Divide(hBPtMin);

   
   hBPt->Rebin(1); 
   
   // D cross-section
//   TFile *infD = new TFile("outputD0_D_pp.root");
   TFile *infD = new TFile(inputDfile.Data());
   TH1D *hDPtMax = (TH1D*)infD->Get("hmaxall");
   TH1D *hDPtMin = (TH1D*)infD->Get("hminall");
   TH1D *hDPt = (TH1D*)infD->Get("hpt");
   hDPt->SetName("hDPt");
   hDPtMax->SetName("hDPtMax");
   hDPtMin->SetName("hDPtMin");
   hDPt->Rebin(1); 

   // ratio of B->D0: not correct85% from PYTHIA
   //hBPt->Scale(0.85);
   hBPt->Scale(0.598);
   
   // c->D (55.7%)
   hDPt->Scale(0.557);

   
   TFile *inf2 = new TFile(ntuplePythia.Data());
//   TFile *inf2 = new TFile("test.root");
   TTree *hi = (TTree*) inf2->Get("ana/hi");

   hi->SetAlias("yD","log((sqrt(1.86484*1.86484+pt*pt*cosh(eta)*cosh(eta))+pt*sinh(eta))/sqrt(1.86484*1.86484+pt*pt))");			    
   hi->SetAlias("yB","log((sqrt(5.3*5.3+pt*pt*cosh(eta)*cosh(eta))+pt*sinh(eta))/sqrt(5.3*5.3+pt*pt))");			    
   hi->SetAlias("yJ","log((sqrt(3.09692*3.09692+pt*pt*cosh(eta)*cosh(eta))+pt*sinh(eta))/sqrt(3.09692*3.09692+pt*pt))");			    


   // 6.5, 8, 10, 13, 30

/*
   TH1D *hBNoCut = (TH1D*)hBPt->Clone("hBNoCut");
   TH1D *hBHasD  = (TH1D*)hBPt->Clone("hBHasD");
   
   hi->Draw("pt>>hBHasD","(abs(pdg)>500&&abs(pdg)<600&&abs(yB)<2.4)&&Sum$(abs(pdg)==421&&abs(yD)<2)>0");
   hi->Draw("pt>>hBNoCut","(abs(pdg)>500&&abs(pdg)<600&&abs(yB)<2.4)");
;

   hBNoCut->Divide(hBHasD);
   hBPt->Divide(hBNoCut);
  */  

// 0-100%
   int npoint = 7;
   	
   double ptBins_npjpsi[8] = {1,3,6.5,8,10,13,30,300};
  
   double raa_npjpsi[7];//      = {1,0.6, 0.52,0.43,0.43,0.34,0.5};  
   double raaStat_npjpsi[7];//  = {1,0.4,0.12,0.08,0.09,0.07,0.5};
   double raaSyst_npjpsi[7];//  = {0,0,0.06,0.05,0.05,0.04,0};

/*
0-10, 10-20, 20-30, 30-40, 40-50, 50-100
double nonPromptJpsiRAA_2012[]           = {0.,0.38,0.43,0.48,0.52,0.65,0.69};
double nonPromptJpsiRAAError_2012[]      = {0.,0.02,0.03,0.03,0.04,0.06,0.07};
double nonPromptJpsiRAAErrorSyst_2012[]  = {0.,0.04,0.05,0.05,0.06,0.07,0.07};
*/

   if (centL==0&&centH==100) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=0.6;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=0.4;	// prelim
      raa_npjpsi[2]=0.52;	raaStat_npjpsi[2]=0.12;		raaSyst_npjpsi[2]=0.06;	// np jpsi pas
      raa_npjpsi[3]=0.43;	raaStat_npjpsi[3]=0.08;		raaSyst_npjpsi[3]=0.05;	// np jpsi pas
      raa_npjpsi[4]=0.43;	raaStat_npjpsi[4]=0.09;		raaSyst_npjpsi[4]=0.05;	// np jpsi pas
      raa_npjpsi[5]=0.34;	raaStat_npjpsi[5]=0.07;		raaSyst_npjpsi[5]=0.04;	// np jpsi pas
      raa_npjpsi[6]=0.5;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.25;	// b-jet
   }
   
   if (centL==0&&centH==10) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.38;	raaStat_npjpsi[2]=0.02;		raaSyst_npjpsi[2]=0.04;	// np jpsi pas
      raa_npjpsi[3]=0.38;	raaStat_npjpsi[3]=0.02;		raaSyst_npjpsi[3]=0.04;	// np jpsi pas
      raa_npjpsi[4]=0.38;	raaStat_npjpsi[4]=0.02;		raaSyst_npjpsi[4]=0.04;	// np jpsi pas
      raa_npjpsi[5]=0.38;	raaStat_npjpsi[5]=0.02;		raaSyst_npjpsi[5]=0.04;	// np jpsi pas
      raa_npjpsi[6]=0.39;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.20;	// b-jet
   
   }

   if (centL==10&&centH==20) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.43;	raaStat_npjpsi[2]=0.03;		raaSyst_npjpsi[2]=0.05;	// np jpsi pas
      raa_npjpsi[3]=0.43;	raaStat_npjpsi[3]=0.03;		raaSyst_npjpsi[3]=0.05;	// np jpsi pas
      raa_npjpsi[4]=0.43;	raaStat_npjpsi[4]=0.03;		raaSyst_npjpsi[4]=0.05;	// np jpsi pas
      raa_npjpsi[5]=0.43;	raaStat_npjpsi[5]=0.03;		raaSyst_npjpsi[5]=0.05;	// np jpsi pas
      raa_npjpsi[6]=0.47;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.24;	// b-jet
   
   }

   if (centL==20&&centH==30) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.48;	raaStat_npjpsi[2]=0.03;		raaSyst_npjpsi[2]=0.05;	// np jpsi pas
      raa_npjpsi[3]=0.48;	raaStat_npjpsi[3]=0.03;		raaSyst_npjpsi[3]=0.05;	// np jpsi pas
      raa_npjpsi[4]=0.48;	raaStat_npjpsi[4]=0.03;		raaSyst_npjpsi[4]=0.05;	// np jpsi pas
      raa_npjpsi[5]=0.48;	raaStat_npjpsi[5]=0.03;		raaSyst_npjpsi[5]=0.05;	// np jpsi pas
      raa_npjpsi[6]=0.47;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.24;	// b-jet
   
   }

   if (centL==30&&centH==40) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.52;	raaStat_npjpsi[2]=0.04;		raaSyst_npjpsi[2]=0.06;	// np jpsi pas
      raa_npjpsi[3]=0.52;	raaStat_npjpsi[3]=0.04;		raaSyst_npjpsi[3]=0.06;	// np jpsi pas
      raa_npjpsi[4]=0.52;	raaStat_npjpsi[4]=0.04;		raaSyst_npjpsi[4]=0.06;	// np jpsi pas
      raa_npjpsi[5]=0.52;	raaStat_npjpsi[5]=0.04;		raaSyst_npjpsi[5]=0.06;	// np jpsi pas
      raa_npjpsi[6]=0.61;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.30;	// b-jet
   
   }

   if (centL==40&&centH==50) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.65;	raaStat_npjpsi[2]=0.06;		raaSyst_npjpsi[2]=0.07;	// np jpsi pas
      raa_npjpsi[3]=0.65;	raaStat_npjpsi[3]=0.06;		raaSyst_npjpsi[3]=0.07;	// np jpsi pas
      raa_npjpsi[4]=0.65;	raaStat_npjpsi[4]=0.06;		raaSyst_npjpsi[4]=0.07;	// np jpsi pas
      raa_npjpsi[5]=0.65;	raaStat_npjpsi[5]=0.06;		raaSyst_npjpsi[5]=0.07;	// np jpsi pas
      raa_npjpsi[6]=0.61;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.30;	// b-jet
   
   }

   if (centL==50&&centH==100) {
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.69;	raaStat_npjpsi[2]=0.07;		raaSyst_npjpsi[2]=0.07;	// np jpsi pas
      raa_npjpsi[3]=0.69;	raaStat_npjpsi[3]=0.07;		raaSyst_npjpsi[3]=0.07;	// np jpsi pas
      raa_npjpsi[4]=0.69;	raaStat_npjpsi[4]=0.07;		raaSyst_npjpsi[4]=0.07;	// np jpsi pas
      raa_npjpsi[5]=0.69;	raaStat_npjpsi[5]=0.07;		raaSyst_npjpsi[5]=0.07;	// np jpsi pas
      raa_npjpsi[6]=0.70;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.35;	// b-jet
   
   }


   if (centL==0&&centH==20) {  //averaged by ncoll
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.4;	raaStat_npjpsi[2]=0.03;		raaSyst_npjpsi[2]=0.05;	// np jpsi pas
      raa_npjpsi[3]=0.4;	raaStat_npjpsi[3]=0.03;		raaSyst_npjpsi[3]=0.05;	// np jpsi pas
      raa_npjpsi[4]=0.4;	raaStat_npjpsi[4]=0.03;		raaSyst_npjpsi[4]=0.05;	// np jpsi pas
      raa_npjpsi[5]=0.4;	raaStat_npjpsi[5]=0.03;		raaSyst_npjpsi[5]=0.05;	// np jpsi pas
      raa_npjpsi[6]=0.42;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.21;	// b-jet
   
   }

   if (centL==10&&centH==30) {  //averaged by ncoll
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.45;	raaStat_npjpsi[2]=0.03;		raaSyst_npjpsi[2]=0.05;	// np jpsi pas
      raa_npjpsi[3]=0.45;	raaStat_npjpsi[3]=0.03;		raaSyst_npjpsi[3]=0.05;	// np jpsi pas
      raa_npjpsi[4]=0.45;	raaStat_npjpsi[4]=0.03;		raaSyst_npjpsi[4]=0.05;	// np jpsi pas
      raa_npjpsi[5]=0.45;	raaStat_npjpsi[5]=0.03;		raaSyst_npjpsi[5]=0.05;	// np jpsi pas
      raa_npjpsi[6]=0.47;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.24;	// b-jet
   
   }


   if (centL==30&&centH==50) {  //averaged by ncoll
      raa_npjpsi[0]=1.0;	raaStat_npjpsi[0]=0.0;		raaSyst_npjpsi[0]=1.0;	// no measurement
      raa_npjpsi[1]=1.0;	raaStat_npjpsi[1]=0.0;		raaSyst_npjpsi[1]=1.0;	// no measurement
      raa_npjpsi[2]=0.57;	raaStat_npjpsi[2]=0.06;		raaSyst_npjpsi[2]=0.07;	// np jpsi pas
      raa_npjpsi[3]=0.57;	raaStat_npjpsi[3]=0.06;		raaSyst_npjpsi[3]=0.07;	// np jpsi pas
      raa_npjpsi[4]=0.57;	raaStat_npjpsi[4]=0.06;		raaSyst_npjpsi[4]=0.07;	// np jpsi pas
      raa_npjpsi[5]=0.57;	raaStat_npjpsi[5]=0.06;		raaSyst_npjpsi[5]=0.07;	// np jpsi pas
      raa_npjpsi[6]=0.61;	raaStat_npjpsi[6]=0.0;		raaSyst_npjpsi[6]=0.30;	// b-jet
   
   }


   TH1D *hNPJpsiRAA = new TH1D("hNPJpsiRAA","",npoint,ptBins_npjpsi);

   for (int i=1;i<=npoint;i++)
   {
      hNPJpsiRAA->SetBinContent(i,raa_npjpsi[i-1]);      
      hNPJpsiRAA->SetBinError(i,sqrt(raaSyst_npjpsi[i-1]*raaSyst_npjpsi[i-1]+raaStat_npjpsi[i-1]*raaStat_npjpsi[i-1]));     }

   TCanvas *cJpsiRAA = new TCanvas("cJpsiRAA","",600,600);
   cJpsiRAA->SetLogx();
   TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
   setex2->Draw();
   hNPJpsiRAA->SetXTitle("Non-prompt J/psi R_{AA} (GeV/c)");
   hNPJpsiRAA->SetXTitle("Non-prompt J/psi p_{T} (GeV/c)");
   hNPJpsiRAA->SetYTitle("R_{AA}");
   hNPJpsiRAA->Draw("e1");


   TCanvas *c = new TCanvas("c","",600,600);   

   TH2D *hJpsi= new TH2D("hJpsi","",hBPt->GetNbinsX(),hBPt->GetBinLowEdge(1),hBPt->GetBinLowEdge(hBPt->GetNbinsX()+1),
                            299*4,1,300);
   TH2D *hD= new TH2D("hD","",hBPt->GetNbinsX(),hBPt->GetBinLowEdge(1),hBPt->GetBinLowEdge(hBPt->GetNbinsX()+1),
                            299*4,1,300);
    hi->Draw("pt:BPt>>hJpsi","pdg==443&&BPt>0&&abs(yJ)<1");	
    hi->Draw("pt:BPt>>hD","abs(pdg)==421&&BPt>0&&abs(yD)<1");	
 
   hJpsi->Sumw2();   
   hD->Sumw2();   
   reweighthisto(hBPt,hD);
   reweighthisto(hBPt,hJpsi);
   hJpsi->ProjectionY()->Draw("hist");
   hD->SetLineColor(4);
   hD->SetMarkerColor(4);
   hD->ProjectionY()->Draw("hist same");
   hBPt->Draw("hist same");
   
   hJpsi->SetXTitle("B p_{T} (GeV/c)");
   hJpsi->SetYTitle("J/#psi p_{T} (GeV/c)");
   hD->SetXTitle("B p_{T} (GeV/c)");
   hD->SetYTitle("D^{0} p_{T} (GeV/c)");
   
   TCanvas *c2= new TCanvas("c2","B RAA band",600,600);
   
   TRandom2 rnd;
//   hJpsi	->ProjectionX()->Draw("hist");
   TH2D *hRAATmp = new TH2D("hRAATmp","",97,3,100,100,0,2);
   hRAATmp->SetXTitle("B p_{T} (GeV/c)");
   hRAATmp->SetYTitle("R_{AA}");
   hRAATmp->Draw();

   TCanvas *c3= new TCanvas("c3","D RAA band",600,600);
   TH2D *hDRAATmp = new TH2D("hDRAATmp","",47,3,50,100,0,2);
   hDRAATmp->SetXTitle("D^{0} p_{T} (GeV/c)");
   hDRAATmp->SetYTitle("R_{AA}");

   hDRAATmp->Draw();

   TCanvas *c4= new TCanvas("c4","B->D fraction band",600,600);
   TH2D *hBtoDTmp = new TH2D("hBtoDTmp","",47,3,50,100,0,2);
   hBtoDTmp->SetXTitle("D^{0} p_{T} (GeV/c)");
   hBtoDTmp->SetYTitle("Non-prompt D fraction");
   hBtoDTmp->Draw();


  TH1D *hDFromBPt= (TH1D*)hD->ProjectionY()->Clone("hDFromBPt");
  TH1D* hDPt_central_rebin = (TH1D*)hDPt->Rebin(nBins,"hDPt_central_rebin",ptBins);
  TH1D* BtoD_central_rebin = (TH1D*)hDFromBPt->Rebin(nBins,"BtoD_central_rebin",ptBins);
  TH1D* hfractionPPnonpromp_rebin = (TH1D*)hDFromBPt->Rebin(nBins,"fraction_central_rebin",ptBins);
  hfractionPPnonpromp_rebin->Divide(hDPt_central_rebin);
  
  TCanvas*cFractPP=new TCanvas("c","c",500,500);
  cFractPP->cd();
  hfractionPPnonpromp_rebin->SetMinimum(0);
  hfractionPPnonpromp_rebin->SetMaximum(0.5);
  hfractionPPnonpromp_rebin->SetXTitle("D^{0} p_{T} (GeV/c)");
  hfractionPPnonpromp_rebin->SetYTitle("Non-prompt D fraction in PP ");
  hfractionPPnonpromp_rebin->Draw();
  cFractPP->SaveAs("cFractPP.pdf");

   
   TH1D *hDFromBPtMax= (TH1D*)hD->ProjectionY()->Clone("hDFromBPtMax");
   TH1D *hDFromBPtMin= (TH1D*)hD->ProjectionY()->Clone("hDFromBPtMin");
   
   setHist(hDFromBPtMax,-1e10);
   setHist(hDFromBPtMin,1e10);
   
   for (int i=0;i<ntest;i++)
   {
       if (i%10==0) cout <<i<<endl;	
       TH1D *hRAASample = (TH1D*)hNPJpsiRAA->Clone(Form(	"hRAASample_%d",i));
       for (int j=1;j<=hRAASample->GetNbinsX();j++) {
          double RAA = (rnd.Rndm()*2-1)*hNPJpsiRAA->GetBinError(j)+hNPJpsiRAA->GetBinContent(j);
	  hRAASample->SetBinContent(j,RAA);
       }
       
       TH2D *hJpsiClone = (TH2D*)hJpsi->Clone(Form("hJpsiClone_%d",i));

       reweighthisto(hBPt,hJpsiClone,hRAASample,1);
       TH1D *hBRAA = hJpsiClone->ProjectionX(Form("hBRAA_%d",i));
       
       c2->cd();
       hBRAA->Divide(hBPt);
       hBRAA->SetLineWidth(3);
       hBRAA->SetLineColor(kGray);
       hBRAA->Rebin(4);
       hBRAA->Scale(1./4.);
       hBRAA->Draw("hist c same");
       
       delete hJpsiClone;
       
       TH2D *hDClone = (TH2D*)hD->Clone(Form("hDClone_%d",i));
       reweighthisto(hBPt,hDClone,hBRAA,0,1);
       
       
       TH1D *hDRAA = hDClone->ProjectionY(Form("hDRAA_%d",i));
       
       getMaximum(hDFromBPtMax,hDRAA);
       getMinimum(hDFromBPtMin,hDRAA);
       
       c3->cd();
       hDRAA->Divide(hDFromBPt);
       hDRAA->SetLineWidth(3);
       hDRAA->SetLineColor(kGray);
       hDRAA->Draw("hist c same");
       
       c4->cd();
       TH1D *hBtoDFrac = hDClone->ProjectionY(Form("hBtoDFrac_%d",i));
       
       hBtoDFrac->Divide(hDPt);
       hBtoDFrac->SetLineWidth(3);
       hBtoDFrac->SetLineColor(kGray);
       hBtoDFrac->Draw("hist same");
       
       delete hDClone;      
//       delete hBRAA;      
//       delete hDRAA;      
       
   }	   
   
   TFile *outf = new TFile(Form("BtoD-%d-%d.root",centL,centH),"recreate");

   TH1D *hDFromBPtCentral=(TH1D*)hDFromBPtMax->Clone("hDFromBPtCentral");
   hDFromBPtCentral->Add(hDFromBPtMin);
   hDFromBPtCentral->Scale(1./2);
   
   hNPJpsiRAA->Write();
   hDFromBPtMax->Write();
   hDFromBPtMin->Write();
   hDFromBPtCentral->Write();
   hDFromBPt->Write();
   hfractionPPnonpromp_rebin->Write();
   hJpsi->Write();
   hD->Write();
   hBPt->Write();
   hDPt->Write();   
   outf->Write();
   

}

int main(int argc, char *argv[])
{
  if((argc != 7))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 7)
    plotFeedDown(argv[1], argv[2], argv[3], atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
  return 0;
}

