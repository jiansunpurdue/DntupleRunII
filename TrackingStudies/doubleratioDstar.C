#include <TF1.h>
#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TLatex.h>

//Dtrk1Pt>0.5&&DRestrk1Pt>0.5&&DRestrk2Pt>0.5&&DRestrk3Pt>0.5&&DRestrk4Pt>0.5
//"(1-exp(-(x-[8])/[0]))*(((x)/[8])**[1]+[2]*(((x)/[8])-1))*[3]+[4]*(TMath::Voigt(x-[5],[6],[7]))"

const int nentries=4;
char *infname = "/data/wangj/Data2015/Dntuple/ntD_DfinderData_pp_20151206_dPt5tkPt1_D0DsDstar3p5p.root";
TString cases[nentries]={ "D03prongsData","D03prongsMC","D05prongsData","D05prongsMC"};
TString trigtree[nentries]={ "ntDD0kpipi","ntDD0kpipi","ntDD0kpipipipi","ntDD0kpipipipi"};
TString selection[nentries]={ "abs(DtktkResmass-1.86486)<0.015&&Dpt>10",
  "abs(DtktkResmass-1.86486)<0.015&&Dpt>10",
  "abs(DtktkResmass-1.86486)<0.015&&Dpt>10",
  "abs(DtktkResmass-1.86486)<0.015&&Dpt>10"};
TString fitfunction[nentries]={ "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))",
  "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))",
  "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))",
  "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))"};

void newfitDstar(int option=3){

  TFile *inf = new TFile(infname);
  TTree *ntmix=(TTree*)inf->Get(trigtree[option]);
  TH1D *h = new TH1D("h","",100,0.139,0.159);
  TCut cutTrk = "";

  TCanvas *c = new TCanvas("c","",750,750);
  ntmix->Draw("Dmass-DtktkResmass>>h",selection[option],"",10000000); 
  h->Sumw2();
  TF1 *f = new TF1("f",fitfunction[option]);
  if(trigtree[option]=="ntDD0kpipipipi"){
    f->SetLineColor(4);
    f->SetParameters(0,0,0,0,0,2e2,1.45491e-1,9e-4,0.1,8e-4);
    f->FixParameter(9,15e-4);
    f->FixParameter(6,0.145491);
    f->FixParameter(7,8e-4);
    f->SetParLimits(8,0,1);
    h->Fit("f","LL");
    h->Fit("f","LL");
    h->Fit("f","LL","",0.142,0.147);
    f->ReleaseParameter(6);
    f->ReleaseParameter(7);
    f->ReleaseParameter(9);
    f->SetParLimits(7,1e-4,9e-4);
    f->SetParLimits(9,1e-4,9e-4);
    h->Fit("f","LL","",0.142,0.148);
    h->Fit("f","LL","",0.142,0.16);
    h->Fit("f","LL","",0.142,0.16);
    h->Fit("f","LL","",0.141,0.16);
    h->Fit("f","LL","",0.141,0.16);
    h->Fit("f","LL","",0.141,0.16);

    h->SetXTitle("M_{K#pi#pi#pi#pi}-M_{K#pi#pi#pi} (GeV/c^{2})");
    h->SetYTitle("Entries");
    h->SetStats(0);
    h->SetAxisRange(1,h->GetMaximum()*1.3,"Y");
    TF1 *f2 = (TF1*)f->Clone("f2");
    f2->SetParameter(5,0);
    f2->SetRange(0.141,0.16);	   
    TF1 *f3 = (TF1*)f->Clone("f3");
    f3->SetParameter(0,0);
    f3->SetParameter(1,0);
    f3->SetParameter(2,0);
    f3->SetParameter(3,0);
    f3->SetParameter(4,0);
  }

  f->SetLineColor(4);
  f2->SetLineColor(4); 
  f2->SetLineStyle(2);
  f3->SetLineStyle(2);
  f2->Draw("same");
  f3->SetLineColor(2);
  f3->SetFillStyle(3004); 
  f3->SetFillColor(2);
  f3->Draw("same");
}
