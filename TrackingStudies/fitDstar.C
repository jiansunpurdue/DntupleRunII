#include <TF1.h>
#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TLatex.h>

void fitDstar(char *infname = "/data/wangj/Data2015/Dntuple/ntD_DfinderData_pp_20151206_dPt5tkPt1_D0DsDstar3p5p.root"){

   TFile *inf = new TFile(infname);
   TTree *ntmix=(TTree*)inf->Get("ntDD0kpipipipi");
   TH1D *h = new TH1D("h","",100,0.139,0.159);
   TCut cutTrk = "";//"trk1PixelHit>=2&&trk1StripHit>=10&&trk1Chi2ndf<5&&trk2PixelHit>=2&&trk2StripHit>=10&&trk2Chi2ndf<5";
   
   TCanvas *c = new TCanvas("c","",750,750);
  
   ntmix->Draw("Dmass-DtktkResmass>>h","abs(DtktkResmass-1.86486)<0.015&&Dpt>10"&&cutTrk,"",10000000);
   h->Sumw2();
   //   TF1 *f = new TF1("f","(1-exp(-(x-[8])/[0]))*(((x)/[8])**[1]+[2]*(((x)/[8])-1))*[3]+[4]*(TMath::Voigt(x-[5],[6],[7]))");
   TF1 *f = new TF1("f","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))");
//   TF1 *f = new TF1("f","(1-exp(-(x-0.13957018)/[0]))*(((x)/0.13957018)**[1]+[2]*(((x)/0.13957018)-1))*[3]+[4]*(TMath::Gaus(x,[5],[6]))");
   //TF1 *f = new TF1("f","(1-exp(-(x-1.86486-0.13957018)/[0]))*(((x-1.86486)/0.13957018)**[1]+[2]*(((x-1.86486)/0.13957018)-1))*[3]+[4]*TMath::Gaus(x,[5],[6])");
   f->SetLineColor(4);
//   f->SetParameters(-2.3825e6,-7.99713e-1,-1.42957,-5.50069e10,5.33573,1.45491e-1,2.78677e-6,1.43145e-3,0.13957018);
//   f->SetParameters(-7.3e5,-2.2e1,5.24e-1,-7.18e9,2e2,1.45491e-1,9e-4,0.1,8e-4);
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
