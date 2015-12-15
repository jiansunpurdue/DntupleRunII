#include "uti.h"

TString weight = "1";

double binwidth3prong=(0.160-0.140)/60.;
double binwidth5prong=(0.160-0.140)/60.;
double minmass3prong=0.142;
double maxmass3prong=0.155;
double minmass5prong=0.142;
double maxmass5prong=0.155;

TString infnameData = "/data/wangj/Data2015/Dntuple/ntD_DfinderData_pp_20151209_dPt5tkPt1_D0DsDstar3p5p.root";
TString infnameMC = "/data/wangj/MC2015/Dntuple/PbPb/ntD_Pythia8_5020GeV_DstarD0kpipipi_755patch3_GEN_SIM_PU_20151120_Dstar5p_tkPt2_20151126_Evt_All.root";

const int nBins=1;  Double_t ptBins[nBins+1]={10.,200.};

void studydoubleratio(TString label="", Bool_t doweight=true)
{
  /*
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
  */
  void clean0 (TH1D* h);  
  TF1* fitDstar3prongs (TTree* nt, TTree* ntMC, double ptmin, double ptmax);
  TF1* fitDstar5prongs (TTree* nt, TTree* ntMC, double ptmin, double ptmax);

  TFile* infData = new TFile(infnameData.Data());
  TFile* infMC = new TFile(infnameMC.Data());

  TTree* ntData3prong = (TTree*)infData->Get("ntDD0kpipi");
  TTree* ntData5prong = (TTree*)infData->Get("ntDD0kpipipipi");
  TTree* ntMC3prong = (TTree*)infMC->Get("ntDD0kpipi");
  TTree* ntMC5prong = (TTree*)infMC->Get("ntDD0kpipipipi");
  
  TTree* HltTreeData = (TTree*) infData->Get("ntHlt");
  ntData3prong->AddFriend(HltTreeData);
  ntData5prong->AddFriend(HltTreeData);
  TTree* HltTreeMC = (TTree*) infMC->Get("ntHlt");
  ntMC3prong->AddFriend(HltTreeMC);
  ntMC5prong->AddFriend(HltTreeMC);

  TH1D* hPtMC3prong = new TH1D("hPtMC3prong","",nBins,ptBins);
  TH1D* hPtMC5prong = new TH1D("hPtMC5prong","",nBins,ptBins);
  TH1D* hPtData3prong = new TH1D("hPtData3prong","",nBins,ptBins);
  TH1D* hPtData5prong = new TH1D("hPtData5prong","",nBins,ptBins);

  for(int i=0;i<nBins;i++)
  {    
    TF1* fMC3prong = fitDstar3prongs(ntMC3prong,ntMC3prong,ptBins[i],ptBins[i+1]);
    TF1* fMC5prong = fitDstar5prongs(ntMC5prong,ntMC5prong,ptBins[i],ptBins[i+1]);
    TF1* fData3prong = fitDstar3prongs(ntData3prong,ntMC3prong,ptBins[i],ptBins[i+1]);
    TF1* fData5prong = fitDstar5prongs(ntData5prong,ntMC5prong,ptBins[i],ptBins[i+1]);

    double yieldMC3prong = fMC3prong->Integral(minmass3prong,maxmass3prong)/binwidth3prong;
    double yieldMC3prongErr = fMC3prong->Integral(minmass3prong,maxmass3prong)/binwidth3prong*fMC3prong->GetParError(0)/fMC3prong->GetParameter(0);
    double yieldData3prong = fData3prong->Integral(minmass3prong,maxmass3prong)/binwidth3prong;
    double yieldData3prongErr = fData3prong->Integral(minmass3prong,maxmass3prong)/binwidth3prong*fData3prong->GetParError(0)/fData3prong->GetParameter(0);
    double yieldMC5prong = fMC5prong->Integral(minmass5prong,maxmass5prong)/binwidth5prong;
    double yieldMC5prongErr = fMC5prong->Integral(minmass5prong,maxmass5prong)/binwidth5prong*fMC5prong->GetParError(0)/fMC5prong->GetParameter(0);
    double yieldData5prong = fData5prong->Integral(minmass5prong,maxmass5prong)/binwidth5prong;
    double yieldData5prongErr = fData5prong->Integral(minmass5prong,maxmass5prong)/binwidth5prong*fData5prong->GetParError(0)/fData5prong->GetParameter(0);

    hPtMC3prong->SetBinContent(i+1,yieldMC3prong/(ptBins[i+1]-ptBins[i]));
    hPtMC3prong->SetBinError(i+1,yieldMC3prongErr/(ptBins[i+1]-ptBins[i]));
    hPtMC5prong->SetBinContent(i+1,yieldData3prong/(ptBins[i+1]-ptBins[i]));
    hPtMC5prong->SetBinError(i+1,yieldData3prongErr/(ptBins[i+1]-ptBins[i]));
    hPtData3prong->SetBinContent(i+1,yieldMC5prong/(ptBins[i+1]-ptBins[i]));
    hPtData3prong->SetBinError(i+1,yieldMC5prongErr/(ptBins[i+1]-ptBins[i]));
    hPtData5prong->SetBinContent(i+1,yieldData5prong/(ptBins[i+1]-ptBins[i]));
    hPtData5prong->SetBinError(i+1,yieldData5prongErr/(ptBins[i+1]-ptBins[i]));
  }

  TFile *outputfile=new TFile("outputfile.root","recreate");
  outputfile->cd();
  hPtMC3prong->Write();
  hPtMC5prong->Write();
  hPtData3prong->Write();
  hPtData5prong->Write();
  outputfile->Close();
}

void clean0(TH1D* h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
  {
    if(h->GetBinContent(i)==0) h->SetBinError(i,1);
  }
}

TF1* fitDstar3prongs(TTree* nt, TTree* ntMC, Double_t ptmin, Double_t ptmax)
{

  TString seldata="abs(DtktkResmass-1.86486)<0.015&&Dpt>10.";
  static int count3p=0;
  count3p++;
  
  TCanvas* c= new TCanvas(Form("c_3p_%d",count3p),"",600,600);
  TH1D* h = new TH1D(Form("h_3p_%d",count3p),"",60,0.140,0.160);
  TF1* f = new TF1(Form("f_3p_%d",count3p),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])/(sqrt(2*3.14159)*[7])+[8]*TMath::Gaus(x,[6],[9])/(sqrt(2*3.14159)*[9]))",minmass3prong,maxmass3prong);
  nt->Project(Form("h_3p_%d",count3p),"Dmass-DtktkResmass",Form("%s*(%s&&Dpt>%f&&Dpt<%f)",weight.Data(),seldata.Data(),ptmin,ptmax));   
    
  f->SetLineColor(4);
  f->SetParameters(0,0,0,0,0,2.e+2,1.45491e-1,9.e-4,0.1,8.e-4);
  f->FixParameter(9,1.5e-3);
  f->FixParameter(6,0.145491);
  f->FixParameter(7,8.e-4);
  f->SetParLimits(8,0,1);
  h->Fit(Form("f_3p_%d",count3p),"LL");
  h->Fit(Form("f_3p_%d",count3p),"LL");
  h->Fit(Form("f_3p_%d",count3p),"LL","",minmass3prong,maxmass3prong);
  f->ReleaseParameter(6);
  f->ReleaseParameter(7);
  f->ReleaseParameter(9);
  f->SetParLimits(7,1.e-4,9.e-4);
  f->SetParLimits(9,1.e-4,9.e-4);
  h->Fit(Form("f_3p_%d",count3p),"LL","",0.142,0.148);
  h->Fit(Form("f_3p_%d",count3p),"LL","",0.142,0.16);
  h->Fit(Form("f_3p_%d",count3p),"LL","",0.142,0.16);
  h->Fit(Form("f_3p_%d",count3p),"LL","",0.141,0.16);
  h->Fit(Form("f_3p_%d",count3p),"LL","",0.141,0.16);
  h->Fit(Form("f_3p_%d",count3p),"LL","",0.141,0.16);

  TF1* background = new TF1(Form("background_3p_%d",count3p),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  background->SetParameter(0,f->GetParameter(0));
  background->SetParameter(1,f->GetParameter(1));
  background->SetParameter(2,f->GetParameter(2));
  background->SetParameter(3,f->GetParameter(3));
  background->SetParameter(4,f->GetParameter(4));
  background->SetLineColor(4);
  background->SetRange(minmass3prong,maxmass3prong);
  background->SetLineStyle(2);  //5=0
  
  TF1* mass = new TF1(Form("fmass_3p_%d",count3p),"[0]*((1-[3])*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+[3]*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
  mass->SetParameters(f->GetParameter(5),f->GetParameter(6),f->GetParameter(7),f->GetParameter(8),f->GetParameter(9));
  mass->SetParError(0,f->GetParError(5));
  mass->SetParError(1,f->GetParError(6));
  mass->SetParError(2,f->GetParError(7));
  mass->SetParError(3,f->GetParError(8));
  mass->SetParError(4,f->GetParError(9));
  mass->SetFillColor(kOrange-3);
  mass->SetFillStyle(3002);
  mass->SetLineColor(kOrange-3);
  mass->SetLineWidth(3);
  mass->SetLineStyle(2);
  
  h->SetXTitle("M_{K#pi#pi}-M_{K#pi} (GeV/c^{2})");
  h->SetYTitle("Entries / (5 MeV/c^{2})");
  h->SetStats(0);
  h->SetAxisRange(1,h->GetMaximum()*1.3,"Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  h->SetStats(0);
  h->Draw("e");

  background->Draw("same");   
  mass->SetRange(0.142,0.152);	
  mass->Draw("same");
  f->Draw("same");
  
  //TH1F* histo_copy_nofitfun = (TH1F*)h->Clone("histo_copy_nofitfun");
  //histo_copy_nofitfun->Draw("esame");
  c->SaveAs(Form("DMass_3prongs-%d.pdf",count3p));

  //h->GetFunction(Form("f%d",count3p))->Delete();
  //delete h;
  return mass;
}


TF1* fitDstar5prongs(TTree* nt, TTree* ntMC, Double_t ptmin, Double_t ptmax)
{
  TString seldata="abs(DtktkResmass-1.86486)<0.015&&Dpt>10.";
  //TString seldata="abs(DtktkResmass-1.86486)<0.015&&Dpt>10.&&Dtrk1Pt>2.&&DRestrk1Pt>2.&&DRestrk2Pt>2.&&DRestrk3Pt>2.&&DRestrk4Pt>2.";
  static int count5p=0;
  count5p++;
  
  TCanvas* c= new TCanvas(Form("c_5p_%d",count5p),"",600,600);
  TH1D* h = new TH1D(Form("h_5p_%d",count5p),"",60,0.140,0.160);
  TF1* f = new TF1(Form("f_5p_%d",count5p),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])/(sqrt(2*3.14159)*[7])+[8]*TMath::Gaus(x,[6],[9])/(sqrt(2*3.14159)*[9]))",minmass3prong,maxmass3prong);
  nt->Project(Form("h_5p_%d",count5p),"Dmass-DtktkResmass",Form("%s*(%s&&Dpt>%f&&Dpt<%f)",weight.Data(),seldata.Data(),ptmin,ptmax));   
    
  f->SetLineColor(4);
  f->SetParameters(0,0,0,0,0,2e2,1.45491e-1,9e-4,0.1,8e-4);
  f->FixParameter(9,15e-4);
  f->FixParameter(6,0.145491);
  f->FixParameter(7,8e-4);
  f->SetParLimits(8,0,1);
  h->Fit(Form("f_5p_%d",count5p),"LL");
  h->Fit(Form("f_5p_%d",count5p),"LL");
  h->Fit(Form("f_5p_%d",count5p),"LL","",minmass3prong,maxmass3prong);
  f->ReleaseParameter(6);
  f->ReleaseParameter(7);
  f->ReleaseParameter(9);
  f->SetParLimits(7,1e-4,9e-4);
  f->SetParLimits(9,1e-4,9e-4);
  h->Fit(Form("f_5p_%d",count5p),"LL","",0.142,0.148);
  h->Fit(Form("f_5p_%d",count5p),"LL","",0.142,0.16);
  h->Fit(Form("f_5p_%d",count5p),"LL","",0.142,0.16);
  h->Fit(Form("f_5p_%d",count5p),"LL","",0.141,0.16);
  h->Fit(Form("f_5p_%d",count5p),"LL","",0.141,0.16);
  h->Fit(Form("f_5p_%d",count5p),"LL","",0.141,0.16);

  TF1* background = new TF1(Form("background_5p_%d",count5p),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]");
  background->SetParameter(0,f->GetParameter(0));
  background->SetParameter(1,f->GetParameter(1));
  background->SetParameter(2,f->GetParameter(2));
  background->SetParameter(3,f->GetParameter(3));
  background->SetParameter(4,f->GetParameter(4));
  background->SetParameter(5,f->GetParameter(5));
  background->SetLineColor(4);
  background->SetRange(minmass3prong,maxmass3prong);
  background->SetLineStyle(2);  //5=0
  
  TF1* mass = new TF1(Form("fmass_5p_%d",count5p),"[0]*((1-[3])*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+[3]*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
  mass->SetParameters(f->GetParameter(5),f->GetParameter(6),f->GetParameter(7),f->GetParameter(8),f->GetParameter(9));
  mass->SetParError(0,f->GetParError(5));
  mass->SetParError(1,f->GetParError(6));
  mass->SetParError(2,f->GetParError(7));
  mass->SetParError(3,f->GetParError(8));
  mass->SetParError(4,f->GetParError(9));
  mass->SetFillColor(kOrange-3);
  mass->SetFillStyle(3002);
  mass->SetLineColor(kOrange-3);
  mass->SetLineWidth(3);
  mass->SetLineStyle(2);
    
  h->SetXTitle("M_{K#pi#pi#pi#pi}-M_{K#pi#pi#pi} (GeV/c^{2})");
  h->SetYTitle("Entries / (5 MeV/c^{2})");
  h->SetStats(0);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  h->SetStats(0);
  h->Draw("e");

  background->Draw("same");   
  mass->SetRange(0.142,0.152);	
  mass->Draw("same");
  f->Draw("same");
  
  //TH1F* histo_copy_nofitfun = ( TH1F * ) h->Clone("histo_copy_nofitfun");
  //histo_copy_nofitfun->Draw("esame");
  
  c->SaveAs(Form("DMass_5prongs-%d.pdf",count5p));
  
  //h->GetFunction(Form("f%d",count5p))->Delete();
  //delete h;
  
  return mass;
}


