#include "uti.h"

Bool_t isMC = false;
TString weight = "1";

const int nentries=4;
TString infname = "/data/wangj/Data2015/Dntuple/ntD_DfinderData_pp_20151206_dPt5tkPt1_D0DsDstar3p5p.root";
TString cases[nentries]={ "D03prongsData","D03prongsMC","D05prongsData","D05prongsMC"};
TString trigtree[nentries]={ "ntDD0kpipi","ntDD0kpipi","ntDD0kpipipipi","ntDD0kpipipipi"};
TString seldata="abs(DtktkResmass-1.86486)<0.015&&Dpt>10";

TString selection[nentries]={ "abs(DtktkResmass-1.86486)<0.015&&Dpt>10",
  "abs(DtktkResmass-1.86486)<0.015&&Dpt>10",
  "abs(DtktkResmass-1.86486)<0.015&&Dpt>10",
  "abs(DtktkResmass-1.86486)<0.015&&Dpt>10"};
TString fitfunction[nentries]={ "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))",
  "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))",
  "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))",
  "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[6],[9]))"};

const int nBins=1; Int_t binsIndex=1;  Double_t ptBins[nBins+1]={10,200};


void studydoubleratio(TString infname="/data/wangj/Data2015/Dntuple/ntD_DfinderData_pp_20151206_dPt5tkPt1_D0DsDstar3p5p.root", TString label="", Bool_t doweight=true,int option=0)
{
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  void clean0 (TH1D* h);
  TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax);

  TFile* inf = new TFile(infname.Data());

  TTree* nt = (TTree*) inf->Get(trigtree[option]);
  TTree* HltTree = (TTree*) inf->Get("ntHlt");
  HltTree->AddFriend(nt);
  nt->AddFriend(HltTree);

  TH1D* hPt = new TH1D("hPt","",nBins,ptBins);

  for(int i=0;i<nBins;i++)
  {
    TF1* f = fit(nt,nt,ptBins[i],ptBins[i+1]);
    double yield = f->Integral(1.7,2.0)/0.005;
    double yieldErr = f->Integral(1.7,2.0)/0.005*f->GetParError(0)/f->GetParameter(0);
    hPt->SetBinContent(i+1,yield/(ptBins[i+1]-ptBins[i]));
    hPt->SetBinError(i+1,yieldErr/(ptBins[i+1]-ptBins[i]));
  }  

  TCanvas* cPt =  new TCanvas("cPt","",600,600);
  cPt->SetLogy();
  hPt->SetXTitle("D^{0} p_{T} (GeV/c)");
  hPt->SetYTitle("Uncorrected dN(D^{0})/dp_{T}");
  hPt->Sumw2();
  hPt->Draw();
}

void clean0(TH1D* h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
  {
    if(h->GetBinContent(i)==0) h->SetBinError(i,1);
  }
}

TF1* fit(TTree* nt, TTree* ntMC, Double_t ptmin, Double_t ptmax)
{
  static int count=0;
  count++;

  TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
  TH1D* h = new TH1D(Form("h-%d",count),"",60,1.7,2.0);
  TF1* f = new TF1(Form("f%d",count),"[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*((1-[8])*TMath::Gaus(x,[6],[7])/(sqrt(2*3.14159)*[7])+[8]*TMath::Gaus(x,[6],[9])/(sqrt(2*3.14159)*[9]))", 0., 0.4);
  nt->Project(Form("h-%d",count),"Dmass",Form("%s*(%s&&Dpt>%f&&Dpt<%f)",weight.Data(),seldata.Data(),ptmin,ptmax));   

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

  f->SetLineColor(4);
  f2->SetLineColor(4); 
  f2->SetLineStyle(2);
  f3->SetLineStyle(2);
  f2->Draw("same");
  f3->SetLineColor(2);
  f3->SetFillStyle(3004); 
  f3->SetFillColor(2);
  f3->Draw("same");

  TF1* mass = new TF1(Form("fmass%d",count),"[5]*((1-[8])*TMath::Gaus(x,[6],[7])/(sqrt(2*3.14159)*[7])+[8]*TMath::Gaus(x,[6],[9])/(sqrt(2*3.14159)*[9]))");
  mass->SetParameters(f->GetParameter(5),f->GetParameter(6),f->GetParameter(7),f->GetParameter(8),f->GetParameter(9));
  mass->SetParError(5,f->GetParError(5));
  mass->SetParError(6,f->GetParError(6));
  mass->SetParError(7,f->GetParError(7));
  mass->SetParError(8,f->GetParError(8));
  mass->SetParError(9,f->GetParError(9));
  mass->SetFillColor(kOrange-3);
  mass->SetFillStyle(3002);
  mass->SetLineColor(kOrange-3);
  mass->SetLineWidth(3);
  mass->SetLineStyle(2);

  Double_t yield = mass->Integral(1.7,2.0)/0.005;
  Double_t yieldErr = mass->Integral(1.7,2.0)/0.005*mass->GetParError(0)/mass->GetParameter(0);

  TLegend* leg = new TLegend(0.65,0.58,0.82,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(h,"Data","pl");
  leg->AddEntry(f,"Fit","l");
  leg->AddEntry(mass,"D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
  leg->Draw("same");

  TLatex Tl;
  Tl.SetNDC();
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.04);
  Tl.SetTextFont(42);
  Tl.DrawLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
  Tl.DrawLatex(0.65,0.93, "pp #sqrt{s_{NN}} = 5.02 TeV");

  TLatex* tex;

  tex = new TLatex(0.22,0.78,Form("%.1f < p_{T} < %.1f GeV/c",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.22,0.83,"|y| < 1.0");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  h->GetFunction(Form("f%d",count))->Delete();
  TH1F* histo_copy_nofitfun = ( TH1F * ) h->Clone("histo_copy_nofitfun");
  histo_copy_nofitfun->Draw("esame");
  //
  if(nBins==1) c->SaveAs("ResultsD0_pp/DMass-inclusive.pdf");
  else c->SaveAs(Form("ResultsD0_pp/DMass-%d.pdf",count));

  return mass;
}

