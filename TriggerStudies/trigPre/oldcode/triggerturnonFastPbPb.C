#include "uti.h"

TString mvatk = "(Dmass>1.75&&Dmass<1.95)";
TString mbtrg="(HLT_HIL1MinimumBiasHF1AND_v1)";
TString prefilter = Form("(DlxyBS/DlxyBSErr)>2.&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dtrk1Pt>9&&Dtrk2Pt>9.&&Dchi2cl>0.10&&TMath::Cos(Dalpha)>0.92&&%s&&%s",mvatk.Data(),mbtrg.Data());
Bool_t isPbPb = true;

void triggerturnonFastPbPb(TString trigger="HLT_HIDmesonHITrackingGlobal_Dpt20_v1")
{
  TH1D* getYield(TTree* nt, TString triggerpass, TString triggername, TString prescale, TString variable, TString varname, TString varlatex, Int_t BIN_NUM, Double_t BIN_MIN, Double_t BIN_MAX, TString addcut="");
  void plotTurnOn(TH1D* hnominator, TH1D* hdenominator, TString triggerlegend, TString triggername, TString varname, TString varlatex, Int_t BIN_NUM, Double_t BIN_MIN, Double_t BIN_MAX);

  TString infname;
  infname = "/data/dmeson2015/PbPbNtuple/ntD_HIMinimumBias2_velicanu_1210.root";
  //infname = "/data/dmeson2015/Dntuple/ntD_HIForestExpress_PbPb_run262620-v6.root";            //From Express

  TFile* infile = new TFile(infname);
  TTree* root = (TTree*)infile->Get("ntDkpi");
  root->AddFriend("ntHlt",infname);

  TH1D* hpp_pt = getYield(root,"","","","Dpt","pt","p_{T} (GeV/c)",8,0,80);
  TH1D* hpp_pt_Hlt = getYield(root,Form("&&%s",trigger.Data()),Form("_%s",trigger.Data()),Form("*%s_Prescl",trigger.Data()),"Dpt","pt","p_{T} (GeV/c)",8,0,80);
  plotTurnOn(hpp_pt_Hlt,hpp_pt,trigger,Form("_%s",trigger.Data()),"pt","p_{T} (GeV/c)",8,0,80);
}

TH1D* getYield(TTree* nt, TString triggerpass, TString triggername, TString prescale, TString variable, TString varname, TString varlatex, Int_t BIN_NUM, Double_t BIN_MIN, Double_t BIN_MAX, TString addcut="")
{
  TH1D* hDistrib = new TH1D(Form("h%s_Distrib_%s",triggername.Data(),varname.Data()),Form(";D %s;Event",varlatex.Data()),BIN_NUM,BIN_MIN,BIN_MAX);
  nt->Project(Form("h%s_Distrib_%s",triggername.Data(),varname.Data()),Form("%s%s",variable.Data(),prescale.Data()),Form("%s%s%s",prefilter.Data(),addcut.Data(),triggerpass.Data()));
  hDistrib->Sumw2();
  TCanvas* cDistrib = new TCanvas(Form("c%s_Distrib_%s",triggername.Data(),varname.Data()),"",500,500);
  hDistrib->Draw();
  hDistrib->SetStats(0);
  if(isPbPb) cDistrib->SaveAs(Form("triggerturnonPlots/data/pbpb/c%s_Distrib_%s.pdf",triggername.Data(),varname.Data()));
  else cDistrib->SaveAs(Form("triggerturnonPlots/data/pp/c%s_Distrib_%s.pdf",triggername.Data(),varname.Data()));

  return hDistrib;
}

void plotTurnOn(TH1D* hnominator, TH1D* hdenominator, TString triggerlegend, TString triggername, TString varname, TString varlatex, Int_t BIN_NUM, Double_t BIN_MIN, Double_t BIN_MAX)
{
  TEfficiency* pEff = new TEfficiency(*hnominator,*hdenominator);
  pEff->SetName(Form("p%s_Eff_%s",triggername.Data(),varname.Data()));
  TCanvas* cEff = new TCanvas(Form("c%s_Eff_%s",triggername.Data(),varname.Data()),"",500,500);
  TH2D* hempty = new TH2D(Form("hempty_%s_Eff_%s",triggername.Data(),varname.Data()),Form(";D^{0} %s;Pass efficiency",varlatex.Data()),BIN_NUM,BIN_MIN,BIN_MAX,10,0,1.2);
  hempty->SetStats(0);
  hempty->Draw();
  pEff->Draw("PSAME");
  TLatex* tex = new TLatex(0.18,0.96,triggerlegend);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->Draw();
  if(isPbPb) cEff->SaveAs(Form("triggerturnonPlots/data/pbpb/c%s_Eff_%s.pdf",triggername.Data(),varname.Data()));
  else cEff->SaveAs(Form("triggerturnonPlots/data/pp/c%s_Eff_%s.pdf",triggername.Data(),varname.Data()));
}
