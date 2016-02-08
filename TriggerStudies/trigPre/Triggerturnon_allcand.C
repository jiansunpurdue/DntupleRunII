#include "uti.h"
using namespace std;

#define Nptbins 29
double ptbins[Nptbins+1] = {0, 4, 6, 8, 10, 13, 15, 18, 20, 23, 26, 28, 30, 33, 36, 38, 40, 43, 46, 48, 50, 53, 56, 58, 60, 65, 70, 75, 80, 100};

//TString mvatk = "Dtrk1highPurity&&Dtrk2highPurity)&&(Dmass>1.75&&Dmass<1.95)";
TString mvatk = "(Dtrk1highPurity&&Dtrk2highPurity)&&(Dtrk1PtErr/Dtrk1Pt)<0.1&&(Dtrk2PtErr/Dtrk2Pt)<0.1&&(Dtrk1PixelHit+Dtrk1StripHit)>10&&(Dtrk2PixelHit+Dtrk2StripHit)>10&&(Dmass>1.75&&Dmass<1.95)";

TString mbtrg = "pPAprimaryVertexFilter&&pBeamScrapingFilter&&(HLT_L1MinimumBiasHF1OR_part0_v1||HLT_L1MinimumBiasHF1OR_part1_v1||HLT_L1MinimumBiasHF1OR_part2_v1||HLT_L1MinimumBiasHF1OR_part3_v1||HLT_L1MinimumBiasHF1OR_part4_v1||HLT_L1MinimumBiasHF1OR_part5_v1||HLT_L1MinimumBiasHF1OR_part6_v1||HLT_L1MinimumBiasHF1OR_part7_v1||HLT_L1MinimumBiasHF1OR_part8_v1||HLT_L1MinimumBiasHF1OR_part9_v1||HLT_L1MinimumBiasHF1OR_part10_v1||HLT_L1MinimumBiasHF1OR_part11_v1||HLT_L1MinimumBiasHF1OR_part12_v1||HLT_L1MinimumBiasHF1OR_part13_v1||HLT_L1MinimumBiasHF1OR_part14_v1||HLT_L1MinimumBiasHF1OR_part15_v1||HLT_L1MinimumBiasHF1OR_part16_v1||HLT_L1MinimumBiasHF1OR_part17_v1||HLT_L1MinimumBiasHF1OR_part18_v1||HLT_L1MinimumBiasHF1OR_part19_v1)";

TString prefilter;

Bool_t isPbPb = false;
Bool_t L1PS1 = false;

void Triggerturnon_allcand(TString infname = "mb.root", TString trigger="HLT_DmesonPPTrackingGlobal_Dpt20_v1", TString L1Seed = "L1_SingleJet28_BptxAND", Bool_t L1Ps1 = false)
{
  TH1::SetDefaultSumw2();

  TH1D* getYield(TTree* nt, TString triggerpass, TString triggername, TString prescale, TString variable, TString varname, TString varlatex, TString addcut="");
  void plotTurnOn(TH1D* hnominator, TH1D* hdenominator, TString triggerlegend, TString triggername, TString varname, TString varlatex);

  //decide whether to require the PS of L1 seed is 1 or not
  L1PS1 = L1Ps1;
  if( L1PS1 )
	  mbtrg = Form("%s&&%s_Prescl==1", mbtrg.Data(), L1Seed.Data() );

  prefilter = Form("(Dtrk1Algo<8&&Dtrk2Algo<8)&&(DlxyBS/DlxyBSErr)>1.&&TMath::Abs(Dy)<1.0&&(DsvpvDistance/DsvpvDisErr)>2.5&&Dtrk1Pt>3.&&Dtrk2Pt>3.&&Dchi2cl>0.10&&Dalpha<0.12&&%s&&%s",mvatk.Data(),mbtrg.Data());

  cout << "prefilter: " << prefilter << endl;

  TFile* infile = new TFile(infname);
  TTree* root = (TTree*)infile->Get("ntDkpi");
  root->AddFriend("ntHlt",infname);
  root->AddFriend("ntSkim",infname);

  TH1D* hpp_pt = getYield(root,"","","","Dpt","pt","p_{T} (GeV/c)");
  TH1D* hpp_pt_Hlt = getYield(root,Form("&&%s",trigger.Data()),Form("_%s",trigger.Data()),Form("%s_Prescl*",L1Seed.Data()),"Dpt","pt","p_{T} (GeV/c)");

  plotTurnOn(hpp_pt_Hlt,hpp_pt,trigger,Form("_%s",trigger.Data()),"pt","p_{T} (GeV/c)");
}

TH1D* getYield(TTree* nt, TString triggerpass, TString triggername, TString prescale, TString variable, TString varname, TString varlatex, TString addcut="")
{
  cout << " triggername: " << triggername << " prescale: " << prescale << endl;
  TH1D* hDistrib = new TH1D(Form("h%s_Distrib_%s",triggername.Data(),varname.Data()),Form(";D %s;Event",varlatex.Data()),Nptbins,ptbins);
  nt->Draw(Form("%s>>h%s_Distrib_%s", variable.Data(), triggername.Data(),varname.Data()),Form("%s(%s%s%s)",prescale.Data(),prefilter.Data(),addcut.Data(),triggerpass.Data()),"goff");
  hDistrib->Sumw2();
  TCanvas* cDistrib = new TCanvas(Form("c%s_Distrib_%s",triggername.Data(),varname.Data()),"",500,500);
  hDistrib->Draw();
  hDistrib->SetStats(0);
  if(isPbPb) cDistrib->SaveAs(Form("triggerturnonPlots/pbpb/c%s_Distrib_%s_%dL1PS1.pdf",triggername.Data(),varname.Data(),L1PS1));
  else cDistrib->SaveAs(Form("triggerturnonPlots/pp/c%s_Distrib_%s_%dL1PS1.pdf",triggername.Data(),varname.Data(),L1PS1));

  return hDistrib;
}

void plotTurnOn(TH1D* hnominator, TH1D* hdenominator, TString triggerlegend, TString triggername, TString varname, TString varlatex)
{
  TGraphAsymmErrors * pEffL1PS1;
  TH1D * pEffL1PS;
  
  if( L1PS1 )
  {
	  pEffL1PS1 = new TGraphAsymmErrors();
	  pEffL1PS1->Divide(hnominator, hdenominator, "cp");
	  pEffL1PS1->SetName(Form("p%s_Eff_%s",triggername.Data(),varname.Data()));
  }
  else
  {
	  pEffL1PS = (TH1D *) hnominator->Clone("pEffL1PS");
	  pEffL1PS->Divide( hnominator, hdenominator, 1.0, 1.0, "B");
	  pEffL1PS->SetName(Form("p%s_Eff_%s",triggername.Data(),varname.Data()));
  }
  

  TCanvas* cEff = new TCanvas(Form("c%s_Eff_%s",triggername.Data(),varname.Data()),"",500,500);
  gPad->SetGridx();
  gPad->SetGridy();
  TH2D* hempty=new TH2D(Form("hempty_%s_Eff_%s",triggername.Data(),varname.Data()),Form(";D^{0} %s;Pass efficiency",varlatex.Data()),20,0,100,40,0,4);
  if( L1PS1 ) hempty->GetYaxis()->SetRangeUser(0.0, 1.1);
  hempty->SetStats(0);
  hempty->Draw();
  if(L1PS1)
	  pEffL1PS1->Draw("PSAME");
  else
	  pEffL1PS->Draw("PSAME");
  TLatex* tex = new TLatex(0.18,0.96,triggerlegend);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->Draw();
  if(isPbPb) cEff->SaveAs(Form("triggerturnonPlots/pbpb/c%s_Eff_%s_%dL1PS1.pdf",triggername.Data(),varname.Data(),L1PS1));
  else cEff->SaveAs(Form("triggerturnonPlots/pp/c%s_Eff_%s_%dL1PS1.pdf",triggername.Data(),varname.Data(),L1PS1));
}
