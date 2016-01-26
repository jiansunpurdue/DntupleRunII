#include "uti.h"
#include "parameters.h"

Double_t setparam0=100.;
Double_t setparam1=1.865;
Double_t setparam2=0.03;
Double_t setparam10=0.005;
Double_t setparam8=0.1;
Double_t setparam9=0.1;
Double_t fixparam1=1.865;
Double_t minhisto=	1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

TString weight = "(1)";
TString seldata;
TString selmc;
TString collisionsystem;
TString cut="Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)&&(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5";
TString cut_recoonly="Dy>-1.&&Dy<1.&&(Dtrk1highPurity&&Dtrk2highPurity)&&Dtrk1Pt>1.5&&Dtrk2Pt>1.5&&abs(Dtrk1Eta)<2.4&&abs(Dtrk2Eta)<2.4";
TString cut_acceptance="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))";
TString selmcgen="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))";

void MCreweighting(int MCsample=1, bool doreweighting=false){

  TString inputFONLL="output_pp_d0meson_5TeV_y1.root";

  TString inputmc;
  TString outputfile;
  
  if(MCsample==0){
    inputmc="/afs/cern.ch/work/w/wangj/public/Dmeson/ntD_20151110_DfinderMC_20151110_EvtMatching_Pythia_D0pt15p0_Pthat15_TuneZ2_5020GeV_GENSIM_75x_1015_20151110_ppGlobaTrackingPPmenuHFlowpuv11_MBseed_twang-Pythia_1107.root";
  }
  if(MCsample==1){
    inputmc="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi/ntD_EvtBase_20160118_Dfinder_20151229_pp_Pythia8_prompt_D0pt15p0_Pthat15_TuneCUETP8M1_5020GeV_evtgen130_GEN_SIM_20151212_dPt1tkPt1_D0Ds.root";
  }
  
  selmc = Form("%s",cut.Data());

  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("ntDkpi");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  ntMC->AddFriend(ntGen);
  
  /*
  TTree* nthi = (TTree*)infMC->Get("ntHi");
  TTree* MCHltTree= (TTree*)infMC->Get("ntHlt"); //ntHlt //HltTree
  ntGen->AddFriend(nthi);
  MCHltTree->AddFriend(ntMC);
  MCHltTree->AddFriend(nthi);
  nthi->AddFriend(ntMC);
  nthi->AddFriend(MCHltTree);
  ntMC->AddFriend(nthi);
  */
  
  TH1D* hPtMC = new TH1D("hPtMC","",nBins,ptBins);
  TH1D* hPtMCrecoonly = new TH1D("hPtMCrecoonly","",nBins,ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",nBins,ptBins);
  TH1D* hPtGenAcc = new TH1D("hPtGenAcc","",nBins,ptBins);
    
  ntMC->Project("hPtMC","Dpt",TCut(weight)*(TCut(selmc.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMC);
  ntMC->Project("hPtMCrecoonly","Dpt",TCut(weight)*(TCut(cut_recoonly.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMCrecoonly);
  ntGen->Project("hPtGen","Gpt",TCut(weight)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGen);
  ntGen->Project("hPtGenAcc","Gpt",TCut(weight)*(TCut(cut_acceptance.Data())));
  divideBinWidth(hPtGenAcc);

  hPtMC->Sumw2();
  TH1D* hEff = (TH1D*)hPtMC->Clone("hEff");
  hEff->Divide(hPtGen);

  TH1D* hEffReco = (TH1D*)hPtMCrecoonly->Clone("hEffReco");
  hEffReco->Sumw2();
  hEffReco->Divide(hPtGen);

  TH1D* hEffAcc = (TH1D*)hPtGenAcc->Clone("hEffAcc");
  hEffAcc->Sumw2();
  hEffAcc->Divide(hEffAcc,hPtGen,1,1,"b");
  
  TH1D* hEffSelection = (TH1D*)hPtMC->Clone("hEffSelection");
  hEffSelection->Sumw2();
  hEffSelection->Divide(hEffSelection,hPtMCrecoonly,1,1,"b");
  
  
  TH1D*hmassReco=new TH1D("hmassReco","hmassReco",nbinsmasshisto,minhisto,maxhisto);
  TH1D*hmassRecoMatched=new TH1D("hmassRecoMatched","hmassReco",nbinsmasshisto,minhisto,maxhisto);
  ntMC->Draw("Dmass>>hmassReco","Dpt>20&&Dpt<30");
  ntMC->Draw("Dmass>>hmassRecoMatched","Dpt>20&&Dpt<30&&(Dgen==23333||Dgen==23344)");

  
  TFile* filePPReference = new TFile(inputFONLL.Data());  
  TGraphAsymmErrors* gaeBplusReference = (TGraphAsymmErrors*)filePPReference->Get("gaeSigmaDzero");
  TH1D* hFONLL = new TH1D("hFONLL","",nBins,ptBins);
  double x,y;
  for(int i=0;i<nBins;i++){
    gaeBplusReference->GetPoint(i,x,y);
    hFONLL->SetBinContent(i+1,y);
  }
  TH1D* hFONLLOverPt=(TH1D*)hFONLL->Clone("hFONLLOverPt");
  TH1D* hFONLLOverPtWeight=(TH1D*)hFONLL->Clone("hFONLLOverPtWeight");

  hFONLLOverPt->Divide(hPtGen);
  TF1 *myfit = new TF1("myfit","[0]+x*[1]+x*x*[2]+x*x*x*[3]", 30, 100);
  hFONLLOverPt->Fit("myfit");
  
  double par0=myfit->GetParameter(0);
  double par1=myfit->GetParameter(1);
  double par2=myfit->GetParameter(2);
  double par3=myfit->GetParameter(3);
  
  TString weightfunction=Form("(%f+%f*Gpt+%f*Gpt*Gpt+%f*Gpt*Gpt*Gpt)",par0,par1,par2,par3);
  TString weightfunctionreco=Form("(%f+%f*Dpt+%f*Dpt*Dpt+%f*Dpt*Dpt*Dpt)",par0,par1,par2,par3);
  TH1D* hPtMCWeight = new TH1D("hPtMCWeight","",nBins,ptBins);
  TH1D* hPtGenWeight = new TH1D("hPtGenWeight","",nBins,ptBins);
  hPtMCWeight->Sumw2();
  hPtGenWeight->Sumw2();

  ntMC->Project("hPtMCWeight","Dpt",TCut(weightfunctionreco)*TCut(weight)*(TCut(selmc.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMCWeight);
  ntGen->Project("hPtGenWeight","Gpt",TCut(weightfunction)*TCut(weight)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGenWeight);
    
  TH1D* hEffWeight = (TH1D*)hPtMCWeight->Clone("hEffWeight");
  hEffWeight->Sumw2();
  hEffWeight->Divide(hPtGenWeight);

  hFONLLOverPtWeight->Divide(hPtGenWeight);
  TH1D* hratioweight = (TH1D*)hEffWeight->Clone("hratioweight");
  hratioweight->Divide(hEff);
  
  TCanvas*canvasEff=new TCanvas("canvasEff","canvasEff",1000.,1000);
  canvasEff->Divide(2,2);
  canvasEff->cd(1);
  
  hEffAcc->SetXTitle("Gen p_{T}");
  hEffAcc->SetYTitle("#alpha");
  hEffAcc->SetMinimum(0);
  hEffAcc->SetMaximum(1.5);
  hEffAcc->SetTitle(";D^{0} p_{T} (GeV/c);Efficiency");
  hEffAcc->Draw();
  canvasEff->cd(2);
  hEffReco->SetMinimum(0);
  hEffReco->SetMaximum(1.5);
  hEffReco->SetTitle(";D^{0} p_{T} (GeV/c);Efficiency");
  hEffReco->SetXTitle("Gen p_{T}");
  hEffReco->SetYTitle("#alpha x #epsilon_{reco}");
  hEffReco->Draw();
  canvasEff->cd(3);
  hEffSelection->SetMinimum(0);
  hEffSelection->SetMaximum(1.5);
  hEffSelection->SetTitle(";D^{0} p_{T} (GeV/c);Efficiency");
  hEffSelection->SetXTitle("Gen p_{T}");
  hEffSelection->SetYTitle("#epsilon_{sel}");
  hEffSelection->Draw();  
  canvasEff->cd(4);
  hEff->SetMinimum(0);
  hEff->SetMaximum(1.5);
  hEff->SetTitle(";D^{0} p_{T} (GeV/c);Efficiency");
  hEff->Sumw2();
  hEff->SetXTitle("Gen p_{T}");
  hEff->SetYTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hEff->Draw();
  canvasEff->SaveAs(Form("canvasEff_%d.pdf",MCsample));

  
  TCanvas*canvas=new TCanvas("canvas","canvas",1000.,1000);
  canvas->SetLogy();
  canvas->Divide(3,4);
  canvas->cd(1);
  hPtGen->SetXTitle("Gen p_{T}");
  hPtGen->SetYTitle("#entries");
  hPtGen->SetMinimum(0.01);  
  hPtGen->Draw();
  canvas->cd(2);
  hPtMC->SetMinimum(0.01);
  hPtMC->SetXTitle("Reco matched p_{T}");
  hPtMC->SetYTitle("#entries");
  hPtMC->Draw();
  canvas->cd(3);
  hEff->SetXTitle("Gen p_{T}");
  hEff->SetYTitle("Efficiency unweighted");
  hEff->Draw();
  canvas->cd(4);
  hFONLL->SetXTitle("p_{T}");
  hFONLL->SetYTitle("FONLL, #entries");
  hFONLL->Draw("p");
  canvas->cd(5);
  hFONLLOverPt->SetXTitle("Gen p_{T}");
  hFONLLOverPt->SetYTitle("FONLL/PYTHIA ");
  hFONLLOverPt->Draw();
  canvas->cd(6);
  hFONLLOverPtWeight->SetMinimum(0.);
  hFONLLOverPtWeight->SetMaximum(3.);  
  hFONLLOverPtWeight->SetXTitle("Gen p_{T}");
  hFONLLOverPtWeight->SetYTitle("FONLL/Gen pt weighed ");
  hFONLLOverPtWeight->Draw();
  canvas->cd(7);
  hPtGenWeight->SetXTitle("Gen p_{T}");
  hPtGenWeight->SetYTitle("#entries weighted");
  hPtGenWeight->SetMinimum(0.01);  
  hPtGenWeight->Draw();
  canvas->cd(8);
  hPtMCWeight->SetXTitle("Reco matched p_{T}");
  hPtMCWeight->SetYTitle("#entries weighted");
  hPtMCWeight->Draw();
  canvas->cd(9);
  hEffWeight->SetMinimum(0);
  hEffWeight->SetMaximum(1.0);
  hEffWeight->SetTitle(";D^{0} p_{T} (GeV/c);Efficiency");
  hEffWeight->SetXTitle("Gen p_{T}");
  hEffWeight->SetYTitle("Efficiency weighted");
  hEffWeight->Draw();  
  canvas->cd(10);
  hratioweight->SetXTitle("Gen p_{T}");
  hratioweight->SetYTitle("weighted/non weighted efficiencies");
  hratioweight->SetMinimum(0.95);
  hratioweight->SetMaximum(1.05);  
  hratioweight->Draw("p");  
  canvas->SaveAs(Form("MCreweighting_%d.pdf",MCsample));
  
  TCanvas*canvasMassMatching=new TCanvas("canvasMassMatching","canvasMassMatching",500,500);
  canvasMassMatching->cd();
  hmassReco->SetLineColor(1);
  hmassRecoMatched->SetLineColor(2);
  hmassReco->Draw("");
  hmassRecoMatched->Draw("same");
  canvasMassMatching->SaveAs(Form("canvasMassMatching_%d.pdf",MCsample));
  
  outputfile=Form("testingWeight_%d.pdf",MCsample);
  TFile* outf = new TFile(outputfile.Data(),"recreate");
  outf->cd();
  hEff->Write();
  hPtGen->Write();
  hPtMC->Write();
  hEffReco->Write();
  myfit->Write();
  outf->Close();

}
