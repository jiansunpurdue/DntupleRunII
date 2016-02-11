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

TString weight = "pthatweight";
TString selmc;

void MCefficiency(TString inputmc="/data/wangj/MC2015/Dntuple/backup/ntD_EvtBase_20160125_Dfinder_20151229_pp_Pythia8_prompt_D0_pthatweight.root", TString selmcgen="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))",TString selmcgenacceptance="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))&&abs(Gtk1eta)<2.0&&abs(Gtk2eta)<2.0&&Gtk1pt>2.0&&Gtk2pt>2.0", TString cut_recoonly="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&Dtrk1Algo>3&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11", TString cut="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&(DsvpvDistance/DsvpvDisErr)>3.5&&(DlxyBS/DlxyBSErr)>1.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&Dtrk1Algo>3&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11",TString label="PP")
{
  selmc = Form("%s",cut.Data());

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);


  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("ntDkpi");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  ntMC->AddFriend(ntGen);
  

  TTree* nthi = (TTree*)infMC->Get("ntHi");
  ntGen->AddFriend(nthi);
  nthi->AddFriend(ntMC);
  ntMC->AddFriend(nthi);
  
  TH1D* hPtMC = new TH1D("hPtMC","",nBins,ptBins);
  TH1D* hPtMCrecoonly = new TH1D("hPtMCrecoonly","",nBins,ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",nBins,ptBins);
  TH1D* hPtGenAcc = new TH1D("hPtGenAcc","",nBins,ptBins);
  TH1D* hPthat = new TH1D("hPthat","",100,0,500);
  TH1D* hPthatweight = new TH1D("hPthatweight","",100,0,500);

  ntMC->Project("hPtMC","Dpt",TCut(weight)*(TCut(selmc.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMC);
  ntMC->Project("hPtMCrecoonly","Dpt",TCut(weight)*(TCut(cut_recoonly.Data())&&"(Dgen==23333)"));
  divideBinWidth(hPtMCrecoonly);
  ntGen->Project("hPtGen","Gpt",TCut(weight)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGen);
  ntGen->Project("hPtGenAcc","Gpt",TCut(weight)*(TCut(selmcgenacceptance.Data())));
  divideBinWidth(hPtGenAcc);

  ntMC->Project("hPthat","pthat","1");
  ntMC->Project("hPthatweight","pthat",TCut(weight));

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
  
  TH2F* hemptyEff=new TH2F("hemptyEff","",50,10.,120.,10.,0,1.5);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  hemptyEff->GetYaxis()->SetTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyEff->GetXaxis()->SetTitleOffset(0.9);
  hemptyEff->GetYaxis()->SetTitleOffset(0.95);
  hemptyEff->GetXaxis()->SetTitleSize(0.05);
  hemptyEff->GetYaxis()->SetTitleSize(0.05);
  hemptyEff->GetXaxis()->SetTitleFont(42);
  hemptyEff->GetYaxis()->SetTitleFont(42);
  hemptyEff->GetXaxis()->SetLabelFont(42);
  hemptyEff->GetYaxis()->SetLabelFont(42);
  hemptyEff->GetXaxis()->SetLabelSize(0.035);
  hemptyEff->GetYaxis()->SetLabelSize(0.035);  
  hemptyEff->SetMaximum(2);
  hemptyEff->SetMinimum(0.);
  hemptyEff->Draw();

  TH2F* hemptyEffAcc=(TH2F*)hemptyEff->Clone("hemptyEffAcc");
  TH2F* hemptyEffReco=(TH2F*)hemptyEff->Clone("hemptyEffReco");
  TH2F* hemptyEffSelection=(TH2F*)hemptyEff->Clone("hemptyEffSelection");
  
  TCanvas*canvasEff=new TCanvas("canvasEff","canvasEff",1100.,1000);
  canvasEff->Divide(2,2);
  canvasEff->cd(1);
  
  hemptyEffAcc->SetYTitle("#alpha");
  hemptyEffAcc->Draw();
  hEffAcc->Draw("same");
  canvasEff->cd(2);
  hEffReco->GetYaxis()->SetTitleOffset(1.2);
  hemptyEffReco->SetYTitle("#alpha x #epsilon_{reco} ");
  hemptyEffReco->Draw();
  hEffReco->Draw("same");
  canvasEff->cd(3);
  hemptyEffSelection->SetYTitle(" #epsilon_{sel}");
  hemptyEffSelection->Draw();  
  hEffSelection->Draw("same");  
  canvasEff->cd(4);
  hemptyEff->Draw();
  hEff->Draw("same");
  canvasEff->SaveAs(Form("canvasEff_study%s.pdf",Form(label.Data())));
  
  
  TH2F* hemptyPthat=new TH2F("hemptyPthat","",50,0.,500.,10,1e-5,1e9);  
  hemptyPthat->GetXaxis()->CenterTitle();
  hemptyPthat->GetYaxis()->CenterTitle();
  hemptyPthat->GetYaxis()->SetTitle("Entries");
  hemptyPthat->GetXaxis()->SetTitle("pthat");
  hemptyPthat->GetXaxis()->SetTitleOffset(0.9);
  hemptyPthat->GetYaxis()->SetTitleOffset(0.95);
  hemptyPthat->GetXaxis()->SetTitleSize(0.05);
  hemptyPthat->GetYaxis()->SetTitleSize(0.05);
  hemptyPthat->GetXaxis()->SetTitleFont(42);
  hemptyPthat->GetYaxis()->SetTitleFont(42);
  hemptyPthat->GetXaxis()->SetLabelFont(42);
  hemptyPthat->GetYaxis()->SetLabelFont(42);
  hemptyPthat->GetXaxis()->SetLabelSize(0.035);
  hemptyPthat->GetYaxis()->SetLabelSize(0.035);  
  hemptyPthat->SetMaximum(2);
  hemptyPthat->SetMinimum(0.);

  TH2F* hemptyPthatWeighted=(TH2F*)hemptyPthat->Clone("hemptyPthatWeighted");
  hemptyPthatWeighted->GetXaxis()->SetTitle("pthat reweighted");
  
  TCanvas*canvasPthat=new TCanvas("canvasPthat","canvasPthat",1000.,500);
  canvasPthat->Divide(2,1);
  canvasPthat->cd(1);
  gPad->SetLogy();
  hemptyPthat->Draw("same");
  hPthat->Draw("same");
  canvasPthat->cd(2);
  gPad->SetLogy();
  hemptyPthatWeighted->Draw();
  hPthatweight->Draw("same");
  canvasPthat->SaveAs(Form("canvasPthat_%s.pdf",Form(label.Data())));

}

int main(int argc, char *argv[])
{
  if((argc != 7))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 7)
    MCefficiency(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
  return 0;
}

