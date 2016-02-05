#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>

TCut massCut="abs(mass-5.279)<0.04";
TString weight = "pthatweight";
TString seldata;
TString selmc;
TString collisionsystem;
TString mass="(Dmass>1.85&&Dmass<1.88)";


void dataMC(TString inputMCprompt="/data/wangj/MC2015/Dntuple/backup/ntD_EvtBase_20160125_Dfinder_20151229_pp_Pythia8_prompt_D0_pthatweight.root", TString inputMCnonprompt="/data/wangj/MC2015/Dntuple/pp/ntD_pp_Dzero_kpi_nonprompt/ntD_EvtBase_20160203_Dfinder_20160201_pp_Pythia8_nonprompt_D0_dPt0tkPt0p5_pthatweight.root", TString inputData="/data/yjlee/dmeson/2015/trigger/mb.root",TString trgselection="1",  TString cut="Dy>-1.&&Dy<1.&&Dtrk1highPurity&&Dtrk2highPurity&&Dtrk1Pt>2.0&&Dtrk2Pt>2.0&&(DsvpvDistance/DsvpvDisErr)>3.5&&(DlxyBS/DlxyBSErr)>1.5&&Dchi2cl>0.05&&Dalpha<0.12&&Dtrk1PtErr/Dtrk1Pt<0.1&&Dtrk2PtErr/Dtrk2Pt<0.1&&abs(Dtrk1Eta)<2.0&&abs(Dtrk2Eta)<2.0&&Dtrk1Algo>3&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11", TString selmcgen="((GisSignal==1||GisSignal==2)&&(Gy>-1&&Gy<1))", TString outputfile="mytest.root"){

    TString va="DsvpvDistance/DsvpvDisErr";
    TString title ="DsvpvDistance/DsvpvDisErr";
    TString titleplot ="displacement";
    bool islog=true;
    
    seldata = Form("%s&&%s&&%s",trgselection.Data(),cut.Data(),mass.Data());
    selmc = Form("%s&&%s&&(Dgen==23333)",cut.Data(),mass.Data());

    int nBin=20;
    double BinL=0;
    double BinH=10;
    
  TFile* infMCprompt = new TFile(inputMCprompt.Data());
  TTree* ntMCprompt = (TTree*)infMCprompt->Get("ntDkpi");
  TTree* nthiMCprompt = (TTree*)infMCprompt->Get("ntHi");
  ntMCprompt->AddFriend(nthiMCprompt);

  TFile* infMCnonprompt = new TFile(inputMCnonprompt.Data());
  TTree* ntMCnonprompt = (TTree*)infMCnonprompt->Get("ntDkpi");
  TTree* nthiMCnonprompt = (TTree*)infMCnonprompt->Get("ntHi");
  ntMCnonprompt->AddFriend(nthiMCnonprompt);

  TFile* infData = new TFile(inputData.Data());
  TTree* ntData = (TTree*)infData->Get("ntDkpi");
  TTree* nthiData = (TTree*)infData->Get("ntHi");
  ntData->AddFriend(nthiMCprompt);

  TH1D *hData = new TH1D("hData","",nBin,BinL,BinH);
  TH1D *hMC = new TH1D("hMC","",nBin,BinL,BinH);
  hData->Sumw2();
  hMC->Sumw2();
  
  ntData->Project("hData","(DsvpvDistance/DsvpvDisErr)",TCut(seldata));
  ntMCprompt->Project("hMC","(DsvpvDistance/DsvpvDisErr)",TCut(selmc));


  hData->Scale(1./hData->Integral(1,20));
  hMC->Scale(1./hMC->Integral(1,20));
  
  cout<<hData->GetEntries()<<endl;
  cout<<hMC->GetEntries()<<endl;
  
  hData->Draw("e");
  hMC->Draw("esame");
  
  TH2F* hempty=new TH2F("hempty","",nBin,BinL-1,BinH+1,10,0.02,.5);  
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetYaxis()->SetTitle("Entries");
  hempty->GetXaxis()->SetTitle(title.Data());
  hempty->GetXaxis()->SetTitleOffset(0.9);
  hempty->GetYaxis()->SetTitleOffset(1.2);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.035);
  hempty->GetYaxis()->SetLabelSize(0.035);  
  
  TCanvas *canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  gPad->SetLogy();
  hempty->Draw();
  hData->SetLineWidth(2);
  hData->SetLineColor(1);
  hData->SetMarkerColor(1);
  hData->Draw("same");
  hMC->SetLineWidth(2);
  hMC->SetLineColor(2);
  hMC->SetMarkerColor(2);
  hMC->Draw("same");
  TLegend *legend=new TLegend(0.4958166,0.7558707,0.7949297,0.9299148,"");
  legend->SetBorderSize(0);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
  legend->SetTextFont(42);
  legend->SetTextSize(0.045);
  TLegendEntry *ent_Data=legend->AddEntry(hData,"Data","pf");
  ent_Data->SetTextFont(42);
  ent_Data->SetLineColor(1);
  TLegendEntry *ent_MC=legend->AddEntry(hMC,"MC","pf");
  ent_MC->SetTextFont(42);
  ent_MC->SetLineColor(2);
  legend->Draw("same");
  canvas->SaveAs(Form("canvas%s.pdf",titleplot.Data()));

}

int main(int argc, char *argv[])
{
  if((argc !=4))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 4)
    dataMC(argv[1],argv[2],argv[3]);
  return 0;
}


/*
.x dataMC.C+("d0/d0Err","d0/#sigma(d0)","d0D0err",40,0,200,1)
.x dataMC.C+("chi2cl","Vertex #chi^{2} Probability","ProbChi2",20,0,1)
.x dataMC.C+("cos(dtheta)","cos(#Delta#theta)","cosdtheta",20,-1,1)
.x dataMC.C+("abs(trk1Dxy/trk1D0Err)","|trk1Dxy/trk1D0Err|","trk1DxyDerr",40,0,200,1)
*/
