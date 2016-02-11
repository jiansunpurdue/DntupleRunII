#include <TH1D.h>
#include <TTree.h>
#include <TLegend.h>
#include <TCut.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>


// Take a tree, a variable and calculate the efficiency
TGraphAsymmErrors* getEfficiency(TTree *t, char *variable, TCut preselection, TCut cut, int nBin, Float_t *bins)
{
   static int count = 0;
   count++;
   TH1D *hPass = new TH1D (Form("hPass%d",count),"",nBin,bins);
   TH1D *hAll = new TH1D (Form("hAll%d",count),"",nBin,bins);
   t->Draw(Form("%s>>hAll%d",variable,count),preselection);
   t->Draw(Form("%s>>hPass%d",variable,count),preselection&&cut);

   TGraphAsymmErrors *g = new TGraphAsymmErrors;
   g->BayesDivide(hPass,hAll);
   return g;
}

void plotTrigger_PbPb(char *infname="/data/jisun/PbPb2015/skim_Dntuple_crab_PbPb_HIMinimumBias2_ForestAOD_Track_AK4CaloJet_D0_tkpt0p9eta1p5_01142016_tkpt8p0fortrigstudy.root")
{

   // ============== Open file and basic settings ===============   
   // Open Dntuple file
   TFile *inf = new TFile(infname);

   TTree *ntDkpi = (TTree*)inf->Get("ntDkpi");
   TTree *ntHlt = (TTree*)inf->Get("ntHlt");
   TTree *ntSkim = (TTree*)inf->Get("ntSkim");
   
   ntDkpi->AddFriend(ntHlt);
   ntDkpi->AddFriend(ntSkim);   

   // Define bin size and bin width for trigger turnon curve histograms
   const int nBin = 15;
   Float_t bins[nBin+1]={0,5,10,15,18,20,25,30,35,40,45,55,60,65,70,80};
 
   // Templates for plotting  
   TH1D *hTmp = new TH1D ("hTmp","",nBin,bins);
   TH1D *hTmp2 = new TH1D ("hTmp2","",nBin,bins);
   
   // ============== Selection Criteria ===============

   //TCut mbCut = "(HLT_HIL1MinimumBiasHF1AND_v1||  \
                  HLT_HIL1MinimumBiasHF2AND_v1||  \
                  HLT_HIL1MinimumBiasHF2AND_part1_v1)";
   TCut mbCut = "(HLT_HIL1MinimumBiasHF1AND_v1)";
   
   // L1 trigger thresholds
   TCut l1CutMBHF1And = "L1_MinimumBiasHF1_AND==1";
   TCut l1Cut28 = "L1_SingleS1Jet28_BptxAND==1";
   TCut l1Cut44 = "L1_SingleJet44_BptxAND==1";

   // D meson selection
   TCut DmassCut             = "(abs(Dmass-1.8696)<0.03)";
   TCut DmesonCut            = "(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12";
   //TCut DmesonDaughterTrkCut = "Dtrk1Pt>8.5&&abs(Dtrk1Eta)<2.0&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=10.5&&Dtrk1PtErr/Dtrk1Pt<0.1&& \
                                Dtrk2Pt>8.5&&abs(Dtrk2Eta)<2.0&&Dtrk2Algo<8&&(Dtrk2PixelHit+Dtrk2StripHit)>=10.5&&Dtrk2PtErr/Dtrk2Pt<0.1";
   //TCut DmesonDaughterTrkCut = "Dtrk1Pt>8.5&&abs(Dtrk1Eta)<2.0&&(Dtrk1Algo<8||Dtrk1Algo==11)&&(Dtrk1PixelHit+Dtrk1StripHit)>=10.5&&Dtrk1PtErr/Dtrk1Pt<0.1&& \
                                Dtrk2Pt>8.5&&abs(Dtrk2Eta)<2.0&&(Dtrk2Algo<8||Dtrk2Algo==11)&&(Dtrk2PixelHit+Dtrk2StripHit)>=10.5&&Dtrk2PtErr/Dtrk2Pt<0.1";
   TCut DmesonDaughterTrkCut = "Dtrk1Pt>8.5&&abs(Dtrk1Eta)<2.0&&(Dtrk1PixelHit+Dtrk1StripHit)>=10.5&&Dtrk1PtErr/Dtrk1Pt<0.1&& \
                                Dtrk2Pt>8.5&&abs(Dtrk2Eta)<2.0&&(Dtrk2PixelHit+Dtrk2StripHit)>=10.5&&Dtrk2PtErr/Dtrk2Pt<0.1";

   // Final selection for D candidates for trigger turnon studies
   TCut DAnaCut = DmassCut && DmesonCut && DmesonDaughterTrkCut;

   // HLT trigger thresholds
   TCut HLTCut20  = "HLT_HIDmesonHITrackingGlobal_Dpt20_v1";
   TCut HLTCut40 = "HLT_HIDmesonHITrackingGlobal_Dpt40_v1";
   TCut HLTCut60 = "HLT_HIDmesonHITrackingGlobal_Dpt60_v1";

   // ============== L1 trigger efficiency study ===============
   TCanvas *c = new TCanvas("c","",600,600);
   
   TGraphAsymmErrors* g20  = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1CutMBHF1And), HLTCut20, nBin, bins);
   TGraphAsymmErrors* g40 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1Cut28), HLTCut40, nBin, bins);
   TGraphAsymmErrors* g60 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1Cut44), HLTCut60, nBin, bins);

   hTmp->Draw();
   hTmp->SetXTitle("D Meson p_{T} (GeV/c)");
   hTmp->SetYTitle("D Meson HLT Efficiency");
   g20->SetLineColor(1);
   g20->SetMarkerColor(1);
   g20->Draw("pl same");
   
   g40->SetLineColor(2);
   g40->SetMarkerColor(2);
   g40->Draw("pl same");

   g60->SetLineColor(kGreen+2);
   g60->SetMarkerColor(kGreen+2);
   g60->Draw("pl same");

   TLegend *leg = new TLegend(0.53,0.2,0.93,0.6);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(g20,"PbPb #sqrt{s} = 5.02 TeV","");
   leg->AddEntry(g20,"HLT D meson 20","pl");
   leg->AddEntry(g40,"HLT D meson 40","pl");
   leg->AddEntry(g60,"HLT D meson 60","pl");
   leg->Draw();

   c->SaveAs("result/Dmeson-HIHLTriggerEfficiency.pdf");
   c->SaveAs("result/Dmeson-HIHLTriggerEfficiency.png");
   c->SaveAs("result/Dmeson-HIHLTriggerEfficiency.C");

   
   // ============== L1 trigger efficiency study ===============
   TCanvas *c2 = new TCanvas("c2","",600,600);
   
   TGraphAsymmErrors* gLMBHF1And = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_MinimumBiasHF1_AND_Prescl==1"), l1CutMBHF1And, nBin, bins);
   TGraphAsymmErrors* gL28 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_SingleS1Jet28_BptxAND_Prescl==1"), l1Cut28, nBin, bins);
   TGraphAsymmErrors* gL44 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_SingleJet44_BptxAND_Prescl==1"), l1Cut44, nBin, bins);
   
   hTmp2->Draw();
   hTmp2->SetXTitle("D meson p_{T} (GeV/c)");
   hTmp2->SetYTitle("L1 Trigger Efficiency");
   
   gLMBHF1And->SetMarkerColor(1);
   gLMBHF1And->SetLineColor(1);
   gLMBHF1And->Draw("pl same");

   gL28->SetMarkerColor(4);
   gL28->SetLineColor(4);
   gL28->Draw("pl same");

   gL44->SetMarkerColor(kGreen+2);
   gL44->SetLineColor(kGreen+2);
   gL44->Draw("pl same");

   TLegend *leg2 = new TLegend(0.53,0.2,0.93,0.6);
   leg2->SetBorderSize(0);
   leg2->SetFillStyle(0);
   leg2->AddEntry(gLMBHF1And,"PbPb #sqrt{s} = 5.02 TeV","");
   leg2->AddEntry(gLMBHF1And,"Level 1 Jet 16","pl");
   leg2->AddEntry(gL28,"Level 1 Jet 28","pl");
   leg2->AddEntry(gL44,"Level 1 Jet 44","pl");
   leg2->Draw();

   c2->SaveAs("result/Dmeson-HIL1TriggerEfficiency.pdf");
   c2->SaveAs("result/Dmeson-HIL1TriggerEfficiency.png");
   c2->SaveAs("result/Dmeson-HIL1TriggerEfficiency.C");
   
   
   // ============== Plot an example D mass distribution ===============
   TCanvas *c3 = new TCanvas("c3","",600,600);
   ntDkpi->Draw("Dmass>>h(100,1.7696,1.9696)",DmesonCut&&DmesonDaughterTrkCut&&mbCut&&l1CutMBHF1And);
   
   // ..done 
}

