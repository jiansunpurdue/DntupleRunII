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

void plotTrigger(char *infname="mb.root")
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
   const int nBin = 12;
   Float_t bins[nBin+1]={0,5,6,8,10,12,15,20,25,30,35,40,70};
 
   // Templates for plotting  
   TH1D *hTmp = new TH1D ("hTmp","",nBin,bins);
   TH1D *hTmp2 = new TH1D ("hTmp2","",nBin,bins);
   
   // ============== Selection Criteria ===============

   // This MB sample has part0;part5;part9; part13; part17 MB triggers
   TCut mbCut = "(HLT_L1MinimumBiasHF1OR_part0_v1||  \
                  HLT_L1MinimumBiasHF1OR_part5_v1||  \
                  HLT_L1MinimumBiasHF1OR_part9_v1||  \
                  HLT_L1MinimumBiasHF1OR_part13_v1|| \
                  HLT_L1MinimumBiasHF1OR_part17_v1)";
   
   // L1 trigger thresholds
   TCut l1Cut16 = "L1_SingleJet16_BptxAND==1";
   TCut l1Cut24 = "L1_SingleJet24_BptxAND==1";
   TCut l1Cut28 = "L1_SingleJet28_BptxAND==1";
   TCut l1Cut40 = "L1_SingleJet40_BptxAND==1";

   // D meson selection
   TCut DmassCut             = "(abs(Dmass-1.8696)<0.03)";
   TCut DmesonCut            = "(DsvpvDistance/DsvpvDisErr)>3.5&&Dchi2cl>0.05&&Dalpha<0.12";
   TCut DmesonDaughterTrkCut = "Dtrk1Pt>2.5&&abs(Dtrk1Eta)<2.0&&Dtrk1Algo<8&&(Dtrk1PixelHit+Dtrk1StripHit)>=11&& \
                                Dtrk2Pt>2.5&&abs(Dtrk2Eta)<2.0&&Dtrk2Algo<8&&(Dtrk2PixelHit+Dtrk2StripHit)>=11";

   // Final selection for D candidates for trigger turnon studies
   TCut DAnaCut = DmassCut && DmesonCut && DmesonDaughterTrkCut;

   // HLT trigger thresholds
   TCut HLTCut8  = "HLT_DmesonPPTrackingGlobal_Dpt8_v1";
   TCut HLTCut15 = "HLT_DmesonPPTrackingGlobal_Dpt15_v1";
   TCut HLTCut20 = "HLT_DmesonPPTrackingGlobal_Dpt20_v1";
   TCut HLTCut30 = "HLT_DmesonPPTrackingGlobal_Dpt30_v1";

   // ============== L1 trigger efficiency study ===============
   TCanvas *c = new TCanvas("c","",600,600);
   
   TGraphAsymmErrors* g8  = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1Cut16), HLTCut8, nBin, bins);
   TGraphAsymmErrors* g15 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1Cut24), HLTCut15, nBin, bins);
   TGraphAsymmErrors* g20 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1Cut28), HLTCut20, nBin, bins);
   TGraphAsymmErrors* g30 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&l1Cut40), HLTCut30, nBin, bins);

   hTmp->Draw();
   hTmp->SetXTitle("D Meson p_{T} (GeV/c)");
   hTmp->SetYTitle("D Meson HLT Efficiency");
   g8->SetLineColor(1);
   g8->SetMarkerColor(1);
   g8->Draw("pl same");
   
   g15->SetLineColor(2);
   g15->SetMarkerColor(2);
   g15->Draw("pl same");

   g20->SetLineColor(4);
   g20->SetMarkerColor(4);
   g20->Draw("pl same");
   
   g30->SetLineColor(kGreen+2);
   g30->SetMarkerColor(kGreen+2);
   g30->Draw("pl same");

   TLegend *leg = new TLegend(0.53,0.2,0.93,0.6);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(g8,"pp #sqrt{s} = 5.02 TeV","");
   leg->AddEntry(g8,"HLT D meson 8","pl");
   leg->AddEntry(g15,"HLT D meson 15","pl");
   leg->AddEntry(g20,"HLT D meson 20","pl");
   leg->AddEntry(g30,"HLT D meson 30","pl");
   leg->Draw();

   c->SaveAs("result/Dmeson-HLTriggerEfficiency.pdf");
   c->SaveAs("result/Dmeson-HLTriggerEfficiency.C");

   
   // ============== HLT trigger efficiency study ===============
   TCanvas *c2 = new TCanvas("c2","",600,600);
   
   TGraphAsymmErrors* gL16 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_SingleJet16_BptxAND_Prescl==1"), l1Cut16, nBin, bins);
   TGraphAsymmErrors* gL24 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_SingleJet24_BptxAND_Prescl==1"), l1Cut24, nBin, bins);
   TGraphAsymmErrors* gL28 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_SingleJet28_BptxAND_Prescl==1"), l1Cut28, nBin, bins);
   TGraphAsymmErrors* gL40 = getEfficiency(ntDkpi,Form("Max$(Dpt*(%s))",DAnaCut.GetTitle()), TCut(DAnaCut&&mbCut&&"L1_SingleJet40_BptxAND_Prescl==1"), l1Cut40, nBin, bins);
   
   hTmp2->Draw();
   hTmp2->SetXTitle("D meson p_{T} (GeV/c)");
   hTmp2->SetYTitle("L1 Trigger Efficiency");
   
   gL16->SetMarkerColor(1);
   gL16->SetLineColor(1);
   gL16->Draw("pl same");

   gL24->SetMarkerColor(2);
   gL24->SetLineColor(2);
   gL24->Draw("pl same");

   gL28->SetMarkerColor(4);
   gL28->SetLineColor(4);
   gL28->Draw("pl same");

   gL40->SetMarkerColor(kGreen+2);
   gL40->SetLineColor(kGreen+2);
   gL40->Draw("pl same");

   TLegend *leg2 = new TLegend(0.53,0.2,0.93,0.6);
   leg2->SetBorderSize(0);
   leg2->SetFillStyle(0);
   leg2->AddEntry(gL16,"pp #sqrt{s} = 5.02 TeV","");
   leg2->AddEntry(gL16,"Level 1 Jet 16","pl");
   leg2->AddEntry(gL24,"Level 1 Jet 24","pl");
   leg2->AddEntry(gL28,"Level 1 Jet 28","pl");
   leg2->AddEntry(gL40,"Level 1 Jet 40","pl");
   leg2->Draw();

   c2->SaveAs("result/Dmeson-L1TriggerEfficiency.pdf");
   c2->SaveAs("result/Dmeson-L1TriggerEfficiency.C");
   
   
   // ============== Plot an example D mass distribution ===============
   TCanvas *c3 = new TCanvas("c3","",600,600);
   ntDkpi->Draw("Dmass>>h(100,1.7696,1.9696)",DmesonCut&&DmesonDaughterTrkCut&&mbCut&&l1Cut16);
   
   // ..done 
}

