#include "Bplusdsigmadpt_1ptbins.h"
#include "iostream"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include "utilities.h"
using namespace std;
void PlotXsec(){
    string fonlldata="../fonllData/B_pp_pt_rap24_2p75_pt0to200.dat";
    TString inputmc;
    TString cut;
    TString seldata_2y;
    TString selmc;
    TString selmcgen;
    inputmc="/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Bmeson_PbPbAna/Samples/BntupleMC_BfinderMC_Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_step3_BuKp_20150526_100kevt_20150604_20150608.root";
//    inputmc="/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Bmeson_PbPbAna/Samples/BntupleMC_BfinderMC_Pyquen_CMSSW742_Unquenched_PbPb_2760GeV_GEN_SIM_PU_BuKp_20150526_100kevt_20150505_20150609.root";

    cut="((abs(mu1eta)<1.0 && mu1pt>3.4) || (abs(mu1eta)>1.0 && abs(mu1eta)<1.5 && mu1pt>(5.8-2.4*abs(mu1eta))) || (abs(mu1eta)>1.5 && abs(mu1eta)<2.4 && mu1pt>(3.4-0.78*abs(mu1eta))))&&((abs(mu2eta)<1.0 && mu2pt>3.4) || (abs(mu2eta)>1.0 && abs(mu2eta)<1.5 && mu2pt>(5.8-2.4*abs(mu2eta))) || (abs(mu2eta)>1.5 && abs(mu2eta)<2.4 && mu2pt>(3.4-0.78*abs(mu2eta))))&&mu1TMOneStationTight && mu2TMOneStationTight && mu1isGlobalMuon && mu2isGlobalMuon && mu1trackerhit>10 && mu2trackerhit>10 && mu1InPixelLayer>1 && mu2InPixelLayer>1 && mu1normchi2<1.8 && mu2normchi2<1.8 && mu1dxyPV<3 && mu2dxyPV<3 && mu1dzPV<30 && mu2dzPV<30 && mumupt>3 && trk1Pt>1 && pt>10 && pt < 60 && mass>5&&mass<6 && abs(y)<2.4 && abs(mumumass-3.096916)<0.15 && isbestchi2 && chi2cl > 0.004345 && cos(dtheta) > 0.256418 && d0/d0Err > 6.679356";//pt10, acc muon, hp tk, SA
//    cut="((abs(mu1eta)<1.0 && mu1pt>3.4) || (abs(mu1eta)>1.0 && abs(mu1eta)<1.5 && mu1pt>(5.8-2.4*abs(mu1eta))) || (abs(mu1eta)>1.5 && abs(mu1eta)<2.4 && mu1pt>(3.4-0.78*abs(mu1eta))))&&((abs(mu2eta)<1.0 && mu2pt>3.4) || (abs(mu2eta)>1.0 && abs(mu2eta)<1.5 && mu2pt>(5.8-2.4*abs(mu2eta))) || (abs(mu2eta)>1.5 && abs(mu2eta)<2.4 && mu2pt>(3.4-0.78*abs(mu2eta))))&&mu1TMOneStationTight && mu2TMOneStationTight && mu1isGlobalMuon && mu2isGlobalMuon && mu1trackerhit>10 && mu2trackerhit>10 && mu1InPixelLayer>1 && mu2InPixelLayer>1 && mu1normchi2<1.8 && mu2normchi2<1.8 && mu1dxyPV<3 && mu2dxyPV<3 && mu1dzPV<30 && mu2dzPV<30 && mumupt>3 && trk1Pt>1 && mass>5&&mass<6 && abs(y)<2.4 && abs(mumumass-3.096916)<0.15 && isbestchi2 && chi2cl > 0.004345 && cos(dtheta) > 0.256418 && d0/d0Err > 6.679356";//
//    cut="((abs(mu1eta)<1.0 && mu1pt>3.4) || (abs(mu1eta)>1.0 && abs(mu1eta)<1.5 && mu1pt>(5.8-2.4*abs(mu1eta))) || (abs(mu1eta)>1.5 && abs(mu1eta)<2.4 && mu1pt>(3.4-0.78*abs(mu1eta))))&&((abs(mu2eta)<1.0 && mu2pt>3.4) || (abs(mu2eta)>1.0 && abs(mu2eta)<1.5 && mu2pt>(5.8-2.4*abs(mu2eta))) || (abs(mu2eta)>1.5 && abs(mu2eta)<2.4 && mu2pt>(3.4-0.78*abs(mu2eta)))) && mu1TMOneStationTight && mu2TMOneStationTight && mumupt>3 && trk1Pt>1 && abs(y)<2.4 && trk1Chi2ndf < 5";
//    cut="mu1TMOneStationTight && mu2TMOneStationTight && mumupt>3 && trk1Pt>1 && abs(y)<2.4 && trk1Chi2ndf < 5";
    selmc=Form("abs(y)<2.4&&gen==23333&&%s",cut.Data());
    selmcgen="abs(y)<2.4&&abs(pdgId)==521&&isSignal==1";
    TFile *infMC = new TFile(inputmc.Data());
    TTree *ntMC = (TTree*)infMC->Get("ntKp");
    TTree *ntGen = (TTree*)infMC->Get("ntGen");
    TH1D *hPtMC = new TH1D("hPtMC","",1,5,100);
    TH1D *hPtGen = new TH1D("hPtGen","",1,5,100);
    
	double diffXsec = 0;
    double BRchain=6.09604e-5;
    double Fraction=0.401;
	
	const int BINS = 1;
	double apt[BINS];
	double asigma[BINS];
	double aptl[BINS];
	double aerrorl[BINS];
	double aerrorh[BINS];
	double aerrorl2[BINS];
	double aerrorh2[BINS];
	double aerrorl3[BINS];
	double aerrorh3[BINS];
	double genB[BINS];
	double effB[BINS];
	double recoB[BINS];
	double bins[BINS+1] = {10,60};
	double PbPb2011[BINS] = {33.1271};
	for(int i = 0; i < BINS; i++){
		Bplusdsigmadpt_1ptbins(false, diffXsec, bins[i], bins[i+1], fonlldata);
		cout<<"diffXsec: "<<diffXsec<<endl;
		apt[i] = (bins[i]+bins[i+1])/2;
		asigma[i] = diffXsec;
		aptl[i] = (bins[i+1]-bins[i])/2;
		aerrorl[i] = 0;
		aerrorh[i] = 0;

		genB[i] = 0.3 * 2 * 166e-6 * 208*208 * diffXsec * (bins[i+1]-bins[i]);
		aerrorl2[i] = 0;
		aerrorh2[i] = 0;
		cout<<"genB: "<<genB[i]<<endl;
		
    	hPtMC = new TH1D(Form("hPtMC%d", i),"",1,bins[i],bins[i+1]);
		hPtGen = new TH1D(Form("hPtGen%d", i),"",1,bins[i],bins[i+1]);
    	ntMC->Project(Form("hPtMC%d", i),"pt",(TCut(selmc.Data())&&"gen==23333"));
	    ntGen->Project(Form("hPtGen%d", i),"pt",(TCut(selmcgen.Data())));
		effB[i] = hPtMC->GetBinContent(1)/hPtGen->GetBinContent(1);
		cout<<hPtMC->GetBinContent(1)<<endl;
		cout<<hPtGen->GetBinContent(1)<<endl;
		cout<<"eff: "<<effB[i]<<endl;
		//if(hPtGen->GetBinContent(1)<500) effB[i] = effB[i-1];//too few statistic 
		recoB[i] = genB[i]*effB[i]*BRchain*Fraction;
		aerrorl3[i] = recoB[i]*sqrt(effB[i]*(1-effB[i])/hPtGen->GetBinContent(1));//binomial error
		aerrorh3[i] = recoB[i]*sqrt(effB[i]*(1-effB[i])/hPtGen->GetBinContent(1));
		cout<<"aerrorl3: "<<aerrorl3[i]<<endl;
		cout<<"recoB: "<<recoB[i]<<endl;
		cout<<"2011 B PbPb: 33.1271"<<endl;
	}
    TCanvas* cr = new TCanvas("cr","cr",600,500);
    gStyle->SetOptTitle(1);
    cr->SetLogy();

    TH2F* hempty=new TH2F("hempty","",1,0,100,10.,10,10000000);
	hempty->SetTitle("average B meson dsigma/dpt");
    hempty->Draw();
    TGraphAsymmErrors* dXsec = new TGraphAsymmErrors(BINS, apt, asigma, aptl, aptl, aerrorl, aerrorh);
	dXsec->SetLineWidth(3);
	dXsec->SetLineColor(4);
	dXsec->Draw("p same");

    TH2F* hempty2=new TH2F("hempty2","",1,0,100,10.,1e2,1e11);
	hempty2->SetTitle("# of produced B meson");
	hempty2->GetXaxis()->SetTitle("B pt");
    hempty2->Draw();
    TGraphAsymmErrors* GenB = new TGraphAsymmErrors(BINS, apt, genB, aptl, aptl, aerrorl2, aerrorh2);
	GenB->SetLineWidth(3);
	GenB->SetLineColor(4);
	GenB->Draw("p same");
	cr->SaveAs("Plots/GenB.png");

    TH2F* hempty3=new TH2F("hempty3","",1,0,100,10.,1e1,1e2);
	hempty3->SetTitle("# of reconstructed B+");
	hempty3->GetXaxis()->SetTitle("B pt");
    hempty3->Draw();
    TGraphAsymmErrors* RecoB = new TGraphAsymmErrors(BINS, apt, recoB, aptl, aptl, aerrorl3, aerrorh3);
	RecoB->SetLineWidth(3);
	RecoB->SetLineColor(4);
	RecoB->Draw("p same");
    TGraphAsymmErrors* B2011 = new TGraphAsymmErrors(BINS, apt, PbPb2011, aptl, aptl, aerrorl3, aerrorh3);
	B2011->SetLineWidth(3);
	B2011->SetLineColor(2);
	B2011->Draw("p same");
    TLegend *leg = myLegend(0.50,0.7,0.86,0.89);
    leg->AddEntry(RecoB,"FONLL calculated","l");
    leg->AddEntry(B2011,"2011 data mesured","l");
	leg->Draw();
	cr->SaveAs("Plots/RecoB.png");
}
