#include "../Bmeson/Bplusdsigmadpt_1ptbins.h"
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
    TString input = "D0_PbPb_acc_eff_ptbin_14_ybin_6_prompt_FONLLweight_cent-0to100_dataptshape_y1_Ncollweight1.root";
    TFile *inf = new TFile(input.Data());
  	TH1D* d0accxeff_pt = (TH1D*)inf->Get("d0accxeff_pt");
	//for(int i = 1; i <= d0accxeff_pt->GetNbinsX(); i++){  
	//	cout<<d0accxeff_pt->GetBinLowEdge(i)<<endl;
	//	cout<<d0accxeff_pt->GetBinLowEdge(i)+d0accxeff_pt->GetBinWidth(i)<<endl;;
	//	cout<<"---"<<endl;
	//}
	//return;
	string fonlldata="../fonllData/D_pp_pt_rap1_2p75_pt0to200.dat";
	double diffXsec = 0;
    double BRchain=3.93e-2+1.399e-4;
    double Fraction=0.557;
	
	const int BINS = 14;
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
    double bins[BINS+1] = {1.5, 2.5, 3.5, 4.5, 5.5, 7, 9, 11, 13, 16, 20, 28, 40, 60, 100};
	double effD[BINS] = {0.0004110786, 0.004868644, 0.01467014, 0.03677896, 0.07008137, 0.1185297, 0.1670131, 0.2038537, 0.2419492, 0.2935908, 0.3535325, 0.4175788, 0.4175788, 0.4175788};
	double effDerr[BINS] = {8.195642e-05, 0.0004130658 , 0.0009189658, 0.001791014, 0.002483458, 0.003386794, 0.004592759, 0.005642633, 0.005676826, 0.006227245, 0.00599039, 0.007349736, 0.007349736, 0.007349736};
	double RAA[BINS] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
	double PbPb2011[BINS-2] = {1628.9, 3459.4, 3915.6, 3172.1, 3383.4, 2191.3, 1010.7, 621.3, 436.0, 281.1, 158.6, 51.9};
    Double_t gaeSigmaDzero_fy3001[12] = {
    1.145624e+08,
    5.264917e+07,
    2.266795e+07,
    1.030325e+07,
    4395853,
    1554302,
    550404.4,
    226614.8,
    89997.46,
    29814.96,
    7003.541,
    1029.863};
	//double pp_datadriven_center_12ptbin[12]={ 191.417,85.7596, 34.7843, 15.1563, 6.07516, 1.95426, 0.631821, 0.248857, 0.100499, 0.0, 0.0, 0.0};

	for(int i = 0; i < BINS; i++){
		Bplusdsigmadpt_1ptbins(false, diffXsec, bins[i], bins[i+1], fonlldata);
		cout<<"diffXsec: "<<diffXsec * Fraction<<endl;
		apt[i] = (bins[i]+bins[i+1])/2;
		asigma[i] = diffXsec * Fraction;
		aptl[i] = (bins[i+1]-bins[i])/2;
		aerrorl[i] = 0;
		aerrorh[i] = 0;

		genB[i] = RAA[i] * 2 * 4.429e-6 * 208*208 * diffXsec * Fraction * (bins[i+1]-bins[i]);//4.429 ub-1 = 30M MB event used (HLT prescale ~ 37) 4.429/165.285*1.16e9 = 31.1M
		aerrorl2[i] = 0;
		aerrorh2[i] = 0;
		cout<<"genB: "<<genB[i]<<endl;
		
        recoB[i] = genB[i]*effD[i]*BRchain;	
		aerrorl3[i] = recoB[i]*effDerr[i];
		aerrorh3[i] = recoB[i]*effDerr[i];
		cout<<"aerrorl3: "<<aerrorl3[i]<<endl;
		cout<<"recoB: "<<recoB[i]<<endl;
	}
    TCanvas* cr = new TCanvas("cr","cr",600,500);
    gStyle->SetOptTitle(1);
    cr->SetLogy();

    TH2F* hempty=new TH2F("hempty","",1,0,100,10.,10,5e9);
	hempty->SetTitle("average B meson dsigma/dpt");
    hempty->Draw();
    TGraphAsymmErrors* dXsec = new TGraphAsymmErrors(BINS, apt, asigma, aptl, aptl, aerrorl, aerrorh);
	dXsec->SetLineWidth(3);
	dXsec->SetLineColor(4);
	dXsec->Draw("p same");
    TGraphAsymmErrors* JiandXsec = new TGraphAsymmErrors(BINS-2, apt, gaeSigmaDzero_fy3001, aptl, aptl, aerrorl, aerrorh);
	JiandXsec->SetLineWidth(3);
	JiandXsec->SetLineColor(2);
	JiandXsec->Draw("p same");
    TLegend *dXecleg = myLegend(0.50,0.7,0.86,0.89);
    dXecleg->AddEntry(dXsec,"FONLL dXsec","l");
    dXecleg->AddEntry(JiandXsec,"HIN-15-005","l");
    dXecleg->Draw();
	cr->SaveAs("Plots/dXsecD.png");

    TH2F* hempty2=new TH2F("hempty2","",1,0,100,10.,1e5,1e11);
	hempty2->SetTitle("# of produced D0 meson");
	hempty2->GetXaxis()->SetTitle("D0 pt");
    hempty2->Draw();
    TGraphAsymmErrors* GenB = new TGraphAsymmErrors(BINS, apt, genB, aptl, aptl, aerrorl2, aerrorh2);
	GenB->SetLineWidth(3);
	GenB->SetLineColor(4);
	GenB->Draw("p same");
	cr->SaveAs("Plots/GenD.png");

    TH2F* hempty3=new TH2F("hempty3","",1,0,100,10.,1e-1,1e5);
	hempty3->SetTitle("# of reconstructed D0");
	hempty3->GetXaxis()->SetTitle("D0 pt");
    hempty3->Draw();
    TGraphAsymmErrors* RecoB = new TGraphAsymmErrors(BINS, apt, recoB, aptl, aptl, aerrorl3, aerrorh3);
	RecoB->SetLineWidth(3);
	RecoB->SetLineColor(4);
	RecoB->Draw("p same");
    TGraphAsymmErrors* D2011 = new TGraphAsymmErrors(BINS-2, apt, PbPb2011, aptl, aptl, aerrorl3, aerrorh3);
    D2011->SetLineWidth(3);
    D2011->SetLineColor(2);
    D2011->Draw("p same");
    TLegend *leg = myLegend(0.50,0.7,0.86,0.89);
    leg->AddEntry(RecoB,"FONLL calculated","l");
    leg->AddEntry(D2011,"HIN-15-005 measured","l");
    leg->Draw();
	cr->SaveAs("Plots/RecoD.png");
}
