void c1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Oct 22 22:06:04 2015) by ROOT version6.02/05
   TCanvas *c1 = new TCanvas("c1", "c1",0,22,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1->Range(-4.603658,-0.07073505,42.34756,0.4733808);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.13);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   Double_t xAxis1[15] = {0, 1.5, 2.5, 3.5, 4.5, 5.5, 7, 9, 11, 13, 16, 20, 28, 40, 100}; 
   
   TH1D *d0accxeff_pt1 = new TH1D("d0accxeff_pt1","d0candpt_matched_cuts",14, xAxis1);
   d0accxeff_pt1->SetBinContent(2,0.0004110786);
   d0accxeff_pt1->SetBinContent(3,0.004868644);
   d0accxeff_pt1->SetBinContent(4,0.01467014);
   d0accxeff_pt1->SetBinContent(5,0.03677896);
   d0accxeff_pt1->SetBinContent(6,0.07008137);
   d0accxeff_pt1->SetBinContent(7,0.1185297);
   d0accxeff_pt1->SetBinContent(8,0.1670131);
   d0accxeff_pt1->SetBinContent(9,0.2038537);
   d0accxeff_pt1->SetBinContent(10,0.2419492);
   d0accxeff_pt1->SetBinContent(11,0.2935908);
   d0accxeff_pt1->SetBinContent(12,0.3535325);
   d0accxeff_pt1->SetBinContent(13,0.4175788);
   d0accxeff_pt1->SetBinError(2,8.195642e-05);
   d0accxeff_pt1->SetBinError(3,0.0004130658);
   d0accxeff_pt1->SetBinError(4,0.0009189658);
   d0accxeff_pt1->SetBinError(5,0.001791014);
   d0accxeff_pt1->SetBinError(6,0.002483458);
   d0accxeff_pt1->SetBinError(7,0.003386794);
   d0accxeff_pt1->SetBinError(8,0.004592759);
   d0accxeff_pt1->SetBinError(9,0.005642633);
   d0accxeff_pt1->SetBinError(10,0.005676826);
   d0accxeff_pt1->SetBinError(11,0.006227245);
   d0accxeff_pt1->SetBinError(12,0.00599039);
   d0accxeff_pt1->SetBinError(13,0.007349736);
   d0accxeff_pt1->SetEntries(442985);
   d0accxeff_pt1->SetLineWidth(2);
   d0accxeff_pt1->GetXaxis()->SetTitle("D0 p_{T} (GeV/c)");
   d0accxeff_pt1->GetXaxis()->SetRange(2,13);
   d0accxeff_pt1->GetXaxis()->SetTitleOffset(1.25);
   d0accxeff_pt1->GetYaxis()->SetTitle("#alpha #times #varepsilon_{reco} #times #varepsilon_{cuts}");
   d0accxeff_pt1->GetYaxis()->SetTitleOffset(1.8);
   d0accxeff_pt1->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.9341608,0.6004698,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(10);
   TText *AText = pt->AddText("d0candpt_matched_cuts");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
