void CrossSection(TString particle="Dzero"){

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TFile*filePPReference=new TFile("../../FONLL/Results/output_pp_d0meson5_5TeV_y1.root");  
  TGraphAsymmErrors*gaeBplusReference=(TGraphAsymmErrors*)filePPReference->Get(Form("gaeSigma%s",particle.Data()));
  
  TFile*filepPb=new TFile(Form("ResultsD0_pp/PtSigma%s.root",particle.Data()));
  TH1F*hSigmapPbStat=(TH1F*)filepPb->Get("hPtSigma");    
  double scalingfactor=1;
    
  for (int i=0;i<hSigmapPbStat->GetNbinsX();i++){
    hSigmapPbStat->SetBinContent(i+1,scalingfactor*hSigmapPbStat->GetBinContent(i+1));
    hSigmapPbStat->SetBinError(i+1,scalingfactor*hSigmapPbStat->GetBinError(i+1));
  }
  
  TCanvas *canvasSigma=new TCanvas("canvasSigma","canvasSigma",500,500);   
  canvasSigma->cd();
  canvasSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
  canvasSigma->SetFillColor(0);
  canvasSigma->SetBorderMode(0);
  canvasSigma->SetBorderSize(2);
  canvasSigma->SetLeftMargin(0.1451613);
  canvasSigma->SetRightMargin(0.05443548);
  canvasSigma->SetTopMargin(0.08474576);
  canvasSigma->SetBottomMargin(0.1165254);
  canvasSigma->SetFrameBorderMode(0);
  canvasSigma->SetFrameBorderMode(0);
  canvasSigma->SetLogy();
  
  TH2F* hemptySigma=new TH2F("hemptySigma","",50,50.,110,10.,1,1e3);  
  hemptySigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptySigma->GetXaxis()->CenterTitle();
  hemptySigma->GetYaxis()->CenterTitle();
  hemptySigma->GetYaxis()->SetTitle("d#sigma / dp_{T}( pb GeV^{-1}c)");
  

  hemptySigma->GetXaxis()->SetTitleOffset(1.);
  hemptySigma->GetYaxis()->SetTitleOffset(1.3);
  hemptySigma->GetXaxis()->SetTitleSize(0.045);
  hemptySigma->GetYaxis()->SetTitleSize(0.045);
  hemptySigma->GetXaxis()->SetTitleFont(42);
  hemptySigma->GetYaxis()->SetTitleFont(42);
  hemptySigma->GetXaxis()->SetLabelFont(42);
  hemptySigma->GetYaxis()->SetLabelFont(42);
  hemptySigma->GetXaxis()->SetLabelSize(0.04);
  hemptySigma->GetYaxis()->SetLabelSize(0.04);  
  hemptySigma->SetMaximum(2);
  hemptySigma->SetMinimum(0.);
  hemptySigma->Draw();
  hSigmapPbStat->Draw("epsame");  
  gaeBplusReference->Draw("2same");  
  
  TLegend *legendSigma=new TLegend(0.5100806,0.5868644,0.8084677,0.7605932,"");
  legendSigma->SetBorderSize(0);
  legendSigma->SetLineColor(0);
  legendSigma->SetFillColor(0);
  legendSigma->SetFillStyle(1001);
  legendSigma->SetTextFont(42);
  legendSigma->SetTextSize(0.045);
  
  TLegendEntry *ent_SigmapPb=legendSigma->AddEntry(hSigmapPbStat,"pp","pf");
  ent_SigmapPb->SetTextFont(42);
  ent_SigmapPb->SetLineColor(1);
  ent_SigmapPb->SetMarkerColor(1);
  
  TLegendEntry *ent_Sigmapp=legendSigma->AddEntry(gaeBplusReference,"FONLL pp ref.","f");
  ent_Sigmapp->SetTextFont(42);
  ent_Sigmapp->SetLineColor(5);
  ent_Sigmapp->SetMarkerColor(1);
  legendSigma->Draw("same");
    
  TLatex * tlatex1=new TLatex(0.1612903,0.8625793,"CMS Preliminary     pp #sqrt{s}= 5.02 TeV");
  tlatex1->SetNDC();
  tlatex1->SetTextColor(1);
  tlatex1->SetTextFont(42);
  tlatex1->SetTextSize(0.045);
  tlatex1->Draw();
    
  TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 9.97 pb^{-1}");
  tlatexlumi->SetNDC();
  tlatexlumi->SetTextColor(1);
  tlatexlumi->SetTextFont(42);
  tlatexlumi->SetTextSize(0.045);
  tlatexlumi->Draw();

  canvasSigma->SaveAs("canvasSigmaDzero.pdf");


}

