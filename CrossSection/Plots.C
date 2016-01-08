void Plots(){
  TFile *fBtoD=new TFile("BtoD-0-100.root");
  TH1F*hBtoD=(TH1F*)fBtoD->Get("hDFromBPt");
  TH1F*hDPt=(TH1F*)fBtoD->Get("hDPt");

  hBtoD->Rebin(46);
  hDPt->Rebin(46);
  hBtoD->Divide(hDPt);
  hBtoD->Draw();


}
