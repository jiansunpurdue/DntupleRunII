#include <TTree.h>

#define MAX_XB 16384
#define MAX_MUON 512
#define MAX_TRACK 4096 //default 2048
#define MAX_GEN 8192 //default 2048
#define MAX_BX 128
#define MAX_Vertices 4096
//#define N_TRIGGER_BOOKINGS 5842

//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Build Branch //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//EvtInfo
Int_t      RunNo;
Int_t      EvtNo;
Int_t      Bsize;
Int_t      Dsize;
Double_t   PVx;
Double_t   PVy;
Double_t   PVz;
Double_t   PVxE;
Double_t   PVyE;
Double_t   PVzE;
Double_t   PVnchi2;
Double_t   PVchi2;
//BInfo
Int_t      Bindex[MAX_XB];
Double_t   Bmass[MAX_XB];
Double_t   Bpt[MAX_XB];
Double_t   Beta[MAX_XB];
Double_t   Bphi[MAX_XB];
Double_t   By[MAX_XB];
Double_t   BvtxX[MAX_XB];
Double_t   BvtxY[MAX_XB];
Double_t   Bd0[MAX_XB];
Double_t   Bd0Err[MAX_XB];
Double_t   Bchi2ndf[MAX_XB];
Double_t   Bchi2cl[MAX_XB];
Double_t   Bdtheta[MAX_XB];
Double_t   Blxy[MAX_XB];
Bool_t     Bisbestchi2[MAX_XB];
//BInfo.MuInfo
Bool_t     Bmu1TrackerMuArbitrated[MAX_XB];
Bool_t     Bmu2TrackerMuArbitrated[MAX_XB];
Bool_t     Bmu1isTrackerMuon[MAX_XB];
Bool_t     Bmu2isTrackerMuon[MAX_XB];
Bool_t     Bmu1isGlobalMuon[MAX_XB];
Bool_t     Bmu2isGlobalMuon[MAX_XB];
Bool_t     Bmu1TMOneStationTight[MAX_XB];
Bool_t     Bmu2TMOneStationTight[MAX_XB];
Int_t      Bmu1InPixelLayer[MAX_XB];
Int_t      Bmu2InPixelLayer[MAX_XB];
Int_t      Bmu1InStripLayer[MAX_XB];
Int_t      Bmu2InStripLayer[MAX_XB];
Int_t      Bmu1InTrackerLayer[MAX_XB];
Int_t      Bmu2InTrackerLayer[MAX_XB];
Double_t   Bmu1eta[MAX_XB];
Double_t   Bmu2eta[MAX_XB];
Double_t   Bmu1phi[MAX_XB];
Double_t   Bmu2phi[MAX_XB];
Double_t   Bmu1y[MAX_XB];
Double_t   Bmu2y[MAX_XB];
Double_t   Bmu1pt[MAX_XB];
Double_t   Bmu2pt[MAX_XB];
Double_t   Bmu1p[MAX_XB];
Double_t   Bmu2p[MAX_XB];
Double_t   Bmu1E[MAX_XB];
Double_t   Bmu2E[MAX_XB];
Double_t   Bmu1dzPV[MAX_XB];
Double_t   Bmu2dzPV[MAX_XB];
Double_t   Bmu1dxyPV[MAX_XB];
Double_t   Bmu2dxyPV[MAX_XB];
Double_t   Bmu1normchi2[MAX_XB];
Double_t   Bmu2normchi2[MAX_XB];
Double_t   Bmu1Chi2ndf[MAX_XB];
Double_t   Bmu2Chi2ndf[MAX_XB];
//BInfo.ujInfo
Double_t   Bmumumass[MAX_XB];
Double_t   Bmumueta[MAX_XB];
Double_t   Bmumuphi[MAX_XB];
Double_t   Bmumuy[MAX_XB];
Double_t   Bmumupt[MAX_XB];
Double_t   Bujmass[MAX_XB];
Double_t   BujvProb[MAX_XB];
Double_t   Bujpt[MAX_XB];
Double_t   Bujeta[MAX_XB];
Double_t   Bujphi[MAX_XB];
Double_t   Bujy[MAX_XB];
Double_t   Bujlxy[MAX_XB];
//BInfo.trkInfo
Double_t   Btrk1Pt[MAX_XB];
Double_t   Btrk2Pt[MAX_XB];
Double_t   Btrk1Eta[MAX_XB];
Double_t   Btrk2Eta[MAX_XB];
Double_t   Btrk1Phi[MAX_XB];
Double_t   Btrk2Phi[MAX_XB];
Double_t   Btrk1Y[MAX_XB];
Double_t   Btrk2Y[MAX_XB];
Double_t   Btrk1Dxy[MAX_XB];
Double_t   Btrk2Dxy[MAX_XB];
Double_t   Btrk1D0Err[MAX_XB];
Double_t   Btrk2D0Err[MAX_XB];
Double_t   Btrk1PixelHit[MAX_XB];
Double_t   Btrk2PixelHit[MAX_XB];
Double_t   Btrk1StripHit[MAX_XB];
Double_t   Btrk2StripHit[MAX_XB];
Double_t   Btrk1Chi2ndf[MAX_XB];
Double_t   Btrk2Chi2ndf[MAX_XB];
//BInfo.tktkInfo
Double_t   Btktkmass[MAX_XB];
Double_t   BtktkmassKK[MAX_XB];
Double_t   BtktkvProb[MAX_XB];
Double_t   Btktkpt[MAX_XB];
Double_t   Btktketa[MAX_XB];
Double_t   Btktkphi[MAX_XB];
Double_t   Btktky[MAX_XB];
Double_t   Bdoubletmass[MAX_XB];
Double_t   Bdoubletpt[MAX_XB];
Double_t   Bdoubleteta[MAX_XB];
Double_t   Bdoubletphi[MAX_XB];
Double_t   Bdoublety[MAX_XB];
//BInfo.genInfo
Double_t   Bgen[MAX_XB];
Int_t      BgenIndex[MAX_XB];
Double_t   Bgenpt[MAX_XB];
Double_t   Bgeneta[MAX_XB];
Double_t   Bgenphi[MAX_XB];
Double_t   Bgeny[MAX_XB];
Int_t      BgenpdgId[MAX_XB];

void buildBBranch(TTree* bnt)
{
  //EvtInfo
  bnt->Branch("RunNo",&RunNo);
  bnt->Branch("EvtNo",&EvtNo);
  bnt->Branch("Bsize",&Bsize);
  bnt->Branch("PVx",&PVx);
  bnt->Branch("PVy",&PVy);
  bnt->Branch("PVz",&PVz);
  bnt->Branch("PVxE",&PVxE);
  bnt->Branch("PVyE",&PVyE);
  bnt->Branch("PVzE",&PVzE);
  bnt->Branch("PVnchi2",&PVnchi2);
  bnt->Branch("PVchi2",&PVchi2);
  //BInfo
  bnt->Branch("Bindex",Bindex,"Bindex[size]/I");
  bnt->Branch("Bmass",Bmass,"Bmass[size]/D");
  bnt->Branch("Bpt",Bpt,"Bpt[size]/D");
  bnt->Branch("Beta",Beta,"Beta[size]/D");
  bnt->Branch("Bphi",Bphi,"Bphi[size]/D");
  bnt->Branch("By",By,"By[size]/D");
  bnt->Branch("BvtxX",BvtxX,"BvtxX[size]/D");
  bnt->Branch("BvtxY",BvtxY,"BvtxY[size]/D");
  bnt->Branch("Bd0",Bd0,"Bd0[size]/D");
  bnt->Branch("Bd0Err",Bd0Err,"Bd0Err[size]/D");
  bnt->Branch("Bchi2ndf",Bchi2ndf,"Bchi2ndf[size]/D");
  bnt->Branch("Bchi2cl",Bchi2cl,"Bchi2cl[size]/D");
  bnt->Branch("Bdtheta",Bdtheta,"Bdtheta[size]/D");
  bnt->Branch("Blxy",Blxy,"Blxy[size]/D");
  bnt->Branch("Bisbestchi2",Bisbestchi2,"Bisbestchi2[size]/O");
  //BInfo.muInfo
  bnt->Branch("Bmu1TrackerMuArbitrated",Bmu1TrackerMuArbitrated,"Bmu1TrackerMuArbitrated[size]/O");
  bnt->Branch("Bmu2TrackerMuArbitrated",Bmu2TrackerMuArbitrated,"Bmu2TrackerMuArbitrated[size]/O");
  bnt->Branch("Bmu1isTrackerMuon",Bmu1isTrackerMuon,"Bmu1isTrackerMuon[size]/O");
  bnt->Branch("Bmu2isTrackerMuon",Bmu2isTrackerMuon,"Bmu2isTrackerMuon[size]/O");
  bnt->Branch("Bmu1isGlobalMuon",Bmu1isGlobalMuon,"Bmu1isGlobalMuon[size]/O");
  bnt->Branch("Bmu2isGlobalMuon",Bmu2isGlobalMuon,"Bmu2isGlobalMuon[size]/O");
  bnt->Branch("Bmu1TMOneStationTight",Bmu1TMOneStationTight,"Bmu1TMOneStationTight[size]/O");
  bnt->Branch("Bmu2TMOneStationTight",Bmu2TMOneStationTight,"Bmu2TMOneStationTight[size]/O");
  bnt->Branch("Bmu1InPixelLayer",Bmu1InPixelLayer,"Bmu1InPixelLayer[size]/I");
  bnt->Branch("Bmu2InPixelLayer",Bmu2InPixelLayer,"Bmu2InPixelLayer[size]/I");
  bnt->Branch("Bmu1InStripLayer",Bmu1InStripLayer,"Bmu1InStripLayer[size]/I");
  bnt->Branch("Bmu2InStripLayer",Bmu2InStripLayer,"Bmu2InStripLayer[size]/I");
  bnt->Branch("Bmu1InTrackerLayer",Bmu1InTrackerLayer,"Bmu1InTrackerLayer[size]/I");
  bnt->Branch("Bmu2InTrackerLayer",Bmu2InTrackerLayer,"Bmu2InTrackerLayer[size]/I");
  bnt->Branch("Bmu1eta",Bmu1eta,"Bmu1eta[size]/D");
  bnt->Branch("Bmu2eta",Bmu2eta,"Bmu2eta[size]/D");
  bnt->Branch("Bmu1phi",Bmu1phi,"Bmu1phi[size]/D");
  bnt->Branch("Bmu2phi",Bmu2phi,"Bmu2phi[size]/D");
  bnt->Branch("Bmu1y",Bmu1y,"Bmu1y[size]/D");
  bnt->Branch("Bmu2y",Bmu2y,"Bmu2y[size]/D");
  bnt->Branch("Bmu1pt",Bmu1pt,"Bmu1pt[size]/D");
  bnt->Branch("Bmu2pt",Bmu2pt,"Bmu2pt[size]/D");
  bnt->Branch("Bmu1p",Bmu1p,"Bmu1p[size]/D");
  bnt->Branch("Bmu2p",Bmu2p,"Bmu2p[size]/D");
  bnt->Branch("Bmu1E",Bmu1E,"Bmu1E[size]/D");
  bnt->Branch("Bmu2E",Bmu2E,"Bmu2E[size]/D");
  bnt->Branch("Bmu1dzPV",Bmu1dzPV,"Bmu1dzPV[size]/D");
  bnt->Branch("Bmu2dzPV",Bmu2dzPV,"Bmu2dzPV[size]/D");
  bnt->Branch("Bmu1dxyPV",Bmu1dxyPV,"Bmu1dxyPV[size]/D");
  bnt->Branch("Bmu2dxyPV",Bmu2dxyPV,"Bmu2dxyPV[size]/D");
  bnt->Branch("Bmu1normchi2",Bmu1normchi2,"Bmu1normchi2[size]/D");
  bnt->Branch("Bmu2normchi2",Bmu2normchi2,"Bmu2normchi2[size]/D");
  bnt->Branch("Bmu1Chi2ndf",Bmu1Chi2ndf,"Bmu1Chi2ndf[size]/D");
  bnt->Branch("Bmu2Chi2ndf",Bmu2Chi2ndf,"Bmu2Chi2ndf[size]/D");
  //BInfo.ujInfo
  bnt->Branch("Bmumumass",Bmumumass,"Bmumumass[size]/D");
  bnt->Branch("Bmumueta",Bmumueta,"Bmumueta[size]/D");
  bnt->Branch("Bmumuphi",Bmumuphi,"Bmumuphi[size]/D");
  bnt->Branch("Bmumuy",Bmumuy,"Bmumuy[size]/D");
  bnt->Branch("Bmumupt",Bmumupt,"Bmumupt[size]/D");  
  bnt->Branch("Bujmass",Bujmass,"Bujmass[size]/D");
  bnt->Branch("BujvProb",BujvProb,"BujvProb[size]/D");
  bnt->Branch("Bujpt",Bujpt,"Bujpt[size]/D");
  bnt->Branch("Bujeta",Bujeta,"Bujeta[size]/D");
  bnt->Branch("Bujphi",Bujphi,"Bujphi[size]/D");
  bnt->Branch("Bujy",Bujy,"Bujy[size]/D");
  bnt->Branch("Bujlxy",Bujlxy,"Bujlxy[size]/D");
  //BInfo.trkInfo
  bnt->Branch("Btrk1Pt",Btrk1Pt,"Btrk1Pt[size]/D");
  bnt->Branch("Btrk2Pt",Btrk2Pt,"Btrk2Pt[size]/D");
  bnt->Branch("Btrk1Eta",Btrk1Eta,"Btrk1Eta[size]/D");  
  bnt->Branch("Btrk2Eta",Btrk2Eta,"Btrk2Eta[size]/D");  
  bnt->Branch("Btrk1Phi",Btrk1Phi,"Btrk1Phi[size]/D");  
  bnt->Branch("Btrk2Phi",Btrk2Phi,"Btrk2Phi[size]/D");  
  bnt->Branch("Btrk1Y",Btrk1Y,"Btrk1Y[size]/D");  
  bnt->Branch("Btrk2Y",Btrk2Y,"Btrk2Y[size]/D");  
  bnt->Branch("Btrk1Dxy",Btrk1Dxy,"Btrk1Dxy[size]/D");
  bnt->Branch("Btrk2Dxy",Btrk2Dxy,"Btrk2Dxy[size]/D");
  bnt->Branch("Btrk1D0Err",Btrk1D0Err,"Btrk1D0Err[size]/D");
  bnt->Branch("Btrk2D0Err",Btrk2D0Err,"Btrk2D0Err[size]/D");
  bnt->Branch("Btrk1PixelHit",Btrk1PixelHit,"Btrk1PixelHit[size]/D");
  bnt->Branch("Btrk2PixelHit",Btrk2PixelHit,"Btrk2PixelHit[size]/D");
  bnt->Branch("Btrk1StripHit",Btrk1StripHit,"Btrk1StripHit[size]/D");
  bnt->Branch("Btrk2StripHit",Btrk2StripHit,"Btrk2StripHit[size]/D");
  bnt->Branch("Btrk1Chi2ndf",Btrk1Chi2ndf,"Btrk1Chi2ndf[size]/D");
  bnt->Branch("Btrk2Chi2ndf",Btrk2Chi2ndf,"Btrk2Chi2ndf[size]/D");
  //BInfo.tktkInfo
  bnt->Branch("Btktkmass",Btktkmass,"Btktkmass[size]/D");
  bnt->Branch("BtktkmassKK",BtktkmassKK,"BtktkmassKK[size]/D");
  bnt->Branch("BtktkvProb",BtktkvProb,"BtktkvProb[size]/D");
  bnt->Branch("Btktkpt",Btktkpt,"Btktkpt[size]/D");
  bnt->Branch("Btktketa",Btktketa,"Btktketa[size]/D");
  bnt->Branch("Btktkphi",Btktkphi,"Btktkphi[size]/D");
  bnt->Branch("Btktky",Btktky,"Btktky[size]/D");
  bnt->Branch("Bdoubletmass",Bdoubletmass,"Bdoubletmass[size]/D");
  bnt->Branch("Bdoubletpt",Bdoubletpt,"Bdoubletpt[size]/D");
  bnt->Branch("Bdoubleteta",Bdoubleteta,"Bdoubleteta[size]/D");  
  bnt->Branch("Bdoubletphi",Bdoubletphi,"Bdoubletphi[size]/D");  
  bnt->Branch("Bdoublety",Bdoublety,"Bdoublety[size]/D");
  //BInfo.genInfo
  bnt->Branch("Bgen",Bgen,"Bgen[size]/D");
  bnt->Branch("BgenIndex",BgenIndex,"BgenIndex[size]/I");
  bnt->Branch("Bgenpt",Bgenpt,"Bgenpt[size]/D");
  bnt->Branch("Bgeneta",Bgeneta,"Bgeneta[size]/D");
  bnt->Branch("Bgenphi",Bgenphi,"Bgenphi[size]/D");
  bnt->Branch("Bgeny",Bgeny,"Bgeny[size]/D");
  bnt->Branch("BgenpdgId",BgenpdgId,"BgenpdgId[size]/I");
}

//DInfo
Int_t      Dindex[MAX_XB];
//Int_t      DisGoodCand[MAX_XB];
Double_t   Dmass[MAX_XB];
Double_t   Dpt[MAX_XB];
Double_t   Deta[MAX_XB];
Double_t   Dphi[MAX_XB];
Double_t   Dy[MAX_XB];
Double_t   DvtxX[MAX_XB];
Double_t   DvtxY[MAX_XB];
Double_t   Dd0[MAX_XB];
Double_t   Dd0Err[MAX_XB];
Double_t   Ddxyz[MAX_XB];
Double_t   DdxyzErr[MAX_XB];
Double_t   Dchi2ndf[MAX_XB];
Double_t   Dchi2cl[MAX_XB];
Double_t   Ddtheta[MAX_XB];
Double_t   Dlxy[MAX_XB];
Bool_t     Disbestchi2[MAX_XB];
//DInfo.b4fitInfo
Double_t   Db4fit_mass[MAX_XB];
Double_t   Db4fit_pt[MAX_XB];
Double_t   Db4fit_eta[MAX_XB];
Double_t   Db4fit_phi[MAX_XB];
//DInfo.trkInfo
Double_t   Dtrk1Pt[MAX_XB];
Double_t   Dtrk2Pt[MAX_XB];
Double_t   Dtrk3Pt[MAX_XB];
Double_t   Dtrk4Pt[MAX_XB];
Double_t   Dtrk1Eta[MAX_XB];
Double_t   Dtrk2Eta[MAX_XB];
Double_t   Dtrk3Eta[MAX_XB];
Double_t   Dtrk4Eta[MAX_XB];
Double_t   Dtrk1Phi[MAX_XB];
Double_t   Dtrk2Phi[MAX_XB];
Double_t   Dtrk3Phi[MAX_XB];
Double_t   Dtrk4Phi[MAX_XB];
Double_t   Dtrk1PtErr[MAX_XB];
Double_t   Dtrk2PtErr[MAX_XB];
Double_t   Dtrk3PtErr[MAX_XB];
Double_t   Dtrk4PtErr[MAX_XB];
Double_t   Dtrk1EtaErr[MAX_XB];
Double_t   Dtrk2EtaErr[MAX_XB];
Double_t   Dtrk3EtaErr[MAX_XB];
Double_t   Dtrk4EtaErr[MAX_XB];
Double_t   Dtrk1PhiErr[MAX_XB];
Double_t   Dtrk2PhiErr[MAX_XB];
Double_t   Dtrk3PhiErr[MAX_XB];
Double_t   Dtrk4PhiErr[MAX_XB];
Double_t   Dtrk1Y[MAX_XB];
Double_t   Dtrk2Y[MAX_XB];
Double_t   Dtrk3Y[MAX_XB];
Double_t   Dtrk4Y[MAX_XB];
Double_t   Dtrk1Dxy[MAX_XB];
Double_t   Dtrk2Dxy[MAX_XB];
Double_t   Dtrk3Dxy[MAX_XB];
Double_t   Dtrk4Dxy[MAX_XB];
Double_t   Dtrk1D0Err[MAX_XB];
Double_t   Dtrk2D0Err[MAX_XB];
Double_t   Dtrk3D0Err[MAX_XB];
Double_t   Dtrk4D0Err[MAX_XB];
Double_t   Dtrk1PixelHit[MAX_XB];
Double_t   Dtrk2PixelHit[MAX_XB];
Double_t   Dtrk3PixelHit[MAX_XB];
Double_t   Dtrk4PixelHit[MAX_XB];
Double_t   Dtrk1StripHit[MAX_XB];
Double_t   Dtrk2StripHit[MAX_XB];
Double_t   Dtrk3StripHit[MAX_XB];
Double_t   Dtrk4StripHit[MAX_XB];
Double_t   Dtrk1nStripLayer[MAX_XB];
Double_t   Dtrk2nStripLayer[MAX_XB];
Double_t   Dtrk3nStripLayer[MAX_XB];
Double_t   Dtrk4nStripLayer[MAX_XB];
Double_t   Dtrk1nPixelLayer[MAX_XB];
Double_t   Dtrk2nPixelLayer[MAX_XB];
Double_t   Dtrk3nPixelLayer[MAX_XB];
Double_t   Dtrk4nPixelLayer[MAX_XB];
Double_t   Dtrk1Chi2ndf[MAX_XB];
Double_t   Dtrk2Chi2ndf[MAX_XB];
Double_t   Dtrk3Chi2ndf[MAX_XB];
Double_t   Dtrk4Chi2ndf[MAX_XB];
Double_t   Dtrk1MassHypo[MAX_XB];
Double_t   Dtrk2MassHypo[MAX_XB];
Double_t   Dtrk3MassHypo[MAX_XB];
Double_t   Dtrk4MassHypo[MAX_XB];
//DInfo.tktkResInfo
Double_t   DtktkResmass[MAX_XB];
Double_t   DtktkRespt[MAX_XB];
Double_t   DtktkReseta[MAX_XB];
Double_t   DtktkResphi[MAX_XB];
//DInfo.genInfo
Double_t   Dgen[MAX_XB];
Int_t      DgenIndex[MAX_XB];
Double_t   Dgenpt[MAX_XB];
Double_t   Dgeneta[MAX_XB];
Double_t   Dgenphi[MAX_XB];
Double_t   Dgeny[MAX_XB];

void buildDBranch(TTree* dnt)
{
  //EvtInfo

  dnt->Branch("RunNo",&RunNo);
  dnt->Branch("EvtNo",&EvtNo);
  dnt->Branch("Dsize",&Dsize);
  dnt->Branch("PVx",&PVx);
  dnt->Branch("PVy",&PVy);
  dnt->Branch("PVz",&PVz);
  dnt->Branch("PVxE",&PVxE);
  dnt->Branch("PVyE",&PVyE);
  dnt->Branch("PVzE",&PVzE);
  dnt->Branch("PVnchi2",&PVnchi2);
  dnt->Branch("PVchi2",&PVchi2);

  /*
  //DInfo
  dnt->Branch("Dindex",Dindex,"Dindex[Dsize]/I");
  dnt->Branch("Dmass",Dmass,"Dmass[Dsize]/D");
  dnt->Branch("Dpt",Dpt,"Dpt[Dsize]/D");
  dnt->Branch("Deta",Deta,"Deta[Dsize]/D");
  dnt->Branch("Dphi",Dphi,"Dphi[Dsize]/D");
  dnt->Branch("Dy",Dy,"Dy[Dsize]/D");
  dnt->Branch("DvtxX",DvtxX,"DvtxX[Dsize]/D");
  dnt->Branch("DvtxY",DvtxY,"DvtxY[Dsize]/D");
  dnt->Branch("Dd0",Dd0,"Dd0[Dsize]/D");
  dnt->Branch("Dd0Err",Dd0Err,"Dd0Err[Dsize]/D");
  dnt->Branch("Ddxyz",Ddxyz,"Ddxyz[Dsize]/D");
  dnt->Branch("DdxyzErr",DdxyzErr,"DdxyzErr[Dsize]/D");
  dnt->Branch("Dchi2ndf",Dchi2ndf,"Dchi2ndf[Dsize]/D");
  dnt->Branch("Dchi2cl",Dchi2cl,"Dchi2cl[Dsize]/D");
  dnt->Branch("Dchi2cl",Dchi2cl,"Dchi2cl[Dsize]/D");
  dnt->Branch("Ddtheta",Ddtheta,"Ddtheta[Dsize]/D");
  dnt->Branch("Dlxy",Dlxy,"Dlxy[Dsize]/D");
  dnt->Branch("Disbestchi2",Disbestchi2,"Disbestchi2[Dsize]/O");
  //DInfo.b4fitInfo
  dnt->Branch("Db4fit_mass",Db4fit_mass,"Db4fit_mass[Dsize]/D");
  dnt->Branch("Db4fit_pt",Db4fit_pt,"Db4fit_pt[Dsize]/D");
  dnt->Branch("Db4fit_eta",Db4fit_eta,"Db4fit_eta[Dsize]/D");
  dnt->Branch("Db4fit_phi",Db4fit_phi,"Db4fit_phi[Dsize]/D");
  //DInfo.trkInfo
  dnt->Branch("Dtrk1Pt",Dtrk1Pt,"Dtrk1Pt[Dsize]/D");
  dnt->Branch("Dtrk2Pt",Dtrk2Pt,"Dtrk2Pt[Dsize]/D");
  dnt->Branch("Dtrk3Pt",Dtrk3Pt,"Dtrk3Pt[Dsize]/D");
  dnt->Branch("Dtrk4Pt",Dtrk4Pt,"Dtrk4Pt[Dsize]/D");
  dnt->Branch("Dtrk1Eta",Dtrk1Eta,"Dtrk1Eta[Dsize]/D");
  dnt->Branch("Dtrk2Eta",Dtrk2Eta,"Dtrk2Eta[Dsize]/D");
  dnt->Branch("Dtrk3Eta",Dtrk3Eta,"Dtrk3Eta[Dsize]/D");
  dnt->Branch("Dtrk4Eta",Dtrk4Eta,"Dtrk4Eta[Dsize]/D");
  dnt->Branch("Dtrk1Phi",Dtrk1Phi,"Dtrk1Phi[Dsize]/D");
  dnt->Branch("Dtrk2Phi",Dtrk2Phi,"Dtrk2Phi[Dsize]/D");
  dnt->Branch("Dtrk3Phi",Dtrk3Phi,"Dtrk3Phi[Dsize]/D");
  dnt->Branch("Dtrk4Phi",Dtrk4Phi,"Dtrk4Phi[Dsize]/D");
  dnt->Branch("Dtrk1PtErr",Dtrk1PtErr,"Dtrk1PtErr[Dsize]/D");
  dnt->Branch("Dtrk2PtErr",Dtrk2PtErr,"Dtrk2PtErr[Dsize]/D");
  dnt->Branch("Dtrk3PtErr",Dtrk3PtErr,"Dtrk3PtErr[Dsize]/D");
  dnt->Branch("Dtrk4PtErr",Dtrk4PtErr,"Dtrk4PtErr[Dsize]/D");
  dnt->Branch("Dtrk1EtaErr",Dtrk1EtaErr,"Dtrk1EtaErr[Dsize]/D");
  dnt->Branch("Dtrk2EtaErr",Dtrk2EtaErr,"Dtrk2EtaErr[Dsize]/D");
  dnt->Branch("Dtrk3EtaErr",Dtrk3EtaErr,"Dtrk3EtaErr[Dsize]/D");
  dnt->Branch("Dtrk4EtaErr",Dtrk4EtaErr,"Dtrk4EtaErr[Dsize]/D");
  dnt->Branch("Dtrk1PhiErr",Dtrk1PhiErr,"Dtrk1PhiErr[Dsize]/D");
  dnt->Branch("Dtrk2PhiErr",Dtrk2PhiErr,"Dtrk2PhiErr[Dsize]/D");
  dnt->Branch("Dtrk3PhiErr",Dtrk3PhiErr,"Dtrk3PhiErr[Dsize]/D");
  dnt->Branch("Dtrk4PhiErr",Dtrk4PhiErr,"Dtrk4PhiErr[Dsize]/D");
  dnt->Branch("Dtrk1Y",Dtrk1Y,"Dtrk1Y[Dsize]/D");
  dnt->Branch("Dtrk2Y",Dtrk2Y,"Dtrk2Y[Dsize]/D");
  dnt->Branch("Dtrk3Y",Dtrk3Y,"Dtrk3Y[Dsize]/D");
  dnt->Branch("Dtrk4Y",Dtrk4Y,"Dtrk4Y[Dsize]/D");
  dnt->Branch("Dtrk1Dxy",Dtrk1Dxy,"Dtrk1Dxy[Dsize]/D");
  dnt->Branch("Dtrk2Dxy",Dtrk2Dxy,"Dtrk2Dxy[Dsize]/D");
  dnt->Branch("Dtrk3Dxy",Dtrk3Dxy,"Dtrk3Dxy[Dsize]/D");
  dnt->Branch("Dtrk4Dxy",Dtrk4Dxy,"Dtrk4Dxy[Dsize]/D");
  dnt->Branch("Dtrk1D0Err",Dtrk1D0Err,"Dtrk1D0Err[Dsize]/D");
  dnt->Branch("Dtrk2D0Err",Dtrk2D0Err,"Dtrk2D0Err[Dsize]/D");
  dnt->Branch("Dtrk3D0Err",Dtrk3D0Err,"Dtrk3D0Err[Dsize]/D");
  dnt->Branch("Dtrk4D0Err",Dtrk4D0Err,"Dtrk4D0Err[Dsize]/D");
  dnt->Branch("Dtrk1PixelHit",Dtrk1PixelHit,"Dtrk1PixelHit[Dsize]/D");
  dnt->Branch("Dtrk2PixelHit",Dtrk2PixelHit,"Dtrk2PixelHit[Dsize]/D");
  dnt->Branch("Dtrk3PixelHit",Dtrk3PixelHit,"Dtrk3PixelHit[Dsize]/D");
  dnt->Branch("Dtrk4PixelHit",Dtrk4PixelHit,"Dtrk4PixelHit[Dsize]/D");
  dnt->Branch("Dtrk1StripHit",Dtrk1StripHit,"Dtrk1StripHit[Dsize]/D");
  dnt->Branch("Dtrk2StripHit",Dtrk2StripHit,"Dtrk2StripHit[Dsize]/D");
  dnt->Branch("Dtrk3StripHit",Dtrk3StripHit,"Dtrk3StripHit[Dsize]/D");
  dnt->Branch("Dtrk4StripHit",Dtrk4StripHit,"Dtrk4StripHit[Dsize]/D");
  dnt->Branch("Dtrk1nStripLayer",Dtrk1nStripLayer,"Dtrk1nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk2nStripLayer",Dtrk2nStripLayer,"Dtrk2nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk3nStripLayer",Dtrk3nStripLayer,"Dtrk3nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk4nStripLayer",Dtrk4nStripLayer,"Dtrk4nStripLayer[Dsize]/D");
  dnt->Branch("Dtrk1nPixelLayer",Dtrk1nPixelLayer,"Dtrk1nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk2nPixelLayer",Dtrk2nPixelLayer,"Dtrk2nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk3nPixelLayer",Dtrk3nPixelLayer,"Dtrk3nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk4nPixelLayer",Dtrk4nPixelLayer,"Dtrk4nPixelLayer[Dsize]/D");
  dnt->Branch("Dtrk1Chi2ndf",Dtrk1Chi2ndf,"Dtrk1Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk2Chi2ndf",Dtrk2Chi2ndf,"Dtrk2Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk3Chi2ndf",Dtrk3Chi2ndf,"Dtrk3Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk4Chi2ndf",Dtrk4Chi2ndf,"Dtrk4Chi2ndf[Dsize]/D");
  dnt->Branch("Dtrk1MassHypo",Dtrk1MassHypo,"Dtrk1MassHypo[Dsize]/D");
  dnt->Branch("Dtrk2MassHypo",Dtrk2MassHypo,"Dtrk2MassHypo[Dsize]/D");
  dnt->Branch("Dtrk3MassHypo",Dtrk3MassHypo,"Dtrk3MassHypo[Dsize]/D");
  dnt->Branch("Dtrk4MassHypo",Dtrk4MassHypo,"Dtrk4MassHypo[Dsize]/D");
  //DInfo.tktkResInfo
  dnt->Branch("DtktkResmass",DtktkResmass,"DtktkResmass[Dsize]/D");
  dnt->Branch("DtktkRespt",DtktkRespt,"DtktkRespt[Dsize]/D");
  dnt->Branch("DtktkReseta",DtktkReseta,"DtktkReseta[Dsize]/D");
  dnt->Branch("DtktkResphi",DtktkResphi,"DtktkResphi[Dsize]/D");
  //DInfo.genInfo
  */
  dnt->Branch("Dgen",Dgen,"Dgen[Dsize]/D");
  /*
  dnt->Branch("DgenIndex",DgenIndex,"DgenIndex[Dsize]/I");
  dnt->Branch("Dgenpt",Dgenpt,"Dgenpt[Dsize]/D");
  dnt->Branch("Dgeneta",Dgeneta,"Dgeneta[Dsize]/D");
  dnt->Branch("Dgenphi",Dgenphi,"Dgenphi[Dsize]/D");
  dnt->Branch("Dgeny",Dgeny,"Dgeny[Dsize]/D");
  */
}

/*
//GenInfo
Int_t Gsize;
Double_t Gy[MAX_GEN];
Double_t Geta[MAX_GEN];
Double_t Gphi[MAX_GEN];
Double_t Gpt[MAX_GEN];
Double_t GpdgId[MAX_GEN];
Int_t GisSignal[MAX_GEN];

void buildGenBranch(TTree* nt)
{
  nt->Branch("Gsize",&Gsize);
  nt->Branch("Gy",Gy,"Gy[Gsize]/D");
  nt->Branch("Geta",Geta,"Geta[Gsize]/D");
  nt->Branch("Gphi",Gphi,"Gphi[Gsize]/D");
  nt->Branch("Gpt",Gpt,"Gpt[Gsize]/D");
  nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/D");
  nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/I");
}
*/

//EvtInfo
Int_t           EvtInfo_RunNo;
Int_t           EvtInfo_EvtNo;
Int_t           EvtInfo_BxNo;
Int_t           EvtInfo_LumiNo;
Int_t           EvtInfo_Orbit;
Bool_t          EvtInfo_McFlag;
Int_t           EvtInfo_nBX;
Int_t           EvtInfo_BXPU[MAX_BX];
Int_t           EvtInfo_nPU[MAX_BX];
Float_t         EvtInfo_trueIT[MAX_BX];
Double_t        EvtInfo_PVx;
Double_t        EvtInfo_PVy;
Double_t        EvtInfo_PVz;
Double_t        EvtInfo_PVxE;
Double_t        EvtInfo_PVyE;
Double_t        EvtInfo_PVzE;
Double_t        EvtInfo_PVnchi2;
Double_t        EvtInfo_PVchi2;

//TrackInfo
Int_t           TrackInfo_size;
Int_t           TrackInfo_index[MAX_TRACK];
Int_t           TrackInfo_handle_index[MAX_TRACK];
Int_t           TrackInfo_charge[MAX_TRACK];
Double_t        TrackInfo_pt[MAX_TRACK];
Double_t        TrackInfo_eta[MAX_TRACK];
Double_t        TrackInfo_phi[MAX_TRACK];
Double_t        TrackInfo_ptErr[MAX_TRACK];
Double_t        TrackInfo_etaErr[MAX_TRACK];
Double_t        TrackInfo_phiErr[MAX_TRACK];
Int_t           TrackInfo_striphit[MAX_TRACK];
Int_t           TrackInfo_pixelhit[MAX_TRACK];
Int_t           TrackInfo_nStripLayer[MAX_TRACK];
Int_t           TrackInfo_nPixelLayer[MAX_TRACK];
Int_t           TrackInfo_fpbarrelhit[MAX_TRACK];
Int_t           TrackInfo_fpendcaphit[MAX_TRACK];
Double_t        TrackInfo_chi2[MAX_TRACK];
Double_t        TrackInfo_ndf[MAX_TRACK];
Double_t        TrackInfo_d0[MAX_TRACK];
Double_t        TrackInfo_d0error[MAX_TRACK];
Double_t        TrackInfo_dzPV[MAX_TRACK];
Double_t        TrackInfo_dxyPV[MAX_TRACK];
//Int_t           TrackInfo_isGoodCand[MAX_TRACK];
Int_t           TrackInfo_geninfo_index[MAX_TRACK];
Int_t           TrackInfo_trackQuality[MAX_TRACK];
Bool_t          TrackInfo_highPurity[MAX_TRACK];

//DInfo
Int_t           DInfo_size;
Int_t           DInfo_index[MAX_XB];
Int_t           DInfo_type[MAX_XB];    
Double_t        DInfo_b4fit_mass[MAX_XB];
Double_t        DInfo_b4fit_pt[MAX_XB];
Double_t        DInfo_b4fit_eta[MAX_XB];
Double_t        DInfo_b4fit_phi[MAX_XB];
Double_t        DInfo_tktkRes_mass[MAX_XB];
Double_t        DInfo_tktkRes_pt[MAX_XB];
Double_t        DInfo_tktkRes_eta[MAX_XB];
Double_t        DInfo_tktkRes_phi[MAX_XB];
Double_t        DInfo_tktkRes_px[MAX_XB];
Double_t        DInfo_tktkRes_py[MAX_XB];
Double_t        DInfo_tktkRes_pz[MAX_XB];
Double_t        DInfo_tktkRes_vtxX[MAX_XB];
Double_t        DInfo_tktkRes_vtxY[MAX_XB];
Double_t        DInfo_tktkRes_vtxZ[MAX_XB];
Double_t        DInfo_tktkRes_vtxXE[MAX_XB];
Double_t        DInfo_tktkRes_vtxYE[MAX_XB];
Double_t        DInfo_tktkRes_vtxZE[MAX_XB];
Double_t        DInfo_tktkRes_vtxdof[MAX_XB];
Double_t        DInfo_tktkRes_vtxchi2[MAX_XB];
Double_t        DInfo_tktkRes_rftk1_pt[MAX_XB];
Double_t        DInfo_tktkRes_rftk1_eta[MAX_XB];
Double_t        DInfo_tktkRes_rftk1_phi[MAX_XB];
Double_t        DInfo_tktkRes_rftk2_pt[MAX_XB];
Double_t        DInfo_tktkRes_rftk2_eta[MAX_XB];
Double_t        DInfo_tktkRes_rftk2_phi[MAX_XB];
Double_t        DInfo_mass[MAX_XB];
Double_t        DInfo_pt[MAX_XB];
Double_t        DInfo_eta[MAX_XB];
Double_t        DInfo_phi[MAX_XB];
Double_t        DInfo_px[MAX_XB];
Double_t        DInfo_py[MAX_XB];
Double_t        DInfo_pz[MAX_XB];
Double_t        DInfo_vtxX[MAX_XB];
Double_t        DInfo_vtxY[MAX_XB];
Double_t        DInfo_vtxZ[MAX_XB];
Double_t        DInfo_vtxXE[MAX_XB];
Double_t        DInfo_vtxYE[MAX_XB];
Double_t        DInfo_vtxZE[MAX_XB];
Double_t        DInfo_vtxdof[MAX_XB];
Double_t        DInfo_vtxchi2[MAX_XB];
Double_t        DInfo_rftk1_px[MAX_XB];
Double_t        DInfo_rftk1_py[MAX_XB];
Double_t        DInfo_rftk1_pz[MAX_XB];
Double_t        DInfo_rftk2_px[MAX_XB];
Double_t        DInfo_rftk2_py[MAX_XB];
Double_t        DInfo_rftk2_pz[MAX_XB];
Double_t        DInfo_rftk3_px[MAX_XB];
Double_t        DInfo_rftk3_py[MAX_XB];
Double_t        DInfo_rftk3_pz[MAX_XB];
Double_t        DInfo_rftk4_px[MAX_XB];
Double_t        DInfo_rftk4_py[MAX_XB];
Double_t        DInfo_rftk4_pz[MAX_XB];
Double_t        DInfo_rftk1_pt[MAX_XB];
Double_t        DInfo_rftk1_eta[MAX_XB];
Double_t        DInfo_rftk1_phi[MAX_XB];
Double_t        DInfo_rftk2_pt[MAX_XB];
Double_t        DInfo_rftk2_eta[MAX_XB];
Double_t        DInfo_rftk2_phi[MAX_XB];
Double_t        DInfo_rftk3_pt[MAX_XB];
Double_t        DInfo_rftk3_eta[MAX_XB];
Double_t        DInfo_rftk3_phi[MAX_XB];
Double_t        DInfo_rftk4_pt[MAX_XB];
Double_t        DInfo_rftk4_eta[MAX_XB];
Double_t        DInfo_rftk4_phi[MAX_XB];
Int_t           DInfo_rftk1_index[MAX_XB];
Int_t           DInfo_rftk2_index[MAX_XB];
Int_t           DInfo_rftk3_index[MAX_XB];
Int_t           DInfo_rftk4_index[MAX_XB];
Int_t           DInfo_rftk1_MassHypo[MAX_XB];
Int_t           DInfo_rftk2_MassHypo[MAX_XB];
Int_t           DInfo_rftk3_MassHypo[MAX_XB];
Int_t           DInfo_rftk4_MassHypo[MAX_XB];

//GenInfo
Int_t           GenInfo_size;
Int_t           GenInfo_index[MAX_GEN];
Int_t           GenInfo_handle_index[MAX_GEN];
Double_t        GenInfo_pt[MAX_GEN];
Double_t        GenInfo_eta[MAX_GEN];
Double_t        GenInfo_phi[MAX_GEN];
Double_t        GenInfo_mass[MAX_GEN];
Int_t           GenInfo_pdgId[MAX_GEN];
Int_t           GenInfo_status[MAX_GEN];
Int_t           GenInfo_nMo[MAX_GEN];
Int_t           GenInfo_nDa[MAX_GEN];
Int_t           GenInfo_mo1[MAX_GEN];
Int_t           GenInfo_mo2[MAX_GEN];
Int_t           GenInfo_da1[MAX_GEN];
Int_t           GenInfo_da2[MAX_GEN];


void setDBranch(TTree *root)
{
  //EvtInfo
  root->SetBranchAddress("EvtInfo.RunNo",&EvtInfo_RunNo);
  root->SetBranchAddress("EvtInfo.EvtNo",&EvtInfo_EvtNo);
  /*
  root->SetBranchAddress("EvtInfo.BxNo",&EvtInfo_BxNo);
  root->SetBranchAddress("EvtInfo.LumiNo",&EvtInfo_LumiNo);
  root->SetBranchAddress("EvtInfo.Orbit",&EvtInfo_Orbit);
  root->SetBranchAddress("EvtInfo.McFlag",&EvtInfo_McFlag);
  root->SetBranchAddress("EvtInfo.nBX",&EvtInfo_nBX);
  root->SetBranchAddress("EvtInfo.BXPU",EvtInfo_BXPU);
  root->SetBranchAddress("EvtInfo.nPU",EvtInfo_nPU);
  root->SetBranchAddress("EvtInfo.trueIT",EvtInfo_trueIT);
  */
  root->SetBranchAddress("EvtInfo.PVx",&EvtInfo_PVx);
  root->SetBranchAddress("EvtInfo.PVy",&EvtInfo_PVy);
  root->SetBranchAddress("EvtInfo.PVz",&EvtInfo_PVz);
  /*
  root->SetBranchAddress("EvtInfo.PVxE",&EvtInfo_PVxE);
  root->SetBranchAddress("EvtInfo.PVyE",&EvtInfo_PVyE);
  root->SetBranchAddress("EvtInfo.PVzE",&EvtInfo_PVzE);
  root->SetBranchAddress("EvtInfo.PVnchi2",&EvtInfo_PVnchi2);
  root->SetBranchAddress("EvtInfo.PVchi2",&EvtInfo_PVchi2);
  */
  //TrackInfo
  /*
  root->SetBranchAddress("TrackInfo.size",&TrackInfo_size);
  root->SetBranchAddress("TrackInfo.index",TrackInfo_index);
  root->SetBranchAddress("TrackInfo.handle_index",TrackInfo_handle_index);
  root->SetBranchAddress("TrackInfo.charge",TrackInfo_charge);
  */
  root->SetBranchAddress("TrackInfo.pt",TrackInfo_pt);
  /*
  root->SetBranchAddress("TrackInfo.eta",TrackInfo_eta);
  root->SetBranchAddress("TrackInfo.phi",TrackInfo_phi);
  root->SetBranchAddress("TrackInfo.ptErr",TrackInfo_ptErr);
  root->SetBranchAddress("TrackInfo.etaErr",TrackInfo_etaErr);
  root->SetBranchAddress("TrackInfo.phiErr",TrackInfo_phiErr);
  root->SetBranchAddress("TrackInfo.striphit",TrackInfo_striphit);
  root->SetBranchAddress("TrackInfo.pixelhit",TrackInfo_pixelhit);
  root->SetBranchAddress("TrackInfo.nStripLayer",TrackInfo_nStripLayer);
  root->SetBranchAddress("TrackInfo.nPixelLayer",TrackInfo_nPixelLayer);
  root->SetBranchAddress("TrackInfo.fpbarrelhit",TrackInfo_fpbarrelhit);
  root->SetBranchAddress("TrackInfo.fpendcaphit",TrackInfo_fpendcaphit);
  root->SetBranchAddress("TrackInfo.chi2",TrackInfo_chi2);
  root->SetBranchAddress("TrackInfo.ndf",TrackInfo_ndf);
  root->SetBranchAddress("TrackInfo.d0",TrackInfo_d0);
  root->SetBranchAddress("TrackInfo.d0error",TrackInfo_d0error);
  root->SetBranchAddress("TrackInfo.dzPV",TrackInfo_dzPV);
  root->SetBranchAddress("TrackInfo.dxyPV",TrackInfo_dxyPV);
  //root->SetBranchAddress("TrackInfo.isGoodCand",TrackInfo_isGoodCand);
  */
  root->SetBranchAddress("TrackInfo.geninfo_index",TrackInfo_geninfo_index);
  /*
  root->SetBranchAddress("TrackInfo.trackQuality",TrackInfo_trackQuality);
  root->SetBranchAddress("TrackInfo.highPurity",TrackInfo_highPurity);
  */
  //DInfo
  root->SetBranchAddress("DInfo.size",&DInfo_size);
  //root->SetBranchAddress("DInfo.index",DInfo_index);
  root->SetBranchAddress("DInfo.type",DInfo_type);
  /*
  root->SetBranchAddress("DInfo.b4fit_mass",DInfo_b4fit_mass);
  root->SetBranchAddress("DInfo.b4fit_pt",DInfo_b4fit_pt);
  root->SetBranchAddress("DInfo.b4fit_eta",DInfo_b4fit_eta);
  root->SetBranchAddress("DInfo.b4fit_phi",DInfo_b4fit_phi);
  root->SetBranchAddress("DInfo.tktkRes_mass",DInfo_tktkRes_mass);
  root->SetBranchAddress("DInfo.tktkRes_pt",DInfo_tktkRes_pt);
  root->SetBranchAddress("DInfo.tktkRes_eta",DInfo_tktkRes_eta);
  root->SetBranchAddress("DInfo.tktkRes_phi",DInfo_tktkRes_phi);
  root->SetBranchAddress("DInfo.tktkRes_px",DInfo_tktkRes_px);
  root->SetBranchAddress("DInfo.tktkRes_py",DInfo_tktkRes_py);
  root->SetBranchAddress("DInfo.tktkRes_pz",DInfo_tktkRes_pz);
  root->SetBranchAddress("DInfo.tktkRes_vtxX",DInfo_tktkRes_vtxX);
  root->SetBranchAddress("DInfo.tktkRes_vtxY",DInfo_tktkRes_vtxY);
  root->SetBranchAddress("DInfo.tktkRes_vtxZ",DInfo_tktkRes_vtxZ);
  root->SetBranchAddress("DInfo.tktkRes_vtxXE",DInfo_tktkRes_vtxXE);
  root->SetBranchAddress("DInfo.tktkRes_vtxYE",DInfo_tktkRes_vtxYE);
  root->SetBranchAddress("DInfo.tktkRes_vtxZE",DInfo_tktkRes_vtxZE);
  root->SetBranchAddress("DInfo.tktkRes_vtxdof",DInfo_tktkRes_vtxdof);
  root->SetBranchAddress("DInfo.tktkRes_vtxchi2",DInfo_tktkRes_vtxchi2);
  root->SetBranchAddress("DInfo.tktkRes_rftk1_pt",DInfo_tktkRes_rftk1_pt);
  root->SetBranchAddress("DInfo.tktkRes_rftk1_eta",DInfo_tktkRes_rftk1_eta);
  root->SetBranchAddress("DInfo.tktkRes_rftk1_phi",DInfo_tktkRes_rftk1_phi);
  root->SetBranchAddress("DInfo.tktkRes_rftk2_pt",DInfo_tktkRes_rftk2_pt);
  root->SetBranchAddress("DInfo.tktkRes_rftk2_eta",DInfo_tktkRes_rftk2_eta);
  root->SetBranchAddress("DInfo.tktkRes_rftk2_phi",DInfo_tktkRes_rftk2_phi);
  */
  root->SetBranchAddress("DInfo.mass",DInfo_mass);
  root->SetBranchAddress("DInfo.pt",DInfo_pt);
  /*
  root->SetBranchAddress("DInfo.eta",DInfo_eta);
  root->SetBranchAddress("DInfo.phi",DInfo_phi);
  root->SetBranchAddress("DInfo.px",DInfo_px);
  root->SetBranchAddress("DInfo.py",DInfo_py);
  root->SetBranchAddress("DInfo.pz",DInfo_pz);
  */
  root->SetBranchAddress("DInfo.vtxX",DInfo_vtxX);
  root->SetBranchAddress("DInfo.vtxY",DInfo_vtxY);
  root->SetBranchAddress("DInfo.vtxZ",DInfo_vtxZ);
  root->SetBranchAddress("DInfo.vtxXE",DInfo_vtxXE);
  root->SetBranchAddress("DInfo.vtxYE",DInfo_vtxYE);
  root->SetBranchAddress("DInfo.vtxZE",DInfo_vtxZE);
  /*
  root->SetBranchAddress("DInfo.vtxdof",DInfo_vtxdof);
  root->SetBranchAddress("DInfo.vtxchi2",DInfo_vtxchi2);
  root->SetBranchAddress("DInfo.rftk1_px",DInfo_rftk1_px);
  root->SetBranchAddress("DInfo.rftk1_py",DInfo_rftk1_py);
  root->SetBranchAddress("DInfo.rftk1_pz",DInfo_rftk1_pz);
  root->SetBranchAddress("DInfo.rftk2_px",DInfo_rftk2_px);
  root->SetBranchAddress("DInfo.rftk2_py",DInfo_rftk2_py);
  root->SetBranchAddress("DInfo.rftk2_pz",DInfo_rftk2_pz);
  root->SetBranchAddress("DInfo.rftk3_px",DInfo_rftk3_px);
  root->SetBranchAddress("DInfo.rftk3_py",DInfo_rftk3_py);
  root->SetBranchAddress("DInfo.rftk3_pz",DInfo_rftk3_pz);
  root->SetBranchAddress("DInfo.rftk4_px",DInfo_rftk4_px);
  root->SetBranchAddress("DInfo.rftk4_py",DInfo_rftk4_py);
  root->SetBranchAddress("DInfo.rftk4_pz",DInfo_rftk4_pz);
  */
  root->SetBranchAddress("DInfo.rftk1_pt",DInfo_rftk1_pt);
  root->SetBranchAddress("DInfo.rftk1_eta",DInfo_rftk1_eta);
  root->SetBranchAddress("DInfo.rftk1_phi",DInfo_rftk1_phi);
  root->SetBranchAddress("DInfo.rftk2_pt",DInfo_rftk2_pt);
  root->SetBranchAddress("DInfo.rftk2_eta",DInfo_rftk2_eta);
  root->SetBranchAddress("DInfo.rftk2_phi",DInfo_rftk2_phi);
  /*
  root->SetBranchAddress("DInfo.rftk3_pt",DInfo_rftk3_pt);
  root->SetBranchAddress("DInfo.rftk3_eta",DInfo_rftk3_eta);
  root->SetBranchAddress("DInfo.rftk3_phi",DInfo_rftk3_phi);
  root->SetBranchAddress("DInfo.rftk4_pt",DInfo_rftk4_pt);
  root->SetBranchAddress("DInfo.rftk4_eta",DInfo_rftk4_eta);
  root->SetBranchAddress("DInfo.rftk4_phi",DInfo_rftk4_phi);
  */
  root->SetBranchAddress("DInfo.rftk1_index",DInfo_rftk1_index);
  root->SetBranchAddress("DInfo.rftk2_index",DInfo_rftk2_index);
  root->SetBranchAddress("DInfo.rftk3_index",DInfo_rftk3_index);
  root->SetBranchAddress("DInfo.rftk4_index",DInfo_rftk4_index);
  root->SetBranchAddress("DInfo.rftk1_MassHypo",DInfo_rftk1_MassHypo);
  root->SetBranchAddress("DInfo.rftk2_MassHypo",DInfo_rftk2_MassHypo);
  root->SetBranchAddress("DInfo.rftk3_MassHypo",DInfo_rftk3_MassHypo);
  root->SetBranchAddress("DInfo.rftk4_MassHypo",DInfo_rftk4_MassHypo);
  //GenInfo
  //root->SetBranchAddress("GenInfo.size",&GenInfo_size);
  root->SetBranchAddress("GenInfo.index",GenInfo_index);
  //root->SetBranchAddress("GenInfo.handle_index",GenInfo_handle_index);
  root->SetBranchAddress("GenInfo.pt",GenInfo_pt);
  root->SetBranchAddress("GenInfo.eta",GenInfo_eta);
  root->SetBranchAddress("GenInfo.phi",GenInfo_phi);
  root->SetBranchAddress("GenInfo.mass",GenInfo_mass);
  root->SetBranchAddress("GenInfo.pdgId",GenInfo_pdgId);
  //root->SetBranchAddress("GenInfo.status",GenInfo_status);
  root->SetBranchAddress("GenInfo.nMo",GenInfo_nMo);
  root->SetBranchAddress("GenInfo.nDa",GenInfo_nDa);
  root->SetBranchAddress("GenInfo.mo1",GenInfo_mo1);
  root->SetBranchAddress("GenInfo.mo2",GenInfo_mo2);
  root->SetBranchAddress("GenInfo.da1",GenInfo_da1);
  root->SetBranchAddress("GenInfo.da2",GenInfo_da2);
}

