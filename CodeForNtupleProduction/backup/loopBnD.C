#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <iostream>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <cmath>
#include "loop.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916

Int_t PION_PDGID = 211;
Int_t KAON_PDGID = 321;
Int_t DZERO_PDGID = 421;
Int_t DPLUS_PDGID = 411;
Int_t DSUBS_PDGID = 431;

Bool_t writeBmeson = false;

//void loop(string infile="/export/d00/scratch/jwang/Bfinder_BoostedMC_20140418_Hijing_PPb502_MinimumBias_HIJINGemb_inclBtoPsiMuMu_5TeV.root", string outfile="/export/d00/scratch/jwang/jpsi.root", bool REAL=0){
int loop(string infile="/mnt/hadoop/cms/store/user/jwang/Bfinder_BoostedMC_20140707_BuJpsiK_pPb.root", string
	  outfile="/export/d00/scratch/jwang/nt_BoostedMC_20140708_BuJpsiK_pPb.root.root", bool REAL=1,bool PbpMC=0,int startEntries=0,int nEntries=0, bool doMuonSelection = 0){

  const char* infname;
  const char* outfname;

  if(REAL) cout<<"--- REAL DATA ---"<<endl;
  else cout<<"--- MC ---"<<endl;

  infname = infile.c_str();
  outfname = outfile.c_str();

  //File type
  TFile *f = new TFile(infname);
  TTree *root = (TTree*)f->Get("demo/root");
  TFile *outf = new TFile(outfname,"recreate");
  setBranch(root);
  
  int isBchannel[7];
  isBchannel[0] = 0; //jpsi+Kp
  isBchannel[1] = 0; //jpsi+pi
  isBchannel[2] = 0; //jpsi+Ks(pi+,pi-)
  isBchannel[3] = 0; //jpsi+K*(K+,pi-)
  isBchannel[4] = 0; //jpsi+K*(K-,pi+)
  isBchannel[5] = 0; //jpsi+phi
  isBchannel[6] = 0; //jpsi+pi pi <= psi', X(3872), Bs->J/psi f0

  int isDchannel[5];
  isDchannel[0] = 1; //
  isDchannel[1] = 1; //
  isDchannel[2] = 0; //
  isDchannel[3] = 0; //
  isDchannel[4] = 0; //

  cout<<"--- Building trees ---"<<endl;
  TTree* ntD1 = new TTree("ntDkPpiM","");       buildDBranch(ntD1);
  TTree* ntD2 = new TTree("ntDkMpiP","");       buildDBranch(ntD2);
  TTree* ntD3 = new TTree("ntDkMpiPpiP","");    buildDBranch(ntD3);
  TTree* ntD4 = new TTree("ntDkPpiMpiM","");    buildDBranch(ntD4);
  TTree* ntD5 = new TTree("ntkMpiPpiPpiM","");  buildDBranch(ntD5);
  if(writeBmeson)
    {
      TTree* ntB0 = new TTree("ntKp","");     buildBBranch(nt0);
      TTree* ntB1 = new TTree("ntpi","");     buildBBranch(nt1);
      TTree* ntB2 = new TTree("ntKs","");     buildBBranch(nt2);
      TTree* ntB3 = new TTree("ntKstar","");  buildBBranch(ntB3);
      TTree* ntB5 = new TTree("ntphi","");    buildBBranch(ntB5);
      TTree* ntB6 = new TTree("ntmix","");    buildBBranch(ntB6);
      TTree* ntGen = new TTree("ntGen","");  buildGenBranch(ntGen);
    }
  cout<<"--- Building trees finished ---"<<endl;

  Long64_t nentries = root->GetEntries();
  //nentries = 10000;
  Long64_t nbytes = 0;
  TVector3* bP = new TVector3;
  TVector3* bVtx = new TVector3;
  TLorentzVector* b4P = new TLorentzVector;
  TLorentzVector* b4Pout = new TLorentzVector;
  TLorentzVector bGen;
  int type,flag;

  for (Long64_t i=startEntries; i<nentries;i++)
    {
      nbytes += root->GetEntry(i);
      if (i%10000==0) cout <<i<<" / "<<nentries<<endl;
      Int_t Dtypesize[5]={0,0,0,0,0};
      //Int_t Dtype1size=0,Dtype2size=0,Dtype3size=0,Dtype4size=0,Dtype5size=0,Dtype6size=0,Dtype7size=0;
      Double_t Dbest=0;
      Int_t Dbestindex=0;

      for(int t=0;t<5;t++)
	{
	  if(isDchannel[t]==1)
	    {
	      best=-1;
	      bestindex=-1;
	      for(int j=0;j<DInfo_size;j++)
		{
		  if(DInfo_type[j]==1)
		    {
		      fillDTree(bP,bVtx,b4P,j,Dtypesize[i],REAL,PbpMC);
		      if(chi2cl[Dtypesize[i]]>best)
			{
			  best = chi2cl[Dtypesize[i]];
			  bestindex = Dtypesize[i];
			}
		      Dtypesize[i]++;
		    }
		}
	      if(bestindex>-1)
		{
		  bestchi2 = bestindex;
		  isbestchi2[bestindex] = true;
		}
	      if(t==0) ntD1->Fill();
	      else if(t==1) ntD2->Fill();
	      else if(t==2) ntD3->Fill();
	      else if(t==3) ntD4->Fill();
	      else if(t==4) ntD5->Fill();
	    }
	}

      //B meson part
      if(writeBmeson)
	{
	  int Btype1size=0,Btype2size=0,Btype3size=0,Btype4size=0,Btype5size=0,Btype6size=0,Btype7size=0;
	  float best,best2,temy;
	  int bestindex,best2index;
	  
	  Bsize=0;
	  best=-1;
	  bestindex=-1;
	  best2=10000.;
	  best2index=-1;
	  for (int j=0;j<BInfo_size;j++) 
	    {
	      if(BInfo_type[j]>7) continue;
	      if(isBchannel[BInfo_type[j]-1]!=1) continue;
	      
	      //muonIdPass{{
	      if (doMuonSelection)
		{
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])) continue;
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096)) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096)) continue;
		}
	      //}}
	      b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	      temy = b4Pout->Rapidity();
	      if(REAL)
		{
		  if(!(((EvtInfo_RunNo>=210498&&EvtInfo_RunNo<=211256&&abs(temy)<2.4)||(EvtInfo_RunNo>=211313&&EvtInfo_RunNo<=211631&&abs(temy)<2.4))||(EvtInfo_RunNo>=211739&&EvtInfo_RunNo<=211831&&abs(temy)<2.4))) continue;
		}
	      else
		{
		  if((PbpMC==0)&&abs(temy)>=2.4) continue;
		  if((PbpMC==1)&&abs(temy)>=2.4) continue;
		}
	      if(BInfo_mass[j]<5 || BInfo_mass[j]>6) continue;
	      if(BInfo_pt[j]<10.) continue;
	      if(BInfo_type[j]==1)
		{
		  fillTree(bP,bVtx,b4P,j,Btype1size,KAON_MASS,0,REAL,PbpMC);
		  if(chi2cl[Btype1size]>best&&trk1Pt[Btype1size]>0.9&&HLT_PAMu3_v1&&abs(mumumass[Btype1size]-3.096916)<0.15&&chi2cl[Btype1size]>1.32e-02&&(d0[Btype1size]/d0Err[Btype1size])>3.41&&cos(dtheta[Btype1size])>-3.46e-01)
		    {
		      best = chi2cl[Btype1size];
		      bestindex = Btype1size;
		    }
		  Btype1size++;
		}
	    }
	  if(bestindex>-1)
	    {
	      bestchi2 = bestindex;
	      isbestchi2[bestindex] = 1;
	    }
	  nt0->Fill();
	  
	  Bsize=0;
	  best=-1;
	  bestindex=-1;
	  best2=10000.;
	  best2index=-1;
	  for (int j=0;j<BInfo_size;j++) 
	    {
	      if(BInfo_type[j]>7) continue;
	      if (isBchannel[BInfo_type[j]-1]!=1) continue;
	      if (doMuonSelection)
		{
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])) continue;
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])<6.) continue;
		  
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096)) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096)) continue;
		}
	      b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	      temy = b4Pout->Rapidity();
	      if(REAL)
		{
		  if(!(((EvtInfo_RunNo>=210498&&EvtInfo_RunNo<=211256&&abs(temy)<2.4)||(EvtInfo_RunNo>=211313&&EvtInfo_RunNo<=211631&&abs(temy)<2.4))||(EvtInfo_RunNo>=211739&&EvtInfo_RunNo<=211831&&abs(temy)<2.4))) continue;
		}
	      else
		{
		  if((PbpMC==0)&&abs(temy)>=2.4) continue;
		  if((PbpMC==1)&&abs(temy)>=2.4) continue;
		}
	      if(BInfo_mass[j]<5 || BInfo_mass[j]>6) continue;
	      if(BInfo_pt[j]<10.) continue;
	      if(BInfo_type[j]==2)
		{
		  fillTree(bP,bVtx,b4P,j,Btype2size,PION_MASS,0,REAL,PbpMC);
		  if(chi2cl[Btype2size]>best)
		    {
		      best = chi2cl[Btype2size];
		      bestindex = Btype2size;
		    }
		  Btype2size++;
		}
	    }
	  if(bestindex>-1)
	    {
	      bestchi2 = bestindex;
	      isbestchi2[bestindex] = 1;
	    }
	  nt1->Fill();
	  
	  Bsize=0;
	  best=-1;
	  bestindex=-1;
	  best2=10000.;
	  best2index=-1;
	  for (int j=0;j<BInfo_size;j++) 
	    {
	      if(BInfo_type[j]>7) continue;
	      if (isBchannel[BInfo_type[j]-1]!=1) continue;
	      
	      if (doMuonSelection)
		{
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])) continue;
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])<6.) continue;
		  
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096)) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096)) continue;
		}
	      b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	      temy = b4Pout->Rapidity();
	      if(REAL)
		{
		  if(!(((EvtInfo_RunNo>=210498&&EvtInfo_RunNo<=211256&&abs(temy)<2.4)||(EvtInfo_RunNo>=211313&&EvtInfo_RunNo<=211631&&abs(temy)<2.4))||(EvtInfo_RunNo>=211739&&EvtInfo_RunNo<=211831&&abs(temy)<2.4))) continue;
		}
	      else
		{
		  if((PbpMC==0)&&abs(temy)>=2.4) continue;
		  if((PbpMC==1)&&abs(temy)>=2.4) continue;
		}
	      if(BInfo_mass[j]<5 || BInfo_mass[j]>6) continue;
	      if(BInfo_pt[j]<10.) continue;
	      if(BInfo_type[j]==3)
		{
		  fillTree(bP,bVtx,b4P,j,Btype3size,PION_MASS,PION_MASS,REAL,PbpMC);
		  if(chi2cl[Btype3size]>best)
		    {
		      best = chi2cl[Btype3size];
		      bestindex = Btype3size;
		    }
		  if(abs(tktkmass[Btype3size]-KSHORT_MASS)<best2)
		    {
		      best2 = abs(tktkmass[Btype3size]-KSHORT_MASS);
		      best2index = Btype3size;
		    }
		  Btype3size++;
		}
	    }
	  if(size>0)
	    {
	      bestchi2 = bestindex;
	      isbestchi2[bestindex] = 1;
	      besttktkmass = best2index;
	      isbesttktkmass[best2index] = 1;
	    }
	  nt2->Fill();
	  
	  Bsize=0;
	  best=-1;
	  bestindex=-1;
	  best2=10000.;
	  best2index=-1;
	  for (int j=0;j<BInfo_size;j++) 
	    {
	      if(BInfo_type[j]>7) continue;
	      if (isBchannel[BInfo_type[j]-1]!=1) continue;
	      
	      //muonIdPass{{
	      if (doMuonSelection)
		{
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])) continue;
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096)) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096)) continue;
		}
	      b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	      temy = b4Pout->Rapidity();
	      if(REAL)
		{
		  if(!(((EvtInfo_RunNo>=210498&&EvtInfo_RunNo<=211256&&abs(temy)<2.4)||(EvtInfo_RunNo>=211313&&EvtInfo_RunNo<=211631&&abs(temy)<2.4))||(EvtInfo_RunNo>=211739&&EvtInfo_RunNo<=211831&&abs(temy)<2.4))) continue;
		}
	      else
		{
		  if((PbpMC==0)&&abs(temy)>=2.4) continue;
		  if((PbpMC==1)&&abs(temy)>=2.4) continue;
		}
	      if(BInfo_mass[j]<5 || BInfo_mass[j]>6) continue;
	      if(BInfo_pt[j]<10.) continue;
	      if(BInfo_type[j]==4 || BInfo_type[j]==5)
		{
		  fillTree(bP,bVtx,b4P,j,Btype4size,KAON_MASS,PION_MASS,REAL,PbpMC);
		  if(chi2cl[Btype4size]>best&&(HLT_PAMu3_v1)&&abs(mumumass[Btype4size]-3.096916)<0.15&&trk1Pt[Btype4size]>0.7&&trk2Pt[Btype4size]>0.7&&chi2cl[Btype4size]>9.94e-02&&(d0[Btype4size]/d0Err[Btype4size])>6.08&&cos(dtheta[Btype4size])>7.93e-01&&abs(tktkmass[Btype4size]-0.89594)<0.10&&tktkmassKK[Btype4size]>1.04)
		    {
		      best = chi2cl[Btype4size];
		      bestindex = Btype4size;
		    }
		  if(abs(tktkmass[Btype4size]-KSTAR_MASS)<best2&&(HLT_PAMu3_v1)&&abs(mumumass[Btype4size]-3.096916)<0.15&&trk1Pt[Btype4size]>0.7&&trk2Pt[Btype4size]>0.7&&chi2cl[Btype4size]>9.94e-02&&(d0[Btype4size]/d0Err[Btype4size])>6.08&&cos(dtheta[Btype4size])>7.93e-01&&abs(tktkmass[Btype4size]-0.89594)<0.10&&tktkmassKK[Btype4size]>1.04)
		    {
		      best2 = abs(tktkmass[Btype4size]-KSTAR_MASS);
		      best2index = Btype4size;
		    }
		  Btype4size++;
		}
	    }
	  if(bestindex>-1)
	    {
	      bestchi2 = bestindex;
	      isbestchi2[bestindex] = 1;
	    }
	  if(best2index>-1)
	    {
	      besttktkmass = best2index;
	      isbesttktkmass[best2index] = 1;
	    }
	  nt3->Fill();
	  
	  Bsize=0;
	  best=-1;
	  bestindex=-1;
	  best2=10000.;
	  best2index=-1;
	  for (int j=0;j<BInfo_size;j++) 
	    {
	      if(BInfo_type[j]>7) continue;
	      if (isBchannel[BInfo_type[j]-1]!=1) continue;
	      
	      if (doMuonSelection)
		{
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])) continue;
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])<6.) continue;
		  
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096)) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096)) continue;
		}
	      b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	      temy = b4Pout->Rapidity();
	      if(REAL)
		{
		  if(!(((EvtInfo_RunNo>=210498&&EvtInfo_RunNo<=211256&&abs(temy)<2.4)||(EvtInfo_RunNo>=211313&&EvtInfo_RunNo<=211631&&abs(temy)<2.4))||(EvtInfo_RunNo>=211739&&EvtInfo_RunNo<=211831&&abs(temy)<2.4))) continue;
		}
	      else
		{
		  if((PbpMC==0)&&abs(temy)>=2.4) continue;
		  if((PbpMC==1)&&abs(temy)>=2.4) continue;
		}
	      if(BInfo_mass[j]<5 || BInfo_mass[j]>6) continue;
	      if(BInfo_pt[j]<10.) continue;
	      if(BInfo_type[j]==6)
		{
		  fillTree(bP,bVtx,b4P,j,Btype6size,KAON_MASS,KAON_MASS,REAL,PbpMC);
		  if(chi2cl[Btype6size]>best&&(HLT_PAMu3_v1)&&abs(mumumass[Btype6size]-3.096916)<0.15&&trk1Pt[Btype6size]>0.7&&trk2Pt[Btype6size]>0.7&&chi2cl[Btype6size]>3.71e-02&&(d0[Btype6size]/d0Err[Btype6size])>3.37&&cos(dtheta[Btype6size])>2.60e-01&&abs(tktkmass[Btype6size]-1.019455)<1.55e-02)
		    {
		      best = chi2cl[Btype6size];
		      bestindex = Btype6size;
		    }
		  if(abs(tktkmass[Btype6size]-PHI_MASS)<best2&&(HLT_PAMu3_v1)&&abs(mumumass[Btype6size]-3.096916)<0.15&&trk1Pt[Btype6size]>0.7&&trk2Pt[Btype6size]>0.7&&chi2cl[Btype6size]>3.71e-02&&(d0[Btype6size]/d0Err[Btype6size])>3.37&&cos(dtheta[Btype6size])>2.60e-01&&abs(tktkmass[Btype6size]-1.019455)<1.55e-02)
		    {
		      best2 = abs(tktkmass[Btype6size]-PHI_MASS);
		      best2index = Btype6size;
		    }
		  Btype6size++;
		}
	    }
	  if(best2index>-1)
	    {
	      besttktkmass = best2index;
	      isbesttktkmass[best2index] = 1;
	    }
	  if(bestindex>-1)
	    {
	      bestchi2 = bestindex;
	      isbestchi2[bestindex] = 1;
	    }
	  nt5->Fill();
	  
	  Bsize=0;
	  best=-1;
	  bestindex=-1;
	  best2=10000.;
	  best2index=-1;
	  for (int j=0;j<BInfo_size;j++) 
	    {
	      if(BInfo_type[j]>7) continue;
	      if (isBchannel[BInfo_type[j]-1]!=1) continue;
	      
	      //muonIdPass{{
	      if (doMuonSelection)
		{
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])) continue;
		  if(!(MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] || MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(abs(MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=3. || abs(MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])>=30.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]<1.) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if(MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>1.8) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]])<6.) continue;
		  if((MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]+MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]])<6.) continue;
		  
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]&4096)) continue;
		  if(!(MuonInfo_muqual[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]&4096)) continue;
		}
	      b4Pout->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);
	      temy = b4Pout->Rapidity();
	      if(REAL)
		{
		  if(!(((EvtInfo_RunNo>=210498&&EvtInfo_RunNo<=211256&&abs(temy)<2.4)||(EvtInfo_RunNo>=211313&&EvtInfo_RunNo<=211631&&abs(temy)<2.4))||(EvtInfo_RunNo>=211739&&EvtInfo_RunNo<=211831&&abs(temy)<2.4))) continue;
		}
	      else
		{
		  if((PbpMC==0)&&abs(temy)>=2.4) continue;
		  if((PbpMC==1)&&abs(temy)>=2.4) continue;
		}
	      if(BInfo_type[j]==7)
		{
		  fillTree(bP,bVtx,b4P,j,Btype7size,PION_MASS,PION_MASS,REAL,PbpMC);
		  if(chi2cl[Btype7size]>best)
		    {
		      best = chi2cl[Btype7size];
		      bestindex = Btype7size;
		    }
		  Btype7size++;
		}
	    }
	  if(size>0)
	    {
	      bestchi2 = bestindex;
	      isbestchi2[bestindex] = 1;
	    }
	  nt6->Fill();
	  
	  if(!REAL)
	    {
	      Gensize = 0;
	      for (int j=0;j<GenInfo_size;j++)
		{
		  bGen.SetPtEtaPhiM(GenInfo_pt[j],GenInfo_eta[j],GenInfo_phi[j],GenInfo_mass[j]);
		  flag=0;
		  for(Btype=1;Btype<8;Btype++)
		    {
		      if (signalGen(type,j))
			{
			  flag=type;
			  break;
			}
		    }
		  Genmu1pt[j] = -1;
		  Genmu1eta[j] = -20;
		  Genmu1phi[j] = -20;
		  Genmu1p[j] = -1;
		  Genmu2pt[j] = -1;
		  Genmu2eta[j] = -20;
		  Genmu2phi[j] = -20;
		  Genmu2p[j] = -1;
		  Gentk1pt[j] = -1;
		  Gentk1eta[j] = -20;
		  Gentk1phi[j] = -20;
		  Gentk2pt[j] = -1;
		  Gentk2eta[j] = -20;
		  Gentk2phi[j] = -20;
		  
		  if(flag!=0)
		    {
		      Genmu1pt[j] = GenInfo_pt[GenInfo_da1[GenInfo_da1[j]]];
		      Genmu1eta[j] = GenInfo_eta[GenInfo_da1[GenInfo_da1[j]]];
		      Genmu1phi[j] = GenInfo_phi[GenInfo_da1[GenInfo_da1[j]]];
		      Genmu1p[j] = Genmu1pt[j]*cosh(Genmu1eta[j]);
		      Genmu2pt[j] = GenInfo_pt[GenInfo_da2[GenInfo_da1[j]]];
		      Genmu2eta[j] = GenInfo_eta[GenInfo_da2[GenInfo_da1[j]]];
		      Genmu2phi[j] = GenInfo_phi[GenInfo_da2[GenInfo_da1[j]]];
		      Genmu2p[j] = Genmu2pt[j]*cosh(Genmu2eta[j]);
		      if(flag==1||flag==2)
			{
			  Gentk1pt[j] = GenInfo_pt[GenInfo_da2[j]];
			  Gentk1eta[j] = GenInfo_eta[GenInfo_da2[j]];
			  Gentk1phi[j] = GenInfo_phi[GenInfo_da2[j]];
			}
		      else
			{
			  Gentk1pt[j] = GenInfo_pt[GenInfo_da1[GenInfo_da2[j]]];
			  Gentk1eta[j] = GenInfo_eta[GenInfo_da1[GenInfo_da2[j]]];
			  Gentk1phi[j] = GenInfo_phi[GenInfo_da1[GenInfo_da2[j]]];
			  Gentk2pt[j] = GenInfo_pt[GenInfo_da2[GenInfo_da2[j]]];
			  Gentk2eta[j] = GenInfo_eta[GenInfo_da2[GenInfo_da2[j]]];
			  Gentk2phi[j] = GenInfo_phi[GenInfo_da2[GenInfo_da2[j]]];
			}
		    }
		  Gensize = GenInfo_size;
		  Geny[j] = bGen.Rapidity();
		  Geneta[j] = bGen.Eta();
		  Genphi[j] = bGen.Phi();
		  Genpt[j] = bGen.Pt();
		  GenpdgId[j] = GenInfo_pdgId[j];
		  GenisSignal[j] = flag;
		}
	      ntGen->Fill();
	      
	    }
	}//writeBmeson
 
     
    }
  
  outf->Write();
  outf->Close();
}


double findMass(int particlePdgId)
{
  if(TMath::Abs(particlePdgId)==211) return PION_MASS;
  if(TMath::Abs(particlePdgId)==321) return KAON_MASS;
}

void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize,  int REAL, int PbpMC)
{
  //EvtInfo
  RunNo = EvtInfo_RunNo+10*PbpMC;
  EvtNo = EvtInfo_EvtNo;
  Dsize = typesize+1;
  PVx = EvtInfo_PVx;
  PVy = EvtInfo_PVy;
  PVz = EvtInfo_PVz;
  PVxE = EvtInfo_PVxE;
  PVyE = EvtInfo_PVyE;
  PVzE = EvtInfo_PVzE;
  PVnchi2 = EvtInfo_PVnchi2;
  PVchi2 = EvtInfo_PVchi2;
  //DInfo
  bP->SetXYZ(DInfo_px[j],DInfo_py[j],DInfo_pz[j]*0);
  bVtx->SetXYZ(DInfo_vtxX[j]-EvtInfo_PVx,
	       DInfo_vtxY[j]-EvtInfo_PVy,
	       DInfo_vtxZ[j]*0-EvtInfo_PVz*0);
  b4P->SetXYZM(DInfo_px[j],DInfo_py[j],DInfo_pz[j],DInfo_mass[j]);
  Dindex[typesize] = typesize;
  DisGoodCand[typesize] = DInfo_isGoodCand[j];
  Dmass[typesize] = DInfo_mass[j];
  Dpt[typesize] = DInfo_pt[j];
  Deta[typesize] = DInfo_eta[j];
  Dphi[typesize] = DInfo_phi[j];
  Dy[typesize] = b4P->Rapidity();
  DvtxX[typesize] = DInfo_vtxX[j] - EvtInfo_PVx;
  DvtxY[typesize] = DInfo_vtxY[j] - EvtInfo_PVy;
  Dd0[typesize] = TMath::Sqrt((DInfo_vtxX[j]-EvtInfo_PVx)*(DInfo_vtxX[j]-EvtInfo_PVx)+(DInfo_vtxY[j]-EvtInfo_PVy)*(DInfo_vtxY[j]-EvtInfo_PVy));
  Dd0Err[typesize] = TMath::Sqrt(DInfo_vtxXE[j]*DInfo_vtxXE[j]+DInfo_vtxYE[j]*DInfo_vtxYE[j]);
  Dchi2ndf[typesize] = DInfo_vtxchi2[j]/DInfo_vtxdof[j];
  Dchi2cl[typesize] = TMath::Prob(DInfo_vtxchi2[j],DInfo_vtxdof[j]);
  Ddtheta[typesize] = bP->Angle(*bVtx);
  Dlxy[typesize] = ((DInfo_vtxX[j]-EvtInfo_PVx)*DInfo_px[j] + (DInfo_vtxY[j]-EvtInfo_PVy)*DInfo_py[j])/DInfo_pt[j];
  Disbestchi2[typesize] = false;
  //DInfo.b4fitInfo
  Db4fit_mass[typesize] = DInfo_b4fit_mass[j];
  Db4fit_pt[typesize] = DInfo_b4fit_pt[j];
  Db4fit_eta[typesize] = DInfo_b4fit_eta[j];
  Db4fit_phi[typesize] = DInfo_b4fit_phi[j];
  //DInfo.trkInfo
  Double_t trk1mass,trk2mass,trk3mass,trk4mass;
  Dtrk1Pt[typesize] = TrackInfo_pt[DInfo_rftk1_index[j]];
  Dtrk2Pt[typesize] = TrackInfo_pt[DInfo_rftk2_index[j]];
  Dtrk1Eta[typesize] = TrackInfo_eta[DInfo_rftk1_index[j]];
  Dtrk2Eta[typesize] = TrackInfo_eta[DInfo_rftk2_index[j]];
  Dtrk1Phi[typesize] = TrackInfo_phi[DInfo_rftk1_index[j]];
  Dtrk2Phi[typesize] = TrackInfo_phi[DInfo_rftk2_index[j]];
  trk1mass = findMass(DInfo_rftk1_MassHypo[j]);
  trk2mass = findMass(DInfo_rftk2_MassHypo[j]);
  b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk1_index[j]],TrackInfo_eta[DInfo_rftk1_index[j]],TrackInfo_phi[DInfo_rftk1_index[j]],trk1mass);
  Dtrk1Y[typesize] = b4P->Rapidity();
  b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk2_index[j]],TrackInfo_eta[DInfo_rftk2_index[j]],TrackInfo_phi[DInfo_rftk2_index[j]],trk2mass);
  Dtrk2Y[typesize] = b4P->Rapidity();
  Dtrk1Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk1_index[j]];
  Dtrk2Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk2_index[j]];
  Dtrk1D0Err[typesize] = TrackInfo_d0error[DInfo_rftk1_index[j]];
  Dtrk2D0Err[typesize] = TrackInfo_d0error[DInfo_rftk2_index[j]];
  Dtrk1PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk1_index[j]];
  Dtrk2PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk2_index[j]];
  Dtrk1StripHit[typesize] = TrackInfo.striphit[DInfo_rftk1_index[j]];
  Dtrk2StripHit[typesize] = TrackInfo.striphit[DInfo_rftk2_index[j]];
  Dtrk1Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk1_index[j]]/TrackInfo_ndf[DInfo_rftk1_index[j]];
  Dtrk2Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk2_index[j]]/TrackInfo_ndf[DInfo_rftk2_index[j]];
  Dtrk1MassHypo[typesize] = DInfo_rftk1_MassHypo[j];
  Dtrk2MassHypo[typesize] = DInfo_rftk2_MassHypo[j];
  if(DInfo_type[j]==1||DInfo_type[j]==2)
    {
      Dtrk3Pt[typesize] = -1;
      Dtrk4Pt[typesize] = -1;
      Dtrk3Eta[typesize] = -20;
      Dtrk4Eta[typesize] = -20;
      Dtrk3Phi[typesize] = -20;
      Dtrk4Phi[typesize] = -20;
      Dtrk3Y[typesize] = -1;
      Dtrk4Y[typesize] = -1;
      Dtrk3Dxy[typesize] = -1;
      Dtrk4Dxy[typesize] = -1;
      Dtrk3D0Err[typesize] = -1;
      Dtrk4D0Err[typesize] = -1;
      Dtrk3PixelHit[typesize] = -1;
      Dtrk4PixelHit[typesize] = -1;
      Dtrk3StripHit[typesize] = -1;
      Dtrk4StripHit[typesize] = -1;
      Dtrk3Chi2ndf[typesize] = -1;
      Dtrk4Chi2ndf[typesize] = -1;
      Dtrk3MassHypo[typesize] = 0;
      Dtrk4MassHypo[typesize] = 0;
      DtktkResmass[typesize] = -1;
      DtktkResvProb[typesize] = -1;
      DtktkRespt[typesize] = -1;
      DtktkReseta[typesize] = -20;
      DtktkResphi[typesize] = -20;
      DtktkResy[typesize] = -1;
    }
  else if(DInfo_type[j]==3||DInfo_type[j]==4)
    {
      Dtrk3Pt[typesize] = TrackInfo_pt[DInfo_rftk3_index[j]];
      Dtrk3Eta[typesize] = TrackInfo_eta[DInfo_rftk3_index[j]];
      Dtrk3Phi[typesize] = TrackInfo_phi[DInfo_rftk3_index[j]];
      trk3mass = findMass(DInfo_rftk3_MassHypo[j]);
      b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk3_index[j]],TrackInfo_eta[DInfo_rftk3_index[j]],TrackInfo_phi[DInfo_rftk3_index[j]],trk3mass);
      Dtrk3Y[typesize] = b4P->Rapidity();
      Dtrk3Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk3_index[j]];
      Dtrk3D0Err[typesize] = TrackInfo_d0error[DInfo_rftk3_index[j]];
      Dtrk3PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk3_index[j]];
      Dtrk3StripHit[typesize] = TrackInfo.striphit[DInfo_rftk3_index[j]];
      Dtrk3Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk3_index[j]]/TrackInfo_ndf[DInfo_rftk3_index[j]];
      Dtrk3MassHypo[typesize] = DInfo_rftk3_MassHypo[j];
      Dtrk4Pt[typesize] = -1;
      Dtrk4Eta[typesize] = -20;
      Dtrk4Phi[typesize] = -20;
      Dtrk4Y[typesize] = -1;
      Dtrk4Dxy[typesize] = -1;
      Dtrk4D0Err[typesize] = -1;
      Dtrk4PixelHit[typesize] = -1;
      Dtrk4StripHit[typesize] = -1;
      Dtrk4Chi2ndf[typesize] = -1;
      Dtrk4MassHypo[typesize] = 0;
      DtktkResmass[typesize] = -1;
      DtktkResvProb[typesize] = -1;
      DtktkRespt[typesize] = -1;
      DtktkReseta[typesize] = -20;
      DtktkResphi[typesize] = -20;
      DtktkResy[typesize] = -1;
    }
  else if(DInfo_type[j]==5)
    {
      Dtrk3Pt[typesize] = TrackInfo_pt[DInfo_rftk3_index[j]];
      Dtrk4Pt[typesize] = TrackInfo_pt[DInfo_rftk4_index[j]];
      Dtrk3Eta[typesize] = TrackInfo_eta[DInfo_rftk3_index[j]];
      Dtrk4Eta[typesize] = TrackInfo_eta[DInfo_rftk4_index[j]];
      Dtrk3Phi[typesize] = TrackInfo_phi[DInfo_rftk3_index[j]];
      Dtrk4Phi[typesize] = TrackInfo_phi[DInfo_rftk4_index[j]];
      trk3mass = findMass(DInfo_rftk3_MassHypo[j]);
      trk4mass = findMass(DInfo_rftk4_MassHypo[j]);
      b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk3_index[j]],TrackInfo_eta[DInfo_rftk3_index[j]],TrackInfo_phi[DInfo_rftk3_index[j]],trk3mass);
      Dtrk3Y[typesize] = b4P->Rapidity();
      b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk4_index[j]],TrackInfo_eta[DInfo_rftk4_index[j]],TrackInfo_phi[DInfo_rftk4_index[j]],trk4mass);
      Dtrk4Y[typesize] = b4P->Rapidity();
      Dtrk3Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk3_index[j]];
      Dtrk4Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk4_index[j]];
      Dtrk3D0Err[typesize] = TrackInfo_d0error[DInfo_rftk3_index[j]];
      Dtrk4D0Err[typesize] = TrackInfo_d0error[DInfo_rftk4_index[j]];
      Dtrk3PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk3_index[j]];
      Dtrk4PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk4_index[j]];
      Dtrk3StripHit[typesize] = TrackInfo.striphit[DInfo_rftk3_index[j]];
      Dtrk4StripHit[typesize] = TrackInfo.striphit[DInfo_rftk4_index[j]];
      Dtrk3Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk3_index[j]]/TrackInfo_ndf[DInfo_rftk3_index[j]];
      Dtrk4Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk4_index[j]]/TrackInfo_ndf[DInfo_rftk4_index[j]];
      Dtrk3MassHypo[typesize] = DInfo_rftk3_MassHypo[j];
      Dtrk4MassHypo[typesize] = DInfo_rftk4_MassHypo[j];
      DtktkResmass[typesize] = -1;
      DtktkResvProb = -1;
      DtktkRespt = -1;
      DtktkReseta = -20;
      DtktkResphi = -20;
      DtktkResy = -1;
    }
  else if(DInfo_type[j]==6||DInfo_type[j]==7)
    {
      Dtrk3Pt[typesize] = TrackInfo_pt[DInfo_rftk3_index[j]];
      Dtrk4Pt[typesize] = TrackInfo_pt[DInfo_rftk4_index[j]];
      Dtrk3Eta[typesize] = TrackInfo_eta[DInfo_rftk3_index[j]];
      Dtrk4Eta[typesize] = TrackInfo_eta[DInfo_rftk4_index[j]];
      Dtrk3Phi[typesize] = TrackInfo_phi[DInfo_rftk3_index[j]];
      Dtrk4Phi[typesize] = TrackInfo_phi[DInfo_rftk4_index[j]];
      trk3mass = findMass(DInfo_rftk3_MassHypo[j]);
      trk4mass = findMass(DInfo_rftk4_MassHypo[j]);
      b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk3_index[j]],TrackInfo_eta[DInfo_rftk3_index[j]],TrackInfo_phi[DInfo_rftk3_index[j]],trk3mass);
      Dtrk3Y[typesize] = b4P->Rapidity();
      b4P->SetPtEtaPhiM(TrackInfo_pt[DInfo_rftk4_index[j]],TrackInfo_eta[DInfo_rftk4_index[j]],TrackInfo_phi[DInfo_rftk4_index[j]],trk4mass);
      Dtrk4Y[typesize] = b4P->Rapidity();
      Dtrk3Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk3_index[j]];
      Dtrk4Dxy[typesize] = TrackInfo_dxyPV[DInfo_rftk4_index[j]];
      Dtrk3D0Err[typesize] = TrackInfo_d0error[DInfo_rftk3_index[j]];
      Dtrk4D0Err[typesize] = TrackInfo_d0error[DInfo_rftk4_index[j]];
      Dtrk3PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk3_index[j]];
      Dtrk4PixelHit[typesize] = TrackInfo_pixelhit[DInfo_rftk4_index[j]];
      Dtrk3StripHit[typesize] = TrackInfo.striphit[DInfo_rftk3_index[j]];
      Dtrk4StripHit[typesize] = TrackInfo.striphit[DInfo_rftk4_index[j]];
      Dtrk3Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk3_index[j]]/TrackInfo_ndf[DInfo_rftk3_index[j]];
      Dtrk4Chi2ndf[typesize] = TrackInfo_chi2[DInfo_rftk4_index[j]]/TrackInfo_ndf[DInfo_rftk4_index[j]];
      Dtrk3MassHypo[typesize] = DInfo_rftk3_MassHypo[j];
      Dtrk4MassHypo[typesize] = DInfo_rftk4_MassHypo[j];
      DtktkResmass[typesize] = DInfo_tktkRes_mass[j];
      DtktkResvProb[typesize] = DInfo_tktkRes_vProb[j];
      DtktkRespt[typesize] = DInfo_tktkRes_pt[j];
      DtktkReseta[typesize] = DInfo_tktkRes_eta[j];
      DtktkResphi[typesize] = DInfo_tktkRes_phi[j];
      DtktkResy[typesize] = DInfo_tktkRes_y[j];
    }

  Int_t hypo=-1,DpdgId;
  if(DInfo_type[j]==1||DInfo_type[j]==2||DInfo_type[j]==5) DpdgId=DZERO_PDGID;
  else if(DInfo_type[j]==3||DInfo_type[j]==4) DpdgId=DPLUS_PDGID;
  else if(DInfo_type[j]==6||DInfo_type[j]==7) DpdgId=DSUBS_PDGID;
  Dgen[typesize] = 0;//gen init
  DgenIndex[typesize] = -1;//gen init
  Dgenpt[typesize] = -1;
  Dgeneta[typesize] = -20;
  Dgenphi[typesize] = -20;
  Dgeny[typesize] = -1;
  Int_t rGenIdxTk1 = -1;
  Int_t rGenIdxTk2 = -1;
  Int_t dGenIdxTk1 = -1;
  Int_t dGenIdxTk2 = -1;
  Int_t dGenIdxTk3 = -1;
  Int_t dGenIdxTk4 = -1;
  Int_t dGenIdxRes = -1;
  if(!real)
    {
      if(DInfo_type[j]==1||DInfo_type[j]==2||DInfo_type[j]==3||DInfo_type[j]==4||DInfo_type[j]==5)
	{
	  if(TrackInfo_geninfo_index[DInfo_rftk1_index[j]]>-1)
	    {
	      int level=0;
	      if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==DInfo_rftk1_MassHypo[j])
		{
		  hypo=-1;
		  if(TMath::Abs(DInfo_rftk1_MassHypo[j])==KAON_PDGID) hypo=0;
		  else if(TMath::Abs(DInfo_rftk1_MassHypo[j])==PION_PDGID) hypo=1;
		  if(hypo==0||hypo==1)
		    {
		      level=1;
		      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]>-1)
			{
			  if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]])==DpdgId)
			    {
			      level=3;
			      dGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]];
			    }
			}
		    }
		}
	    }
	  Dgen[typesize]=level;
	  if(TrackInfo_geninfo_index[DInfo_rftk2_index[j]]>-1)
	    {
	      int level =0;
	      if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==DInfo_rftk2_MassHypo[j])
		{
		  if((hypo==0&&TMath::Abs(DInfo_rftk2_MassHypo[j])==PION_PDGID)||(hypo==1&&TMath::Abs(DInfo_rftk2_MassHypo[j])==KAON_PDGID))
		    {
		      if(hypo==1&&TMath::Abs(DInfo_rftk2_MassHypo[j])==KAON_PDGID) hypo=2;
		      level=1;
		      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]>-1)
                        {
                          if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]])==DpdgId)
                            {
                              level=3;
                              dGenIdxTk2=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]];
                            }
                        }
		    }
		  else if(hypo==1&&TMath::Abs(DInfo_rftk2_MassHypo[j])==PION_PDGID&&(DInfo_type[j]==3||DInfo_type[j]==4||DInfo_type[j]==5))
		    {
		      hypo=3;
		      level=1;
		      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]>-1)
			{
			  if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]])==DpdgId)
                            {
			      level=3;
                              dGenIdxTk2=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]];
			    }
			}
		    }
		}
	      Dgen[typesize]+=(level*10);
	    }
	  if(DInfo_type[j]==1||DInfo_type[j]==2)
	    {
	      Dgen[typesize]+=3300;
	    }
	  else
	    {
	      if(TrackInfo_geninfo_index[DInfo_rftk3_index[j]]>-1)
		{
		  int level=0;
		  if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]])==DInfo_rftk3_MassHypo[j])
		    {
		      if(((hypo==0||hypo==2)&&TMath::Abs(DInfo_rftk3_MassHypo[j])==PION_PDGID)||(hypo==3&&TMath::Abs(DInfo_rftk3_MassHypo[j])==KAON_PDGID))
			{
			  if(hypo==3&&TMath::Abs(DInfo_rftk3_MassHypo[j])==KAON_PDGID) hypo=4;
			  level=1;
			  if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]]>-1)
			    {
			      if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]]])==DpdgId)
				{
				  level=3;
				  dGenIdxTk3=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]];
				}
			    }
			  else if(hypo==3&&&&TMath::Abs(DInfo_rftk3_MassHypo[j])==PION_PDGID&&DInfo_type[j]==5)
			    {
			      hypo=5;
			      level=1;
			      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]]>-1)
				{
				  if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]]])==DpdgId)
				    {
				      level=3;
				      dGenIdxTk3=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk3_index[j]]];
				    }
				}
			    }
			}
		    }
		  Dgen[typesize]+=(level*100);
		}
	      if(DInfo_type[j]==3||DInfo_type[j]==4)
		{
		  Dgen[typesize]+=3000;
		}
	      else
		{
		  if(TrackInfo_geninfo_index[DInfo_rftk4_index[j]]>-1)
		    {
		      int level=0;
		      if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk4_index[j]]])==DInfo_rftk4_MassHypo[j])
			{
			  if(((hypo==0||hypo==2||hypo==4)&&TMath::Abs(DInfo_rftk4_MassHypo[j])==PION_PDGID)||(hypo==5&&TMath::Abs(DInfo_rftk4_MassHypo[j])==KAON_PDGID))
			    {
			      level=1;
			      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk4_index[j]]]>-1)
				{
				  if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk4_index[j]]]])==DpdgId)
				    {
				      level=3;
				      dGenIdxTk4=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk4_index[j]]];
				    }
				}
			    }
			}
		      Dgen[typesize]+=(level*1000);
		    }
		}
	    }
	  if(Dgen[typesize]==3000&&dGenIdxTk1==dGenIdxTk2)
	    {
	      if(DInfo_type[j]==1||DInfo_type[j]==2)
		{
		  Dgen[typesize]+=20000;
		}
	      else if(dGenIdxTk1==dGenIdxTk3)
		{
		  if(DInfo_type[j]==3||DInfo_type[j]==4)
		    {
		      Dgen[typesize]+=20000;
		    }
		  else if(dGenIdxTk1==dGenIdxTk4)
		    {
		      Dgen[typesize]+=20000;
		    }
		}	    
	    }
	}//if(DInfo_type[j]==1||DInfo_type[j]==2||DInfo_type[j]==3||DInfo_type[j]==4||DInfo_type[j]==5)
      if(Dgen==23333)
	{
	  DgenIndex[typesize] = dGenIdxTk1;
	  Dgenpt[typesize] = GenInfo_pt[DgenIndex[typesize]];
	  Dgeneta[typesize] = GenInfo_eta[DgenIndex[typesize]];
	  Dgenphi[typesize] = GenInfo_phi[DgenIndex[typesize]];
	  b4P->SetXYZM(GenInfo_pt[DgenIndex[typesize]]*cos(GenInfo_phi[DgenIndex[typesize]]),
		       GenInfo_pt[DgenIndex[typesize]]*sin(GenInfo_phi[DgenIndex[typesize]]),
		       GenInfo_pt[DgenIndex[typesize]]*sinh(GenInfo_eta[DgenIndex[typesize]]),
		       GenInfo_mass[DgenIndex[typesize]]);
	  Dgeny[typesize] = b4P->Rapidity();
	}
    }//if(!real)
}//fillDtree




void fillBTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, int j, int typesize, float track_mass1, float track_mass2, int REAL, int PbpMC)
{
  //Event Info
  RunNo = EvtInfo_RunNo+10*PbpMC;
  EvtNo = EvtInfo_EvtNo;
  Bsize = typesize+1;
  PVx = EvtInfo_PVx;
  PVy = EvtInfo_PVy;
  PVz = EvtInfo_PVz;
  PVxE = EvtInfo_PVxE;
  PVyE = EvtInfo_PVyE;
  PVzE = EvtInfo_PVzE;
  PVnchi2 = EvtInfo_PVnchi2;
  PVchi2 = EvtInfo_PVchi2;

  //BInfo
  bP->SetXYZ(BInfo_px[j],BInfo_py[j],BInfo_pz[j]*0);
  bVtx->SetXYZ(BInfo_vtxX[j]-EvtInfo_PVx,
	       BInfo_vtxY[j]-EvtInfo_PVy,
	       BInfo_vtxZ[j]*0-EvtInfo_PVz*0);
  b4P->SetXYZM(BInfo_px[j],BInfo_py[j],BInfo_pz[j],BInfo_mass[j]);

  Bindex[typesize] = typesize;
  Bmass[typesize] = BInfo_mass[j];
  Bpt[typesize] = BInfo_pt[j];
  Beta[typesize] = BInfo_eta[j];
  Bphi[typesize] = BInfo_phi[j];
  By[typesize] = b4P->Rapidity();
  BvtxX[typesize] = BInfo_vtxX[j]-EvtInfo_PVx;
  BvtxY[typesize] = BInfo_vtxY[j]-EvtInfo_PVy;
  Bd0[typesize] = TMath::Sqrt((BInfo_vtxX[j]-EvtInfo_PVx)*(BInfo_vtxX[j]-EvtInfo_PVx)+(BInfo_vtxY[j]-EvtInfo_PVy)*(BInfo_vtxY[j]-EvtInfo_PVy));
  Bd0Err[typesize] = TMath::Sqrt(BInfo_vtxXE[j]*BInfo_vtxXE[j]+BInfo_vtxYE[j]*BInfo_vtxYE[j]);
  Bchi2ndf[typesize] = BInfo_vtxchi2[j]/BInfo_vtxdof[j];
  Bchi2cl[typesize] = TMath::Prob(BInfo_vtxchi2[j],BInfo_vtxdof[j]);
  Bdtheta[typesize] = bP->Angle(*bVtx);
  Blxy[typesize] = ((BInfo_vtxX[j]-EvtInfo_PVx)*BInfo_px[j]+(BInfo_vtxY[j]-EvtInfo_PVy)*BInfo_py[j])/BInfo_pt[j];
  Bisbestchi2[typesize] = false;  
  //BInfo_MuInfo
  Bmu1TrackerMuArbitrated[typesize] = MuonInfo_TrackerMuonArbitrated[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2TrackerMuArbitrated[typesize] = MuonInfo_TrackerMuonArbitrated[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1isTrackerMuon[typesize] = MuonInfo_isTrackerMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2isTrackerMuon[typesize] = MuonInfo_isTrackerMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1isGlobalMuon[typesize] = MuonInfo_isGlobalMuon[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2isGlobalMuon[typesize] = MuonInfo_isGlobalMuon[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1TMOneStationTight[typesize] = MuonInfo_TMOneStationTight[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2TMOneStationTight[typesize] = MuonInfo_TMOneStationTight[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1InPixelLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2InPixelLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1InTrackerLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]] + MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2InTrackerLayer[typesize] = MuonInfo_i_nPixelLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]] + MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1InStripLayer[typesize] = MuonInfo_i_nStripLayer[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2InStripLayer[typesize] = MuonInfo_i_nStripLayer[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Double_t mu1px,mu1py,mu1pz,mu1E;
  mu1px = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*cos(MuonInfo_phi[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  mu1py = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*sin(MuonInfo_phi[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  mu1pz = MuonInfo_pt[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]*sinh(MuonInfo_eta[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]);
  b4P->SetXYZM(mu1px,mu1py,mu1pz,MUON_MASS);
  Bmu1eta[typesize] = b4P->Eta();
  Bmu1phi[typesize] = b4P->Phi();
  Bmu1y[typesize] = b4P->Rapidity();
  Bmu1pt[typesize] = b4P->Pt();
  Bmu1p[typesize] = b4P->P();
  Bmu1E[typesize] = b4P->E();
  Double_t mu2px,mu2py,mu2pz,mu2E;  
  mu2px = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*cos(MuonInfo_phi[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  mu2py = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*sin(MuonInfo_phi[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  mu2pz = MuonInfo_pt[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]*sinh(MuonInfo_eta[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]);
  b4P->SetXYZM(mu2px,mu2py,mu2pz,MUON_MASS);
  Bmu2eta[typesize] = b4P->Eta();
  Bmu2phi[typesize] = b4P->Phi();
  Bmu2y[typesize] = b4P->Rapidity();
  Bmu2pt[typesize] = b4P->Pt();
  Bmu2E[typesize] = b4P->E();
  Bmu1dzPV[typesize] = MuonInfo_dzPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2dzPV[typesize] = MuonInfo_dzPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1dxyPV[typesize] = MuonInfo_dxyPV[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2dxyPV[typesize] = MuonInfo_dxyPV[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1normchi2[typesize] = MuonInfo_normchi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2normchi2[typesize] = MuonInfo_normchi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  Bmu1Chi2ndf[typesize] = MuonInfo_i_chi2[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]/MuonInfo_i_ndf[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]];
  Bmu2Chi2ndf[typesize] = MuonInfo_i_chi2[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]/MuonInfo_i_ndf[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]];
  //BInfo_ujInfo
  b4P->SetPxPyPzE(mu1px+mu2px,
		  mu1py+mu2py,
		  mu1pz+mu2pz,
		  mu1E+mu2E);
  Bmumumass[typesize] = b4P->Mag();
  Bmumueta[typesize] = b4P->Eta();
  Bmumuphi[typesize] = b4P->Phi();
  Bmumuy[typesize] = b4P->Rapidity();
  Bmumupt[typesize] = b4P->Pt();
  Bujmass[typesize] = BInfo_uj_mass[BInfo_rfuj_index[j]];
  BujvProb[typesize] = TMath::Prob(BInfo_uj_vtxchi2[BInfo_rfuj_index[j]],BInfo_uj_vtxdof[BInfo_rfuj_index[j]]);
  b4P->SetXYZM(BInfo_uj_px[BInfo_rfuj_index[j]],
	       BInfo_uj_py[BInfo_rfuj_index[j]],
	       BInfo_uj_pz[BInfo_rfuj_index[j]],
	       BInfo_uj_mass[BInfo_rfuj_index[j]]);
  Bujpt[typesize] = b4P->Pt();
  Bujeta[typesize] = b4P->PseudoRapidity();
  Bujphi[typesize] = b4P->Phi();
  Bujy[typesize] = b4P->Rapidity();
  Bujlxy[typesize] = ((BInfo_uj_vtxX[BInfo_rfuj_index[j]]-EvtInfo_PVx)*BInfo_uj_px[BInfo_rfuj_index[j]]+(BInfo_uj_vtxY[BInfo_rfuj_index[j]]-EvtInfo_PVy)*BInfo_uj_py[BInfo_rfuj_index[j]])/ujpt[typesize];
  //BInfo_trkInfo
  Double_t tk1px,tk1py,tk1pz,tk1E;
  Double_t tk2px,tk2py,tk2pz,tk2E;
  if(BInfo_type[j]==1||BInfo_type[j]==2)
    {
      b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass1);
      Btrk1Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
      Btrk1Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
      Btrk1Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
      Btrk1Y[typesize] = b4P->Rapidity();
      Btrk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
      Btrk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
      Btrk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
      Btrk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
      Btrk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
      Btrk2Pt[typesize] = -1;
      Btrk2Eta[typesize] = -20;
      Btrk2Phi[typesize] = -20;
      Btrk2Y[typesize] = -1;
      Btrk2Dxy[typesize] = -1;
      Btrk2D0Err[typesize] = -1;
      Btrk2PixelHit[typesize] = -1;
      Btrk2StripHit[typesize] = -1;
      Btrk2Chi2ndf[typesize] = -1;
      Btktkmass[typesize] = -1;
      BtktkvProb[typesize] = -1;
      Btktkpt[typesize] = -1;
      Btktketa[typesize] = -20;
      Btktkphi[typesize] = -20;
      Btktky[typesize] = -1;
      Bdoubletmass[typesize] = -1;
      Bdoubletpt[typesize] = -1;
      Bdoubleteta[typesize] = -20;
      Bdoubletphi[typesize] = -20;
      Bdoublety[typesize] = -1;
    }
  else
    {
      if(BInfo_type[j]==5)
	{
	  b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],track_mass1);
	  Btrk1Pt[typesize] = TrackInfo_pt[BInfo_rftk2_index[j]];
	  Btrk1Eta[typesize] = TrackInfo_eta[BInfo_rftk2_index[j]];
	  Btrk1Phi[typesize] = TrackInfo_phi[BInfo_rftk2_index[j]];
	  Btrk1Y[typesize] = b4P->Rapidity();
	  Btrk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
	  Btrk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk2_index[j]];
	  Btrk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
	  Btrk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk2_index[j]];
	  Btrk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
	  tk1px = b4P->Px();
	  tk1py = b4P->Py();
	  tk1pz = b4P->Pz();
	  tk1E = b4P->E();
	  b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass2);
	  Btrk2Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
	  Btrk2Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
	  Btrk2Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
	  Btrk2Y[typesize] = b4P->Rapidity();
	  Btrk2Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
	  Btrk2D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
	  Btrk2PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
	  Btrk2StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
	  Btrk2Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
	  tk2px = b4P->Px();
	  tk2py = b4P->Py();
	  tk2pz = b4P->Pz();
	  tk2E = b4P->E();
	}
      else
	{
	  b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk1_index[j]],TrackInfo_eta[BInfo_rftk1_index[j]],TrackInfo_phi[BInfo_rftk1_index[j]],track_mass1);
	  Btrk1Pt[typesize] = TrackInfo_pt[BInfo_rftk1_index[j]];
	  Btrk1Eta[typesize] = TrackInfo_eta[BInfo_rftk1_index[j]];
	  Btrk1Phi[typesize] = TrackInfo_phi[BInfo_rftk1_index[j]];
	  Btrk1Y[typesize] = b4P->Rapidity();
	  Btrk1Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk1_index[j]];
	  Btrk1D0Err[typesize] = TrackInfo_d0error[BInfo_rftk1_index[j]];
	  Btrk1PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk1_index[j]];
	  Btrk1StripHit[typesize] = TrackInfo_striphit[BInfo_rftk1_index[j]];
	  Btrk1Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk1_index[j]]/TrackInfo_ndf[BInfo_rftk1_index[j]];
	  tk1px = b4P->Px();
	  tk1py = b4P->Py();
	  tk1pz = b4P->Pz();
	  tk1E = b4P->E();
	  b4P->SetPtEtaPhiM(TrackInfo_pt[BInfo_rftk2_index[j]],TrackInfo_eta[BInfo_rftk2_index[j]],TrackInfo_phi[BInfo_rftk2_index[j]],track_mass2);
	  Btrk2Pt[typesize] = TrackInfo_pt[BInfo_rftk2_index[j]];
	  Btrk2Eta[typesize] = TrackInfo_eta[BInfo_rftk2_index[j]];
	  Btrk2Phi[typesize] = TrackInfo_phi[BInfo_rftk2_index[j]];
	  Btrk2Y[typesize] = b4P->Rapidity();
	  Btrk2Dxy[typesize] = TrackInfo_dxyPV[BInfo_rftk2_index[j]];
	  Btrk2D0Err[typesize] = TrackInfo_d0error[BInfo_rftk2_index[j]];
	  Btrk2PixelHit[typesize] = TrackInfo_pixelhit[BInfo_rftk2_index[j]];
	  Btrk2StripHit[typesize] = TrackInfo_striphit[BInfo_rftk2_index[j]];
	  Btrk2Chi2ndf[typesize] = TrackInfo_chi2[BInfo_rftk2_index[j]]/TrackInfo_ndf[BInfo_rftk2_index[j]];
	  Btk2px = b4P->Px();
	  Btk2py = b4P->Py();
	  Btk2pz = b4P->Pz();
	  Btk2E = b4P->E();
	}
      b4P->SetPxPyPzE(tk1px+tk2px,
		      tk1py+tk2py,
		      tk1pz+tk2pz,
		      tk1E+tk2E);
      Btktkmass[typesize] = b4P->Mag();
      Btktketa[typesize] = b4P->Eta();
      Btktkphi[typesize] = b4P->Phi();
      Btktky[typesize] = b4P->Rapidity();
      Btktkpt[typesize] = b4P->Pt();
      BtktkvProb[typesize] = TMath::Prob(BInfo_tktk_vtxchi2[j],BInfo_tktk_vtxdof[j]);
      Bdoubletmass[typesize] = BInfo_tktk_mass[j];
      b4P->SetXYZM(BInfo_tktk_px[j],BInfo_tktk_py[j],BInfo_tktk_pz[j],BInfo_tktk_mass[j]);
      Bdoubletpt[typesize] = b4P->Pt();
      Bdoubleteta[typesize] = b4P->PseudoRapidity();
      Bdoubletphi[typesize] = b4P->Phi();
      Bdoublety[typesize] = b4P->Rapidity();
    }

  //gen info judgement
  if(!REAL)
    {
      Bgen[typesize] = 0;//gen init
      BgenIndex[typesize] = -1;//gen init
      Bgenpt[typesize] = -1;
      Bgeneta[typesize] = -20;
      Bgenphi[typesize] = -20;
      Bgeny[typesize] = -1;     
      int mGenIdxTk1=-1;
      int mGenIdxTk2=-1;
      int bGenIdxTk1=-1;
      int bGenIdxTk2=-1;
      int bGenIdxMu1=-1;
      int bGenIdxMu2=-1;
      int ujGenIdxMu1=-1;
      int ujGenIdxMu2=-1;
      
      float BId,MId,tk1Id,tk2Id;
      //tk1:positive, tk2:negtive
      if(BInfo_type[j]==1)
	{
	  BId = 521;//B+-
	  MId = -1;
	  tk1Id = 321;//K+-
	  tk2Id = -1;
	}
      if(BInfo_type[j]==2)
	{
	  BId = 521;//B+-
	  MId = -1;
	  tk1Id = 211;//pi+-
	  tk2Id = -1;
	}
      if(BInfo_type[j]==3)
	{
	  BId = 511;//B0
	  MId = 310;//Ks
	  tk1Id = 211;//pi+
	  tk2Id = 211;//pi-
	}
      if(BInfo_type[j]==4)
	{
	  BId = 511;//B0
	  MId = 313;//K*0
	  tk1Id = 321;//K+
	  tk2Id = 211;//pi-
	}
      if(BInfo_type[j]==5)
	{
	  BId = 511;//B0
	  MId = 313;//K*0
	  tk1Id = 211;//pi+
	  tk2Id = 321;//K-
	}
      if(BInfo_type[j]==6)
	{
	  BId = 531;//Bs
	  MId = 333;//phi
	  tk1Id = 321;//K+
	  tk2Id = 321;//K-
	}      
      int twoTks,kStar,flagkstar=0;
      if(BInfo_type[j]==1 || BInfo_type[j]==2) twoTks=0;
      else twoTks=1;
      if(BInfo_type[j]==4 || BInfo_type[j]==5) kStar=1;
      else kStar=0;
      int nonprompt=0,prompt=0;

      // tk1
      if(TrackInfo_geninfo_index[BInfo_rftk1_index[j]]>-1)
	{
	  int level =0;
	  if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]])==tk1Id)
	    {
	      level = 1;
	      if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]>-1)
		{
		  if(!twoTks)//one trk channel
		    {
		      mGenIdxTk1=0;
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==BId)
			{
			  level = 3;
			  bGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
			}		  
		    }
		  else//two trk channel
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==MId)
			{
			  level = 2;
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]])==BId)
				{
				  level = 3;
				  bGenIdxTk1=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]];
				}
			    }
			  mGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
			}
		    }
		}
	    }
	  Bgen[typesize]=level;
	}
      
      //tk2
      if(!twoTks)//one trk channel
	{
	  Bgen[typesize]+=30;
	  mGenIdxTk2=0;
	  bGenIdxTk2=0;
	}
      else//two trk channel
	{
	  if(TrackInfo_geninfo_index[BInfo_rftk2_index[j]]>-1)
	    {
	      int level =0;
	      if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]])==tk2Id)
		{
		  level = 1;
		  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]>-1)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]])==MId)
			{
			  level = 2;
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]])==BId)
				{
				  level = 3;
				  bGenIdxTk2 = GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]];
				}
			    }
			  mGenIdxTk2 = GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]];
			}
		    }
		}
	      Bgen[typesize]+=(level*10);
	    }
	}      
      //mu1
      //cout<<MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]<<endl;
      if(MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]>-1)
	{
	  int level =0;
	  if(abs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]])==13)
	    {
	      level=1;
	      if(GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]>-1)
		{
		  if(GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]==443)
		    {
		      ujGenIdxMu1 = GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]];
		      level=2;
		      if(GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]>-1)
			{
			  if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]]])==BId)
			    {
			      //nonprompt=1;
			      level = 3;
			      bGenIdxMu1=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu1_index[BInfo_rfuj_index[j]]]]];
			      flagkstar++;///////////////////////////////////////////////=1
			    }
			}
		      else 
			{
			  prompt=1;
			}
		    } 
		}
	    }
	  Bgen[typesize]+=(level*100);
	}
      
      //mu2
      if(MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]>-1)
	{  
	  int level =0;
	  if(abs(GenInfo_pdgId[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]])==13)
	    {
	      level = 1;
	      if(GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]>-1)
		{
		  if(GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]==443)
		    {
		      ujGenIdxMu2 = GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]];
		      level = 2;
		      if(GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]>-1)
			{
			  if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]]])==BId)
			    {
			      level = 3;
			      bGenIdxMu2=GenInfo_mo1[GenInfo_mo1[MuonInfo_geninfo_index[BInfo_uj_rfmu2_index[BInfo_rfuj_index[j]]]]];
			      flagkstar++;///////////////////////////////////////////////////=2
			    }
			}
		    }
		}
	    }
	  Bgen[typesize]+=(level*1000);
	}

      int level=0;
      if(mGenIdxTk1!=-1 && mGenIdxTk2!=-1)
	{
	  if(!twoTks) level=1;
	  else
	    {
	      if(mGenIdxTk1==mGenIdxTk2) level=1;
	    }
	}
      if(bGenIdxMu1!=-1 && bGenIdxMu1==bGenIdxMu2 && bGenIdxMu1==bGenIdxTk1)
	{
	  if(!twoTks)
	    {
	      level=2;
	      BgenIndex[typesize] = bGenIdxMu1;
	    }
	  else if(bGenIdxMu1==bGenIdxTk2)
	    {
	      level=2;
	      BgenIndex[typesize] = bGenIdxMu1;
	    }
	}
      Bgen[typesize]+=(level*10000);

      //kstar#############################################################################
      if(kStar)
	{
	  //tk1
	  if(TrackInfo_geninfo_index[BInfo_rftk1_index[j]]>-1)
	    {
	      if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]])==tk2Id)
		{
		  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]>-1)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]])==MId)
			{
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]]])==BId)
				{
				  flagkstar++;//////////////////////////////////////////////=3
				  bGenIdxTk1=GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]]];
				}
			    }
			  mGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk1_index[j]]];
			}
		    }
		}
	    }
	  
	  //tk2
	  if(TrackInfo_geninfo_index[BInfo_rftk2_index[j]]>-1)
	    {
	      if(abs(GenInfo_pdgId[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]])==tk1Id)
		{
		  if(GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]>-1)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]])==MId)
			{
			  if(GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]>-1)
			    {
			      if(abs(GenInfo_pdgId[GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]]])==BId)
				{
				  flagkstar++;////////////////////////////////////////////////////=4
				  bGenIdxTk2 = GenInfo_mo1[GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]]];
				}
			    }
			  mGenIdxTk2 = GenInfo_mo1[TrackInfo_geninfo_index[BInfo_rftk2_index[j]]];
			}
		    }
		}
	    }
	  if(flagkstar==4)
	    {
	      if((bGenIdxMu1!=-1) 
		 && (bGenIdxMu1==bGenIdxMu2)
		 && (bGenIdxMu1==bGenIdxTk1)
		 && (bGenIdxMu1==bGenIdxTk2)
		 )
		{
		  Bgen[typesize]=41000;
		}
	    }
	}//kstar End#############################################################################

      int tgenIndex=BgenIndex[typesize];
      if(Bgen[typesize]==23333 || Bgen[typesize]==41000)
	{
	  BgenpdgId[typesize] = GenInfo_pdgId[tgenIndex];
	  Bgenpt[typesize] = GenInfo_pt[tgenIndex];
	  Bgeneta[typesize] = GenInfo_eta[tgenIndex];
	  Bgenphi[typesize] = GenInfo_phi[tgenIndex];
	  b4P->SetXYZM(GenInfo_pt[tgenIndex]*cos(GenInfo_phi[tgenIndex]),
		       GenInfo_pt[tgenIndex]*sin(GenInfo_phi[tgenIndex]),
		       GenInfo_pt[tgenIndex]*sinh(GenInfo_eta[tgenIndex]),
		       GenInfo_mass[tgenIndex]);
	  Bgeny[typesize] = b4P->Rapidity();
	}
    }
}











int signalGen(int Btype, int j)
{
  float BId,MId,tk1Id,tk2Id;
  int twoTks;
  //tk1:positive, tk2:negtive
  if(Btype==1)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 321;//K+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==2)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 211;//pi+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==3)
    {
      BId = 511;//B0
      MId = 310;//Ks
      tk1Id = 211;//pi+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==4)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 321;//K+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==5)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = -321;//pi+
      tk2Id = 211;//K-
      twoTks = 1;
    }
  if(Btype==6)
    {
      BId = 531;//Bs
      MId = 333;//phi
      tk1Id = 321;//K+
      tk2Id = -321;//K-
      twoTks = 1;
    }

  int flag=0;
  if (abs(GenInfo_pdgId[j])==BId&&GenInfo_nDa[j]==2&&GenInfo_da1[j]!=-1&&GenInfo_da2[j]!=-1)
    {
      if (abs(GenInfo_pdgId[GenInfo_da1[j]]==443))//jpsi
	{
	  if(GenInfo_da1[GenInfo_da1[j]]!=-1&&GenInfo_da2[GenInfo_da1[j]]!=-1)
	    {
	      if(abs(GenInfo_pdgId[GenInfo_da1[GenInfo_da1[j]]])==13&&abs(GenInfo_pdgId[GenInfo_da2[GenInfo_da1[j]]])==13)
		{
		  if(!twoTks)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_da2[j]])==tk1Id) flag++;
		    }
		  else
		    {
		      if (abs(GenInfo_pdgId[GenInfo_da2[j]])==MId) 
			{
			  if(GenInfo_da1[GenInfo_da2[j]]!=-1 && GenInfo_da2[GenInfo_da2[j]]!=-1)
			    {
			      if(GenInfo_pdgId[GenInfo_da1[GenInfo_da2[j]]]==tk1Id && GenInfo_pdgId[GenInfo_da2[GenInfo_da2[j]]]==tk2Id) flag++;
			    }
			}
		    }
		}
	    }
	}
    }
  return flag;
}

