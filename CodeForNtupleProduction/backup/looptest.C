#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TString.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "looptest.h"
#include <iomanip>

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

#define FILE_NUM 500

int looptest(TString infile="/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple/DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat0_TuneZ2_Unquenched_2760GeV_KpPim_20150902/Bfinder_PbPb_all_300_2_aU6.root", TString outfile="comp0.root", bool REAL=false, int startEntries=0,
	     TString kp="KmPip",TString pthat="0")
{
  double findMass(Int_t particlePdgId);
  void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, Int_t j, Int_t typesize, Bool_t REAL);
  bool isDsignalGen(Int_t Dtype, Int_t j);

  if(REAL) cout<<"--- REAL DATA ---"<<endl;
  else cout<<"--- MC ---"<<endl;
  
  //File type
  //TFile *f = new TFile(infile);
  ifstream readfilelist(Form("filelist_pthat%s_%s.txt",pthat.Data(),kp.Data()));
  int filename0[FILE_NUM],filename1[FILE_NUM];
  TString filename2[FILE_NUM];
  int ireadfile=0;
  for(ireadfile=0;ireadfile<FILE_NUM;ireadfile++)
    {
      readfilelist>>filename0[ireadfile];
      readfilelist>>filename1[ireadfile];
      readfilelist>>filename2[ireadfile];
      if(filename1[ireadfile]!=1 && filename1[ireadfile]!=2) break;
    }
  cout<<ireadfile<<endl;
  TChain* root = new TChain("Dfinder/root");
  for(int ifile=0;ifile<ireadfile;ifile++)
    {
      root->Add(Form("/mnt/hadoop/cms/store/user/twang/HI_DfinderNtuple/DinderMC_Pyquen_D0tokaonpion_D0pt1p0_Pthat%s_TuneZ2_Unquenched_2760GeV_%s_20150902/Bfinder_PbPb_all_%i_%i_%s.root",pthat.Data(),kp.Data(),filename0[ifile],filename1[ifile],filename2[ifile].Data()));
    }
  setDBranch(root);

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
  //TTree* ntGen = new TTree("ntGen","");         buildGenBranch(ntGen);
  cout<<"--- Building trees finished ---"<<endl;

  Long64_t nentries = root->GetEntries();
  //nentries = 10000;
  Long64_t nbytes = 0;
  TVector3* bP = new TVector3;
  TVector3* bVtx = new TVector3;
  TLorentzVector* b4P = new TLorentzVector;
  TLorentzVector* bGen = new TLorentzVector;

  Int_t itest=0;
  for (Int_t i=startEntries;i<nentries;i++)
    {
      if (i%10000==0) cout <<i<<" / "<<nentries<<endl;
      nbytes += root->GetEntry(i);
      Int_t Dtypesize[5]={0,0,0,0,0};
      Double_t Dbest=0;
      Int_t Dbestindex=0;
      for(int t=0;t<5;t++)
	{
	  Dsize=0;
	  if(isDchannel[t]==1)
	    {
	      Dbest=-1;
	      Dbestindex=-1;
	      for(int j=0;j<DInfo_size;j++)
		{
		  if(DInfo_pt[j]<=2.) continue;
		  if(TrackInfo_pt[DInfo_rftk1_index[j]]<=1.5&& (t==0||t==1)) continue;
		  if(TrackInfo_pt[DInfo_rftk2_index[j]]<=1.5&& (t==0||t==1)) continue;
		  
		  double d03D=TMath::Sqrt((DInfo_vtxX[j]-EvtInfo_PVx)*(DInfo_vtxX[j]-EvtInfo_PVx)+(DInfo_vtxY[j]-EvtInfo_PVy)*(DInfo_vtxY[j]-EvtInfo_PVy)+(DInfo_vtxZ[j]-EvtInfo_PVz)*(DInfo_vtxZ[j]-EvtInfo_PVz));
		  double d03Derr=TMath::Sqrt(DInfo_vtxXE[j]*DInfo_vtxXE[j]+DInfo_vtxYE[j]*DInfo_vtxYE[j]+DInfo_vtxZE[j]*DInfo_vtxZE[j]);
  		  if(d03D/d03Derr<2.5 && (t==0||t==1)) continue;
		  
		  if(DInfo_type[j]==(t+1))
		    {
		      fillDTree(bP,bVtx,b4P,j,Dtypesize[t],REAL);
		      if(Dchi2cl[Dtypesize[t]]>Dbest)
		      	{
		      	  Dbest = Dchi2cl[Dtypesize[t]];
		      	  Dbestindex = Dtypesize[t];
		      	}
		      Dtypesize[t]++;
		    }
		}
	      if(Dbestindex>-1)
	      	{
	      	  Disbestchi2[Dbestindex] = true;
	      	}
	      if(t==0) ntD1->Fill();
	      else if(t==1) ntD2->Fill();
	      else if(t==2) ntD3->Fill();
	      else if(t==3) ntD4->Fill();
	      else if(t==4) ntD5->Fill();
	    }
	}
    }
  cout<<itest<<endl;

  return 1;
}


double findMass(Int_t particlePdgId)
{
  if(TMath::Abs(particlePdgId)==211) return PION_MASS;
  if(TMath::Abs(particlePdgId)==321) return KAON_MASS;
  else
    {
      cout<<"ERROR: find particle mass falied >> Particle pdgId: "<<particlePdgId<<endl;
      return 0;
    }
}

void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, Int_t j, Int_t typesize, Bool_t REAL)
{

  Int_t hypo=-1,DpdgId=0,swap=-1;
  int level=0;
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
  if(!REAL)
    {
      if(DInfo_type[j]==1||DInfo_type[j]==2||DInfo_type[j]==3||DInfo_type[j]==4||DInfo_type[j]==5)
	{
	  if(TrackInfo_geninfo_index[DInfo_rftk1_index[j]]>-1)
	    {
	      level=0;
	      if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==DInfo_rftk1_MassHypo[j])
		{
		  hypo=-1;
		  swap=0;
		  if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==KAON_PDGID) hypo=0;
		  else if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==PION_PDGID) hypo=1;
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
	      else if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==DInfo_rftk2_MassHypo[j] && (DInfo_type[j]==1||DInfo_type[j]==2))
		{
		  hypo=-1;
		  swap=1;
                  if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==KAON_PDGID) hypo=0;
                  else if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])==PION_PDGID) hypo=1;
                  if(hypo==0||hypo==1)
                    {
                      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]>-1)
                        {
                          if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]])==DpdgId)
                            {
                              level=4;
                              dGenIdxTk1=GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]];
                            }
                        }
                    }
		}
	      Dgen[typesize]=level;
	    }
	  if(TrackInfo_geninfo_index[DInfo_rftk2_index[j]]>-1)
	    {
	      level=0;
	      if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==DInfo_rftk2_MassHypo[j] && swap==0)
		{
		  if((hypo==0&&TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==PION_PDGID)||(hypo==1&&TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==KAON_PDGID))
		    {
		      level=1;
		      if(hypo==1&&TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==KAON_PDGID) hypo=2;
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
	      else if(TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==DInfo_rftk1_MassHypo[j] && (DInfo_type[j]==1||DInfo_type[j]==2) && swap==1)
                {
		  if((hypo==0&&TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==PION_PDGID)||(hypo==1&&TMath::Abs(GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])==KAON_PDGID))
                    {
                      if(GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]>-1)
                        {
                          if(TMath::Abs(GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]])==DpdgId)
                            {
                              level=4;
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
		  level=0;
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
			  else if(hypo==3&&TMath::Abs(DInfo_rftk3_MassHypo[j])==PION_PDGID&&DInfo_type[j]==5)
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
		      level=0;
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
	  if((Dgen[typesize]==3333||Dgen[typesize]==3344)&&dGenIdxTk1==dGenIdxTk2)
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
      if(Dgen[typesize]==23333)
	{
	  if(DInfo_mass[j]<1.80)
	    {
	      TVector2 tk1v2,tk2v2,dv2;
	      tk1v2.SetMagPhi(GenInfo_pt[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]],GenInfo_phi[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]);
	      tk2v2.SetMagPhi(GenInfo_pt[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]],GenInfo_phi[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]);
	      dv2=tk1v2+tk2v2;
	      TLorentzVector tk1v4,tk2v4,dv4;
	      tk1v4.SetPtEtaPhiM(DInfo_rftk1_pt[j],DInfo_rftk1_eta[j],DInfo_rftk1_phi[j],GenInfo_mass[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]);
	      tk2v4.SetPtEtaPhiM(DInfo_rftk2_pt[j],DInfo_rftk2_eta[j],DInfo_rftk2_phi[j],GenInfo_mass[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]);
	      dv4.SetPxPyPzE(tk1v4.Px()+tk2v4.Px(),tk1v4.Py()+tk2v4.Py(),tk1v4.Pz()+tk2v4.Pz(),tk1v4.Energy()+tk2v4.Energy());
	      cout<<setiosflags(std::ios::left)<<"--"
		  <<setw(20)<<Form(" TkpdgId(%i,%i)",GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]],GenInfo_pdgId[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])
		  <<setw(14)<<Form(" Tkmo(%i,%i)",GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]],GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])
		  <<setw(22)<<Form(" TkmopdgId(%i,%i)",GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]],GenInfo_pdgId[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]])
		  <<setw(16)<<Form(" NoOfdau(%i,%i)",GenInfo_nDa[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]],GenInfo_nDa[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]]])
		  <<setw(24)<<Form(" DpTvscom(%.2f,%.2f)",GenInfo_pt[GenInfo_mo1[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]]],dv2.Mod())
		  <<setw(22)<<Form(" Tk1pTRvsG(%.2f,%.2f)",TrackInfo_pt[DInfo_rftk1_index[j]],GenInfo_pt[TrackInfo_geninfo_index[DInfo_rftk1_index[j]]])
		  <<setw(22)<<Form(" Tk2pTRvsG(%.2f,%.2f)",TrackInfo_pt[DInfo_rftk2_index[j]],GenInfo_pt[TrackInfo_geninfo_index[DInfo_rftk2_index[j]]])
		  <<setw(26)<<Form(" Tk1pTBfvsAf(%.2f,%.2f)",TrackInfo_pt[DInfo_rftk1_index[j]],DInfo_rftk1_pt[j])
		  <<setw(26)<<Form(" Tk2pTBfvsAf(%.2f,%.2f)",TrackInfo_pt[DInfo_rftk2_index[j]],DInfo_rftk2_pt[j])
		  <<setw(24)<<Form(" Dmsvscom(%.2f,%.2f)",DInfo_mass[j],dv4.M())
		  <<endl;
	    }
	}
    }//if(!real)
}//fillDtree

//
bool isDsignalGen(Int_t Dtype, Int_t j)
{
  bool flag=false;
  if(Dtype==1||Dtype==2)
    {
      if(abs(GenInfo_pdgId[j])==DZERO_PDGID&&GenInfo_nDa[j]==2&&GenInfo_da1[j]!=-1&&GenInfo_da2[j]!=-1)
	{
	  if((TMath::Abs(GenInfo_pdgId[GenInfo_da1[j]])==KAON_PDGID&&TMath::Abs(GenInfo_pdgId[GenInfo_da2[j]])==PION_PDGID) || 
	     (TMath::Abs(GenInfo_pdgId[GenInfo_da1[j]])==PION_PDGID&&TMath::Abs(GenInfo_pdgId[GenInfo_da2[j]])==KAON_PDGID))
	    {
	      flag=true;
	    }
	}
    }
  return flag;
}

