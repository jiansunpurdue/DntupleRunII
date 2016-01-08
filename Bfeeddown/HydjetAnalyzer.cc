// -*- C++ -*-
//
// Package:    HydjetAnalyzer
// Class:      HydjetAnalyzer
// 
/**\class HydjetAnalyzer HydjetAnalyzer.cc yetkin/HydjetAnalyzer/src/HydjetAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Yetkin Yilmaz
//         Created:  Tue Dec 18 09:44:41 EST 2007
// $Id: HydjetAnalyzer.cc,v 1.24 2011/01/21 15:45:18 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// root include file
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

static const int MAXPARTICLES = 5000000;
static const int MAXVTX = 1000;
static const int ETABINS = 3; // Fix also in branch string

//
// class decleration
//

struct HydjetEvent{

   int event;
   float b;
   float npart;
   float ncoll;
   float nhard;
   float phi0;
   float scale;

   int n[ETABINS];
   float ptav[ETABINS];

   int mult;
   float pt[MAXPARTICLES];
   float eta[MAXPARTICLES];
   float phi[MAXPARTICLES];
   int pdg[MAXPARTICLES];
   int chg[MAXPARTICLES];

   int momPdg[MAXPARTICLES];
   float momPt[MAXPARTICLES];
   float momEta[MAXPARTICLES];
   float momPhi[MAXPARTICLES];

   int gmomPdg[MAXPARTICLES];
   float gmomPt[MAXPARTICLES];
   float gmomEta[MAXPARTICLES];
   float gmomPhi[MAXPARTICLES];

   int ggmomPdg[MAXPARTICLES];
   float ggmomPt[MAXPARTICLES];
   float ggmomEta[MAXPARTICLES];
   float ggmomPhi[MAXPARTICLES];

   int BPdg[MAXPARTICLES];
   float BPt[MAXPARTICLES];
   float BEta[MAXPARTICLES];
   float BPhi[MAXPARTICLES];

   float vx;
   float vy;
   float vz;
   float vr;

};

class HydjetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HydjetAnalyzer(const edm::ParameterSet&);
      ~HydjetAnalyzer();


   private:
      virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   std::ofstream out_b;
   std::string fBFileName;

   std::ofstream out_n;
   std::string fNFileName;

   std::ofstream out_m;
   std::string fMFileName;

  
   TTree* hydjetTree_;
   HydjetEvent hev_;

   TNtuple *nt;

   std::string output;           // Output filename
 
   bool doAnalysis_;
   bool printLists_;
   bool doCF_;
   bool doVertex_;
  bool useHepMCProduct_;
  bool doHI_;
  bool doParticles_;

   double etaMax_;
   double ptMin_;
  edm::InputTag src_;
  edm::InputTag genParticleSrc_;
  edm::InputTag genHIsrc_;
   edm::InputTag simVerticesTag_;

   edm::ESHandle < ParticleDataTable > pdt;
   edm::Service<TFileService> f;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HydjetAnalyzer::HydjetAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   fBFileName = iConfig.getUntrackedParameter<std::string>("output_b", "b_values.txt");
   fNFileName = iConfig.getUntrackedParameter<std::string>("output_n", "n_values.txt");
   fMFileName = iConfig.getUntrackedParameter<std::string>("output_m", "m_values.txt");
   doAnalysis_ = iConfig.getUntrackedParameter<bool>("doAnalysis", true);
   useHepMCProduct_ = iConfig.getUntrackedParameter<bool>("useHepMCProduct", false);
   printLists_ = iConfig.getUntrackedParameter<bool>("printLists", false);
   doCF_ = iConfig.getUntrackedParameter<bool>("doMixed", false);
   doVertex_ = iConfig.getUntrackedParameter<bool>("doVertex", false);
   if (doVertex_) {
      simVerticesTag_ = iConfig.getParameter<edm::InputTag>("simVerticesTag");
   }
   etaMax_ = iConfig.getUntrackedParameter<double>("etaMax", 2);
   ptMin_ = iConfig.getUntrackedParameter<double>("ptMin", 0);
   src_ = iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("generator"));
   genParticleSrc_ =iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("genParticles"));
   genHIsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("heavyIon"));
   doParticles_ = iConfig.getUntrackedParameter<bool>("doParticles", true);

}


HydjetAnalyzer::~HydjetAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HydjetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace HepMC;
  
   hev_.event = iEvent.id().event();
   for(int ieta = 0; ieta < ETABINS; ++ieta){
      hev_.n[ieta] = 0;
      hev_.ptav[ieta] = 0;
   }
   hev_.mult = 0;
      
   double phi0 = 0;
   double b = -1;
   double scale = -1;
   int npart = -1;
   int ncoll = -1;
   int nhard = -1;
   double vx = -99;
   double vy = -99;
   double vz = -99;
   double vr = -99;
   const GenEvent* evt;
  
   int nmix = -1;
   int np = 0;
   int sig = -1;
   int src = -1;

   if(useHepMCProduct_){
   if(doCF_){
     Handle<CrossingFrame<HepMCProduct> > cf;
     iEvent.getByLabel(InputTag("mix","source"),cf);
     MixCollection<HepMCProduct> mix(cf.product());
     nmix = mix.size();
     cout<<"Mix Collection Size: "<<mix<<endl;

     MixCollection<HepMCProduct>::iterator mbegin = mix.begin();
     MixCollection<HepMCProduct>::iterator mend = mix.end();
     
     for(MixCollection<HepMCProduct>::iterator mixit = mbegin; mixit != mend; ++mixit){
       const GenEvent* subevt = (*mixit).GetEvent();
       int all = subevt->particles_size();
       np += all;
       HepMC::GenEvent::particle_const_iterator begin = subevt->particles_begin();
       HepMC::GenEvent::particle_const_iterator end = subevt->particles_end();
       for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
//	 if((*it)->status() == 1){
	   int pdg_id = (*it)->pdg_id();
	   float eta = (*it)->momentum().eta();
	   float phi = (*it)->momentum().phi();
	   float pt = (*it)->momentum().perp();
	   const ParticleData * part = pdt->particle(pdg_id );
	   int charge = static_cast<int>(part->charge());

	   hev_.pt[hev_.mult] = pt;
	   hev_.eta[hev_.mult] = eta;
	   hev_.phi[hev_.mult] = phi;
	   hev_.pdg[hev_.mult] = pdg_id;
	   hev_.chg[hev_.mult] = charge;

	   eta = fabs(eta);
	   int etabin = 0;
	   if(eta > 0.5) etabin = 1;
	   if(eta > 1.) etabin = 2;
	   if(eta < 2.){
	     hev_.ptav[etabin] += pt;
	     ++(hev_.n[etabin]);
	   }
	   ++(hev_.mult);
	   }
//       }
     }
   }else{
   
   Handle<HepMCProduct> mc;
   iEvent.getByLabel(src_,mc);
   evt = mc->GetEvent();
   scale = evt->event_scale();
      
   const HeavyIon* hi = evt->heavy_ion();
   if(hi){
      b = hi->impact_parameter();
      npart = hi->Npart_proj()+hi->Npart_targ();
      ncoll = hi->Ncoll();
      nhard = hi->Ncoll_hard();
      phi0 = hi->event_plane_angle();
      
      if(printLists_){
	 out_b<<b<<endl;
	 out_n<<npart<<endl;
      }
   }
   
   src = evt->particles_size();
   
   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
     if((*it)->status() == 1){
       int pdg_id = (*it)->pdg_id();
       float eta = (*it)->momentum().eta();
       float phi = (*it)->momentum().phi();
       float pt = (*it)->momentum().perp();
       const ParticleData * part = pdt->particle(pdg_id );
       int charge = static_cast<int>(part->charge());
       
       hev_.pt[hev_.mult] = pt;
       hev_.eta[hev_.mult] = eta;
       hev_.phi[hev_.mult] = phi;
       hev_.pdg[hev_.mult] = pdg_id;
       hev_.chg[hev_.mult] = charge;
       
       eta = fabs(eta);
       int etabin = 0;
       if(eta > 0.5) etabin = 1; 
       if(eta > 1.) etabin = 2;
       if(eta < 2.){
	 hev_.ptav[etabin] += pt;
	 ++(hev_.n[etabin]);
	 }
       ++(hev_.mult);
     }
   }
   }
   }else{
     edm::Handle<reco::GenParticleCollection> parts;
     iEvent.getByLabel(genParticleSrc_,parts);
     for(unsigned int i = 0; i < parts->size(); ++i){
       const reco::GenParticle& p = (*parts)[i];
       hev_.pt[hev_.mult] = p.pt();
       hev_.eta[hev_.mult] = p.eta();
       hev_.phi[hev_.mult] = p.phi();
       hev_.pdg[hev_.mult] = p.pdgId();
       hev_.chg[hev_.mult] = p.charge();
       hev_.BPdg[hev_.mult] = 0;
       hev_.BPt[hev_.mult] = -1;
       hev_.BEta[hev_.mult] = 0;
       hev_.BPhi[hev_.mult] = 0;
       hev_.momPdg[hev_.mult] = 0;
       hev_.momPt[hev_.mult] = -1;
       hev_.momEta[hev_.mult] = 0;
       hev_.momPhi[hev_.mult] = 0;
       hev_.gmomPdg[hev_.mult] = 0;
       hev_.gmomPt[hev_.mult] = -1;
       hev_.gmomEta[hev_.mult] = 0;
       hev_.gmomPhi[hev_.mult] = 0;
       hev_.ggmomPdg[hev_.mult] = 0;
       hev_.ggmomPt[hev_.mult] = -1;
       hev_.ggmomEta[hev_.mult] = 0;
       hev_.ggmomPhi[hev_.mult] = 0;
       if (fabs(p.pdgId())!=421&&fabs(p.pdgId())!=443) continue;

       if (p.mother()!=0){
          hev_.momPdg[hev_.mult] = p.mother()->pdgId();
          hev_.momPt[hev_.mult] = p.mother()->pt();
          hev_.momEta[hev_.mult] = p.mother()->eta();
          hev_.momPhi[hev_.mult] = p.mother()->phi();
          if ((fabs(hev_.momPdg[hev_.mult])>=500&&fabs(hev_.momPdg[hev_.mult])<600)
            ||(fabs(hev_.momPdg[hev_.mult])>=5000&&fabs(hev_.momPdg[hev_.mult])<6000)) {
             hev_.BPdg[hev_.mult] = p.mother()->pdgId();
             hev_.BPt[hev_.mult] = p.mother()->pt();
             hev_.BEta[hev_.mult] = p.mother()->eta();
             hev_.BPhi[hev_.mult] = p.mother()->phi();
          }
        
          if (p.mother()->mother()!=0){
             hev_.gmomPdg[hev_.mult] = p.mother()->mother()->pdgId();
             hev_.gmomPt[hev_.mult] = p.mother()->mother()->pt();
             hev_.gmomEta[hev_.mult] = p.mother()->mother()->eta();
             hev_.gmomPhi[hev_.mult] = p.mother()->mother()->phi();
             if ((fabs(hev_.gmomPdg[hev_.mult])>=500&&fabs(hev_.gmomPdg[hev_.mult])<600)
               ||(fabs(hev_.gmomPdg[hev_.mult])>=5000&&fabs(hev_.gmomPdg[hev_.mult])<6000)) {
                hev_.BPdg[hev_.mult] = p.mother()->mother()->pdgId();
                hev_.BPt[hev_.mult] = p.mother()->mother()->pt();
                hev_.BEta[hev_.mult] = p.mother()->mother()->eta();
                hev_.BPhi[hev_.mult] = p.mother()->mother()->phi();
             }

             if (p.mother()->mother()->mother()!=0){
                hev_.ggmomPdg[hev_.mult] = p.mother()->mother()->mother()->pdgId();
                hev_.ggmomPt[hev_.mult] = p.mother()->mother()->mother()->pt();
                hev_.ggmomEta[hev_.mult] = p.mother()->mother()->mother()->eta();
                hev_.ggmomPhi[hev_.mult] = p.mother()->mother()->mother()->phi();
                if ((fabs(hev_.ggmomPdg[hev_.mult])>=500&&fabs(hev_.ggmomPdg[hev_.mult])<600)
                  ||(fabs(hev_.ggmomPdg[hev_.mult])>=5000&&fabs(hev_.ggmomPdg[hev_.mult])<6000)) {
                   hev_.BPdg[hev_.mult] = p.mother()->mother()->mother()->pdgId();
                   hev_.BPt[hev_.mult] = p.mother()->mother()->mother()->pt();
                   hev_.BEta[hev_.mult] = p.mother()->mother()->mother()->eta();
                   hev_.BPhi[hev_.mult] = p.mother()->mother()->mother()->phi();
                }
   
             }
          }
       }
       double eta = fabs(p.eta());

       int etabin = 0;
       if(eta > 0.5) etabin = 1;
       if(eta > 1.) etabin = 2;
       if(eta < 2.){
         hev_.ptav[etabin] += p.pt();
         ++(hev_.n[etabin]);
       }
       ++(hev_.mult);
     }
     if(doHI_){
       edm::Handle<GenHIEvent> higen;
       iEvent.getByLabel(genHIsrc_,higen);
     }
   }

   if(doVertex_){
      edm::Handle<edm::SimVertexContainer> simVertices;
      iEvent.getByLabel<edm::SimVertexContainer>(simVerticesTag_, simVertices);
      
      if (! simVertices.isValid() ) throw cms::Exception("FatalError") << "No vertices found\n";
      int inum = 0;
      
      edm::SimVertexContainer::const_iterator it=simVertices->begin();
      SimVertex vertex = (*it);
      cout<<" Vertex position "<< inum <<" " << vertex.position().rho()<<" "<<vertex.position().z()<<endl;
      vx = vertex.position().x();
      vy = vertex.position().y();
      vz = vertex.position().z();
      vr = vertex.position().rho();
   }
   
   for(int i = 0; i<3; ++i){
      hev_.ptav[i] = hev_.ptav[i]/hev_.n[i];
   }
   
   hev_.b = b;
   hev_.scale = scale;
   hev_.npart = npart;
   hev_.ncoll = ncoll;
   hev_.nhard = nhard;
   hev_.phi0 = phi0;
   hev_.vx = vx;
   hev_.vy = vy;
   hev_.vz = vz;
   hev_.vr = vr;

   nt->Fill(nmix,np,src,sig);

   hydjetTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
HydjetAnalyzer::beginRun(const edm::Run&, const edm::EventSetup& iSetup) 
{
   iSetup.getData(pdt);
}

void 
HydjetAnalyzer::beginJob()
{

   if(printLists_){
      out_b.open(fBFileName.c_str());
      if(out_b.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fBFileName;
      out_n.open(fNFileName.c_str());
      if(out_n.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fNFileName;
      out_m.open(fMFileName.c_str());
      if(out_m.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fMFileName;
   }   
   
   if(doAnalysis_){
      nt = f->make<TNtuple>("nt","Mixing Analysis","mix:np:src:sig");

      hydjetTree_ = f->make<TTree>("hi","Tree of Hydjet Events");
      hydjetTree_->Branch("event",&hev_.event,"event/I");
      hydjetTree_->Branch("b",&hev_.b,"b/F");
      hydjetTree_->Branch("npart",&hev_.npart,"npart/F");
      hydjetTree_->Branch("ncoll",&hev_.ncoll,"ncoll/F");
      hydjetTree_->Branch("nhard",&hev_.nhard,"nhard/F");
      hydjetTree_->Branch("phi0",&hev_.phi0,"phi0/F");
      hydjetTree_->Branch("scale",&hev_.scale,"scale/F");

      hydjetTree_->Branch("n",hev_.n,"n[3]/I");
      hydjetTree_->Branch("ptav",hev_.ptav,"ptav[3]/F");

      if(doParticles_){
	
	hydjetTree_->Branch("mult",&hev_.mult,"mult/I");
	hydjetTree_->Branch("pt",hev_.pt,"pt[mult]/F");
	hydjetTree_->Branch("eta",hev_.eta,"eta[mult]/F");
	hydjetTree_->Branch("phi",hev_.phi,"phi[mult]/F");
	hydjetTree_->Branch("pdg",hev_.pdg,"pdg[mult]/I");
	hydjetTree_->Branch("chg",hev_.chg,"chg[mult]/I");

	hydjetTree_->Branch("momEta",hev_.momEta,"momEta[mult]/F");
	hydjetTree_->Branch("momPhi",hev_.momPhi,"momPhi[mult]/F");
	hydjetTree_->Branch("momPdg",hev_.momPdg,"momPdg[mult]/I");
	hydjetTree_->Branch("momPt",hev_.momPt,"momPt[mult]/F");

	hydjetTree_->Branch("gmomEta",hev_.gmomEta,"gmomEta[mult]/F");
	hydjetTree_->Branch("gmomPhi",hev_.gmomPhi,"gmomPhi[mult]/F");
	hydjetTree_->Branch("gmomPdg",hev_.gmomPdg,"gmomPdg[mult]/I");
	hydjetTree_->Branch("gmomPt",hev_.gmomPt,"gmomPt[mult]/F");

	hydjetTree_->Branch("ggmomEta",hev_.ggmomEta,"ggmomEta[mult]/F");
	hydjetTree_->Branch("ggmomPhi",hev_.ggmomPhi,"ggmomPhi[mult]/F");
	hydjetTree_->Branch("ggmomPdg",hev_.ggmomPdg,"ggmomPdg[mult]/I");
	hydjetTree_->Branch("ggmomPt",hev_.ggmomPt,"ggmomPt[mult]/F");

	hydjetTree_->Branch("BEta",hev_.BEta,"BEta[mult]/F");
	hydjetTree_->Branch("BPhi",hev_.BPhi,"BPhi[mult]/F");
	hydjetTree_->Branch("BPdg",hev_.BPdg,"BPdg[mult]/I");
	hydjetTree_->Branch("BPt",hev_.BPt,"BPt[mult]/F");
	
	hydjetTree_->Branch("vx",&hev_.vx,"vx/F");
	hydjetTree_->Branch("vy",&hev_.vy,"vy/F");
	hydjetTree_->Branch("vz",&hev_.vz,"vz/F");
	hydjetTree_->Branch("vr",&hev_.vr,"vr/F");
      }
   }
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HydjetAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HydjetAnalyzer);
