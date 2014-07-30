#include <iostream>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2D.h"
#include <TMath.h>
#include <Math/VectorUtil.h>
#include <TRandom.h>
#include "TLorentzVector.h"

#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

class GENAnalyzer : public edm::EDAnalyzer
{

   public:
   
      //
      explicit GENAnalyzer( const edm::ParameterSet& ) ;
      virtual ~GENAnalyzer() {} // no need to delete ROOT stuff
                                   // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endRun( const edm::Run&, const edm::EventSetup& ) ;
      virtual void endJob() ;

   private:
   
     TH1D*       fHistNeutMass ;
   //  TH1D*       fHist4muMass ;
   //  TH1D*       fHistZZMass ;
     
     TH1D*       timeNeut ;
     TH1D*       betaNeut ;
     TH1D*       ptNeut ;
     TH1D*       massNeut ;
     TH1D*       eNeut ;
     
     TH1D*       pho_time ;
     TH1D*       gamma_time ;
     TH1D*       pho_T ;
     TH1D*       pho_pt ;
     TH1D*       pho_mass ;
     TH1D*       pho_N ;
     TH1D*       gravpt_hist ;
     TH1D*       N_neut,*N_grav;
     TH2D*       probVspt,*probVsdr,*probVsptVsdr;
     TH1D*       NeutCtau ;
}; 

using namespace edm;
using namespace std;

GENAnalyzer::GENAnalyzer( const ParameterSet& pset )
  : fHistNeutMass(0),timeNeut(0),betaNeut(0),ptNeut(0),massNeut(0),eNeut(0),pho_time(0),gamma_time(0),pho_T(0),pho_pt(0),pho_mass(0),pho_N(0),gravpt_hist(0),N_neut(0),N_grav(0),probVspt(0),probVsdr(0),probVsptVsdr(0), NeutCtau(0)
{
// actually, pset is NOT in use - we keep it here just for illustratory putposes
}

void GENAnalyzer::beginJob()
{
  
  Service<TFileService> fs;
  fHistNeutMass = fs->make<TH1D>(  "HistNeutMass", "Photon-Gravitino inv. mass", 100,  0., 1500. ) ;
  //fHist4muMass = fs->make<TH1D>(  "Hist4muMass", "4-mu inv. mass", 100, 170., 210. ) ;
  //fHistZZMass  = fs->make<TH1D>(  "HistZZMass",  "ZZ inv. mass",   100, 170., 210. ) ;    
    
  NeutCtau = fs->make<TH1D>(  "NeutCtau", "c#tau_{#chi^{0}_{1}} = |dr|/#gamma#beta", 100, -1.0, 20000. ) ;
  timeNeut = fs->make<TH1D>(  "TimeNeut", "Neutralino c#tau(mm)", 100,  -1.0, 20000. ) ;
  betaNeut = fs->make<TH1D>(  "BetaNeut", "Neutralino #beta = |p|/E", 100, 0., 1.0 ) ;
  ptNeut = fs->make<TH1D>(  "PtNeut", "Neutralino Pt(GeV/c)", 100,  0., 1500. ) ;
  massNeut = fs->make<TH1D>(  "MassNeut", "Neutralino Mass (GeV/c^2)", 100, 0., 500. ) ;
  eNeut = fs->make<TH1D>(  "EnergyNeut", "Neutralino Energy (GeV)", 100, 0., 1000. ) ;
  pho_time = fs->make<TH1D>(  "Photon_Time", "T_{#gamma}[ns]:|dr|/#beta", 100, -1., 100. ) ;
  gamma_time = fs->make<TH1D>(  "Gamma_Time", "T_{#gamma}[ns]:t_pho-t_Neut", 100, -1., 200. ) ;
  pho_T = fs->make<TH1D>(  "Direct_pho_time", "Photon MC Time[ns]:position.t()", 100, -1., 200. ) ;
  pho_pt = fs->make<TH1D>(  "Photon_Pt", "Photon Pt (GeV/c", 100, 0., 1500. ) ;
  pho_mass = fs->make<TH1D>(  "Photon_Mass", "Photon Mass (GeV/c^2)", 100, 0., 1000. ) ;
  pho_N = fs->make<TH1D>(  "Photon_N", "Number Of Photons", 100, 0, 10 ) ;
  gravpt_hist = fs->make<TH1D>(  "Gravitino_Pt", "Gravitino Pt (GeV/c", 100, 0., 1500. ) ;
  N_neut = fs->make<TH1D>(  "Neutralino_N", "Number Of Neutralinos", 100, 0, 10 ) ;
  N_grav = fs->make<TH1D>(  "Gravitino_N", "Number Of Gravitinos", 100, 0, 10 ) ;
  probVspt = fs->make<TH2D>(  "probVspt", "Decay Probability Vs Pt", 5000, 0.0, 500.0, 10, -0.05, 1.1 ) ;
  probVsdr = fs->make<TH2D>(  "probVsdr", "Decay Probability Vs drT[mm]", 20000, 0.0, 10000.0, 10, -0.05, 1.1 ) ;
  probVsptVsdr = fs->make<TH2D>(  "probVsptVsdr", "Decay Probability Vs Pt, Varying drT", 5000, 0.0, 500.0, 10, -0.05, 1.1 ) ;
  return ;
  
}

void GENAnalyzer::analyze( const Event& e, const EventSetup& )
{
  
  // here's an example of accessing GenEventInfoProduct
  Handle< GenEventInfoProduct > GenInfoHandle;
  e.getByLabel( "generator", GenInfoHandle );
  double qScale = GenInfoHandle->qScale();
  double pthat = ( GenInfoHandle->hasBinningValues() ? 
                  (GenInfoHandle->binningValues())[0] : 0.0);
  cout << " qScale = " << qScale << " pthat = " << pthat << endl;
  //
  // this (commented out) code below just exemplifies how to access certain info 
  //
  //double evt_weight1 = GenInfoHandle->weights()[0]; // this is "stanrd Py6 evt weight;
                                                    // corresponds to PYINT1/VINT(97)
  //double evt_weight2 = GenInfoHandle->weights()[1]; // in case you run in CSA mode or otherwise
                                                    // use PYEVWT routine, this will be weight
						    // as returned by PYEVWT, i.e. PYINT1/VINT(99)
  //std::cout << " evt_weight1 = " << evt_weight1 << std::endl;
  //std::cout << " evt_weight2 = " << evt_weight2 << std::endl;
  //double weight = GenInfoHandle->weight();
  //std::cout << " as returned by the weight() method, integrated event weight = " << weight << std::endl;
  
  // here's an example of accessing particles in the event record (HepMCProduct)
  //
  Handle< HepMCProduct > EvtHandle ;
  
  // find initial (unsmeared, unfiltered,...) HepMCProduct
  //
  e.getByLabel( "generator", EvtHandle ) ;
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  
   HepMC::GenVertex*   NeutDecVtx = 0 ;
  // HepMC::GenVertex*   NeutProdVtx = 0 ;
   HepMC::FourVector   NeutDecPos  ;
   HepMC::FourVector   NeutProdPos  ;
   HepMC::FourVector   NeutDecP4  ;
  // HepMC::GenVertex*   phoVtx = 0 ;
   HepMC::FourVector   phoP4 ;
   HepMC::FourVector   phoPos ;
   HepMC::FourVector   GravP4 ;
   HepMC::FourVector   GravPos ;
   double t_pho = -99.;
   double t_Neut = -99.;
   double ctau_neut = -99.;
   double beta = -99.;
   double t_gamma = -99.;
   double Gamma  = -99.;
   int  n_photons = -99999; 
   int  n_neut = -99999; 
   double ptN=0, mass=0 ;  // Nuetralino Pt & Mass
   vector<HepMC::GenParticle*> Neutralinos;
   vector<HepMC::GenParticle*> photons;
   vector<HepMC::GenParticle*> gravitinos;
  // find the 1st vertex with outgoing Higgs 
  // and get Higgs decay vertex from there;
  //
  // in principal, one can look for the vertex 
  // with incoming Higgs as well...
  //
for ( HepMC::GenEvent::vertex_const_iterator  vit=Evt->vertices_begin(); vit!=Evt->vertices_end(); vit++ )
  {
   for ( HepMC::GenVertex::particles_out_const_iterator pout=(*vit)->particles_out_const_begin();
            pout!=(*vit)->particles_out_const_end(); pout++ )
      {
        double Px, Py, Pz, e, m, Pt;
        if ( (*pout)->pdg_id() == 1000022 && (*pout)->status() == 3) 
          { 
           if ( (*pout)->end_vertex() != 0 )
              {
                NeutDecVtx = (*pout)->end_vertex() ;
                NeutProdPos = (*pout)->production_vertex()->position() ;
                NeutDecPos = (*pout)->end_vertex()->position() ;
                NeutDecP4 = (*pout)->momentum() ;

                Px = NeutDecP4.px() ;
                Py = NeutDecP4.py() ;
                Pz = NeutDecP4.pz() ;
                e  = NeutDecP4.e() ;
                //double  Dx = ((*pout)->end_vertex()->position().x() - (*pout)->production_vertex()->position().x() ) ;
                //double  Dy = ((*pout)->end_vertex()->position().y() - (*pout)->production_vertex()->position().y() ) ;
                //double  Dz = ((*pout)->end_vertex()->position().z() - (*pout)->production_vertex()->position().z() ) ;
                //double  drT  = sqrt(((Dx*Dx) + (Dy*Dy)));      
                //m  = NeutDecP4.m(); // (*it)->generatedMass() ;
                m  = (*pout)->generatedMass() ;
		mass = NeutDecP4.m();
                Pt = sqrt(((Px*Px) + (Py*Py) ));
		ptN = Pt;
                //double P = sqrt(((Px*Px) + (Py*Py) +(Pz*Pz)));
                beta  = (e == 0 )? 999 : sqrt(((Px*Px) + (Py*Py) + (Pz*Pz)))/e ;
                Gamma = e/m;
		double prob  = 1 - TMath::Exp(-((m*1.3)/(6*Pt)) );
		Neutralinos.push_back(*pout);
		n_neut++;
                // Fill Neutralino Info
                betaNeut->Fill(beta);
                ptNeut->Fill(Pt);
                massNeut->Fill(m); 
                eNeut->Fill(e); 
		probVspt->Fill(Pt, prob);

              }
  
     HepMC::FourVector Mom2part ;
     double            XMass2part = 0.;
  // double            XMass4part = 0.;
  // double            XMass2pairs  = 0.;
     vector< HepMC::FourVector > Mom2partCont ;
      HepMC::GenVertex* neut_p = NeutDecVtx;
for ( HepMC::GenVertex::particles_out_const_iterator des = neut_p->particles_out_const_begin();
	 des != neut_p->particles_out_const_end(); des++ )
  {
    // kinematics of Neutralino/ Photon
    // Now Loop through Stable Neutralino Descendants!
    // skip other than Photon && Gravitino
    // for Photon
    HepMC::GenVertex* des_p = (*des)->end_vertex();
    if ( (abs((*des)->pdg_id()) != 22 ) && (abs((*des)->pdg_id()) != 1000039) ) continue ; 
    for( HepMC::GenVertex::particles_out_const_iterator p_des = des_p->particles_out_const_begin();
                                                       p_des != des_p->particles_out_const_end(); p_des++)
     {
        if (abs((*des)->pdg_id()) == 22  && ((*des)->end_vertex()!=0 ))
            {
                // phoVtx = StableNeutDesc[i]->end_vertex();
                  //if ( e.id().event() == 1 ) (*des)->print() ;
                 if ( (*p_des)->status() != 1 ) continue ; 
                 phoP4  = HepMC::FourVector( (*p_des)->momentum().px(),(*p_des)->momentum().py(), (*p_des)->momentum().pz(), (*p_des)->momentum().e() ) ;
                 phoPos = HepMC::FourVector( des_p->position().x(), des_p->position().y(), des_p->position().z(), des_p->position().t() );                        
                 double   px = (*p_des)->momentum().px(); 
                 double   py = (*p_des)->momentum().py(); 
                 double   mp  = (*p_des)->momentum().m(); //generatedMass() ;
                // double   pt = Sqrt( ((*des)->momentum().px()*(*des)->momentum().px()) + ((*des)->momentum().py()*(*des)->momentum().py())); 
                 double   pt = sqrt(((px*px) + (py*py)));        
                 //double   t_p = (*des)->end_vertex()->position().t()/300 ; // 300(mm/ns) //- NeutDecPos.t(); //phoPos.t();     
                 double   t_p = des_p->position().t()/300 ; // 300(mm/ns) //- NeutDecPos.t(); //phoPos.t();     
                 double  dx = (des_p->position().x() - neut_p->position().x() ) ;
		 double  dy = (des_p->position().y() - neut_p->position().y() ) ;
		 double  dz = (des_p->position().z() - neut_p->position().z() );

                 //double  dr  = sqrt(((dx*dx) + (dy*dy) + (dz*dz)));     
                 double  drT  = sqrt(((dx*dx) + (dy*dy) ));      
		 double xprob = 1 - TMath::Exp(-((mass*drT)/(6000*ptN)) );
                 t_pho  = sqrt(((dx*dx) + (dy*dy) + (dz*dz)))/(beta*300);     
                 ctau_neut  = sqrt(((dx*dx) + (dy*dy) + (dz*dz)))/(beta*Gamma); // mm     
                 t_Neut = (des_p->position().t() -  neut_p->position().t())/Gamma;  //ctau[mm]
                 //t_Neut = (des_p->position().t())/Gamma;  //ctau[mm]
                 //t_Neut = neut_p->position().t()/Gamma;  //ctau[mm]
                 t_gamma = ( des_p->position().t() - neut_p->position().t() )/(300);  //t[ns]
                 //Fill Photon info
                 pho_T->Fill(t_p); // Raw Photon time
                 pho_pt->Fill(pt);
                 pho_mass->Fill(mp);
		 NeutCtau->Fill(ctau_neut);
                 timeNeut->Fill(t_Neut);
		 probVsdr->Fill(drT, xprob);
		 probVsptVsdr->Fill(ptN, xprob);
                 pho_time->Fill(t_pho); // dr/beta*C photon time;
                 gamma_time->Fill(t_gamma);
		 photons.push_back(*des);
                 n_photons++; 
		 }
                  
            
           if (abs((*des)->pdg_id()) == 1000039 && ((*des)->end_vertex()!=0 ))
              {
                //if ( (*des)->status() != 1 ) continue ;
                 GravP4 = HepMC::FourVector( (*des)->momentum().px(), (*des)->momentum().py(), (*des)->momentum().pz(), (*des)->momentum().e() );
              //  GravP4 = (*des)->momentum();
                GravPos = HepMC::FourVector( des_p->position().x(), des_p->position().y(), des_p->position().z(), des_p->position().t() ); 
		gravpt_hist->Fill( sqrt( (( (*des)->momentum().px()*(*des)->momentum().px()) + ( (*des)->momentum().py()*(*des)->momentum().py() ) ) ));
		gravitinos.push_back(*des);
              }
                      
	      Mom2part = HepMC::FourVector(( GravP4.px() + phoP4.px()),
	                                ( GravP4.py() + phoP4.py()),
					( GravP4.pz() + phoP4.pz()),
					( GravP4.e() + phoP4.e())) ;
             XMass2part = Mom2part.m() ;
             fHistNeutMass->Fill( XMass2part ) ;
           }
     pho_N->Fill(photons.size());
     N_neut->Fill(Neutralinos.size());
     N_grav->Fill(gravitinos.size());
        }
                  
        break ;
       }
    } // eof of Loop over all Out Particles from vertext
 } 
/*
// Check info of first event
if ( e.id().event() == 1 )
  {
      cout << " " << endl ;
      cout << " We do some example printouts in the event 1 " << endl ;
      cout << "Neutralino Decay found at the vertex " << NeutDecVtx->barcode() <<" (barcode)" << endl ;
      vector<HepMC::GenParticle*> NeutChildren;
       
    // Get Neutralino Children 
    for ( HepMC::GenVertex::particles_out_const_iterator Neut_daughters =  NeutDecVtx->particles_out_const_begin(); 
	      Neut_daughters != NeutDecVtx->particles_out_const_end(); Neut_daughters++ ) 
      { 
         NeutChildren.push_back(*Neut_daughters);
      }
      
     cout << " Number of Neutralino (immediate) children = " << NeutChildren.size() << endl ;
      for (unsigned int ic=0; ic<NeutChildren.size(); ic++ )
      {
          NeutChildren[ic]->print() ;   
      }
  }
*/

   
/*   
for ( HepMC::GenVertex::particle_iterator
         des=HiggsDecVtx->particles_begin(HepMC::descendants);
	 des!=HiggsDecVtx->particles_end(HepMC::descendants); des++ )
   {
      if ( e.id().event() == 1 ) (*des)->print() ;
      if ( (*des)->status() == 1 ) StableHiggsDesc.push_back(*des) ;
   }
 */   
   
   // make 4-part inv.mass
   //
/*
   double px4, py4, pz4, e4;
   px4=py4=pz4=e4=0. ;
   if ( StableHiggsDesc.size() == 4 )
   {
       for ( unsigned int i=0; i<StableHiggsDesc.size(); i++ )
       {
          px4 += StableHiggsDesc[i]->momentum().px();
          py4 += StableHiggsDesc[i]->momentum().py();
          pz4 += StableHiggsDesc[i]->momentum().pz();
          e4  += StableHiggsDesc[i]->momentum().e();
       }
       XMass4part = HepMC::FourVector(px4,py4,pz4,e4).m() ;
       fHist4muMass->Fill( XMass4part ) ;
   }
   //cout << " 4-part inv. mass = " << XMass4part << endl ;
   // make 2-pairs (ZZ) inv.mass
   //
   //cout << " selected Z-candidates in this event : " << Mom2partCont.size() << endl ;
   for ( unsigned int i=0; i<Mom2partCont.size(); i++ )
   {
      for ( unsigned int j=i+1; j<Mom2partCont.size(); j++ )
      {
         // Mom2pairs = Mom2partCont[i] + Mom2partCont[j] ;
	 XMass2pairs = HepMC::FourVector((Mom2partCont[i].px()+Mom2partCont[j].px()),
	                                 (Mom2partCont[i].py()+Mom2partCont[j].py()),
					 (Mom2partCont[i].pz()+Mom2partCont[j].pz()),
					 (Mom2partCont[i].e() +Mom2partCont[j].e())).m() ;
	 fHistZZMass->Fill( XMass2pairs ) ;
         //cout << " 2-pairs (ZZ) inv. mass = " << XMass2pairs << endl ;
      }
   }
    */
   return ;
   
}

void GENAnalyzer::endRun( const edm::Run& r, const edm::EventSetup& )
{

   return;

}


void GENAnalyzer::endJob()
{
   
   return ;
}
 
typedef GENAnalyzer  MyGenAnalyzer;
DEFINE_FWK_MODULE(MyGenAnalyzer);
