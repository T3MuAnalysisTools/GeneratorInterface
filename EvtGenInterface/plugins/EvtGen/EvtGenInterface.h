// Class Based on EvtGenInterface(LHC).
//
// Created March 2014
//
// This class is a modification of the original EvtGenInterface which was developed for EvtGenLHC 9.1.
// The modifications for EvtGen 1.3.0 are implemented by Ian M. Nugent
// I would like to thank the EvtGen developers, in particular John Black, and Mikhail Kirsanov for their assistance.
//


#ifndef gen_EvtGenInterface_h
#define gen_EvtGenInterface_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <memory>
#include <string>
#include <vector>

#include "EvtGenBase/EvtParticle.hh"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "GeneratorInterface/EvtGenInterface/interface/EvtGenInterfaceBase.h"

class myEvtRandomEngine;

namespace HepMC {
  class GenParticle;
  class GenEvent;
}				

class EvtId; 
class EvtGen;

namespace gen {
   class EvtGenInterface : public EvtGenInterfaceBase {
     public:
    // ctor & dtor
    EvtGenInterface( const edm::ParameterSet& );
    ~EvtGenInterface() override;
    
    void init() override;
    const std::vector<int>& operatesOnParticles() override { return m_PDGs; }      
    HepMC::GenEvent* decay( HepMC::GenEvent* ) override;
    void setRandomEngine(CLHEP::HepRandomEngine* v) override;
    static double flat();
    
  private:
    bool addToHepMC(HepMC::GenParticle* partHep,const EvtId &idEvt, HepMC::GenEvent* theEvent, bool del_daug); 
    void update_particles(HepMC::GenParticle* partHep,HepMC::GenEvent* theEvent,HepMC::GenParticle* p);
    void SetDefault_m_PDGs();
    bool findLastinChain(HepMC::GenParticle* &p);    
    bool hasnoDaughter(HepMC::GenParticle* p);
    void go_through_daughters(EvtParticle* part);
    bool find_decay(HepMC::GenParticle* sourceparticle, int pdgid);
    bool filter_acceptance(std::vector<HepMC::GenParticle*>  particles);
    bool filter_acceptance2(std::vector<TLorentzVector> particles);
    std::vector<TLorentzVector> SortedPtMuons(std::vector<TLorentzVector>  p);
    bool CheckEvtParticle(EvtParticle* p);
    EvtGen *m_EvtGen;                // EvtGen main  object

    std::vector<EvtId> forced_id;     // EvtGen Id's of particles  which are to be forced by EvtGen
    std::vector<int> forced_pdgids;    // PDG Id's of particles which are to be forced by EvtGen
    
    std::vector<int> ignore_pdgids;  // HepId's of particles  which are to be ignroed by EvtGen
    
    // Adding parameters for polarization of spin-1/2 particles
    std::vector<int> polarize_ids;
    std::vector<double> polarize_pol;
    std::map<int, float> polarizations;
    std::vector<int> forced_channels;
    std::vector<std::string> forced_parent_meson_name;
    int nredecays_of_parent_particle;
    int min_inv_mass;
    int max_inv_mass;

    int BmixingOption = 1;        
    edm::ParameterSet* fPSet;

    static CLHEP::HepRandomEngine* fRandomEngine;
    myEvtRandomEngine* the_engine;
  };
}
#endif
