#ifndef CustomThreeMuFilter_checkAncestor_dev_h
#define CustomThreeMuFilter_checkAncestor_dev_h

// -*- C++ -*-
//
// Package:    GeneratorInterface/genFilters
// Class:      genFilters
//
/**\class genFilters genFilters.cc GeneratorInterface/genFilters/plugins/genFilters.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vladimir Cherepanov
//         Created:  Thu, 05 Mar 2020 13:45:13 GMT
//
//



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TLorentzVector.h"

namespace edm {
  class HepMCProduct;
}

//
// class declaration
//

class CustomThreeMuFilter_checkAncestor_dev : public edm::EDFilter {
 public:
  explicit CustomThreeMuFilter_checkAncestor_dev(const edm::ParameterSet&);
  ~CustomThreeMuFilter_checkAncestor_dev() override;
  
 private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  bool isAncestor(HepMC::GenParticle* particle, int IDtoMatch) const;
  // ----------member data ---------------------------
  
  edm::InputTag src_;              // input tag
  edm::EDGetTokenT<edm::HepMCProduct> token_;
  int numRequired_;                // number of particles required to pass filter
  bool acceptMore_;                // if true (default), accept numRequired or more.
                                   // if false, accept events with exactly equal to numRequired.
  std::vector<int> particleID_;    // vector of particle IDs to look for
  // the four next variables can either be a vector of length 1 (in which case the same
  // value is used for all particle IDs) or of length equal to the length of ParticleID (in which
  // case the corresponding value is used for each).
  std::vector<int> motherID_;      // mother ID of particles (optional)
  std::vector<double> ptMin_;      // minimum Pt of particles
  std::vector<double> etaMax_;     // maximum fabs(eta) of particles
  std::vector<int> status_;        // status of particles
  float invMassMin_;               // min invariant mass of all particles
  float invMassMax_;               // max invariant mass of all particles
  float maxDr_;                    // max dR between any of the particles
  int totalEvents_;                // counters
  int passedEvents_;
};


#endif
