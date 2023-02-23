#ifndef ThreeMuonsSameOrigin_ztt_taue_Extended_h
#define ThreeMuonsSameOrigin_ztt_taue_Extended_h

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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TLorentzVector.h"
namespace edm {
  class HepMCProduct;
}

//
// class declaration
//

  class IsL1 {
  public:
      // returns true if dimuon L1 trigger satisfied
      bool operator()(float pt1, float pt2, float eta1, float eta2) { 
      if ( fabs(eta1-eta2)<1.81&&((pt1>9.9||pt2>9.9)||(fabs(eta1)<1.61&&fabs(eta2)<1.61)) ) return 1;
      return 0;
      }
  };
  
  class IsdR {
  public:
      // returns true if two muons have dR less than maxdR
      bool operator()(float phi1, float phi2, float eta1, float eta2, float maxdR) {
          float dphi = phi1 - phi2;
          if(dphi >= TMath::Pi()) dphi = dphi-2*TMath::Pi();// makes dphi between -pi and pi
          if(dphi <=-TMath::Pi()) dphi = dphi+2*TMath::Pi();
          float dR = sqrt(pow(dphi,2) + pow(eta1 - eta2,2));
          if ( dR >= maxdR ) return 1;
          return 0;
      }
  };
  
  class IsHLT {
  public:
      // returns true if HLT trigger satisfied
      bool operator()(float pt1, float pt2, float maxhlt1, float maxhlt2) { 
      if ( pt1>maxhlt1&&pt2>maxhlt2) return 1;
      return 0;
      }
  };
  
  class HasCommonVertex3Mu {
  public:
      // returns true if the 3 muons have a common ancestor vertex.
      // Inputs ancestor barcodes and the maximum ancestor level
      // max level 1 corresponds to parents, 2 corresponds to grand-parents, etc
      bool operator()(std::vector<int> vec_A, std::vector<int> vec_B, std::vector<int> vec_C, unsigned int max_l) {
  
          bool found(false);
  
          for (unsigned int A = 0; A < vec_A.size() &&(A<max_l); A++){// all possible combinations of barcodes
              for (unsigned int B = 0; B < vec_B.size() &&(B<max_l); B++){
                  for (unsigned int C = 0; C < vec_C.size() &&(C<max_l); C++){
                      if((vec_A.at(A)==vec_B.at(B))&&(vec_B.at(B)==vec_C.at(C))) found=true;
                  }
              }
          }
          if (found) return 1;
          return 0;
      }
  };
  
  class HasCommonVertex2Mu {
  public:
      // returns true if 2 of the 3 muons have a common ancestor vertex.
      // Inputs ancestor barcodes and the maximum ancestor level
      // max level 1 corresponds to parents, 2 corresponds to grand-parents, etc
      bool operator()(std::vector<int> vec_A, std::vector<int> vec_B, unsigned int max_l) {
  
          bool found(false);
  
          for (unsigned int A = 0; A < vec_A.size() &&(A<max_l); A++){// all possible combinations of barcodes
              for (unsigned int B = 0; B < vec_B.size() &&(B<max_l); B++){
                  if(vec_A.at(A)==vec_B.at(B)) found=true;
              }
          }
          if (found) return 1;
          return 0;
      }
  };

class ThreeMuonsSameOrigin_ztt_taue_Extended : public edm::EDFilter {
 public:
  explicit ThreeMuonsSameOrigin_ztt_taue_Extended(const edm::ParameterSet&);
  ~ThreeMuonsSameOrigin_ztt_taue_Extended() override;
  
 private:
  bool filter(edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  
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
