#include "GeneratorInterface/genFilters/interface/CustomThreeParticleFilter.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

CustomThreeParticleFilter::CustomThreeParticleFilter(const edm::ParameterSet& iConfig) :
  src_(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag(std::string("generator"),"unsmeared"))),
  token_(consumes<edm::HepMCProduct>(src_)),
  numRequired_(iConfig.getParameter<int>("NumRequired")),
  particleID_(iConfig.getParameter< std::vector<int> >("ParticleID")),
  ptMin_(iConfig.getParameter< std::vector<double> >("PtMin")),
  etaMax_(iConfig.getParameter< std::vector<double> >("EtaMax")),
  status_(iConfig.getParameter< std::vector<int> >("Status")),
  invMassMin_(iConfig.getParameter<double>("invMassMin")),
  invMassMax_(iConfig.getParameter<double>("invMassMax")),
  maxDr_(iConfig.getParameter<double>("maxDr")),

  totalEvents_(0), passedEvents_(0)
{
  //here do whatever other initialization is needed

  // default pt, eta, status cuts to "don't care"
  std::vector<double> defptmin(1, 0);
  std::vector<double> defetamax(1, 999.0);
  std::vector<int> defstat(1, 0);
  std::vector<int> defmother;
  //  double defmassmin(0);
  //  double defmassmax(1000);
 
  defmother.push_back(0);
  motherID_ = iConfig.getUntrackedParameter< std::vector<int> >("MotherID", defstat);

  // check for same size
  if ( (ptMin_.size() > 1 &&  particleID_.size() != ptMin_.size()) 
       ||  (etaMax_.size() > 1 && particleID_.size() != etaMax_.size()) 
       ||  (status_.size() > 1 && particleID_.size() != status_.size()) 
       ||  (motherID_.size() > 1 && particleID_.size() != motherID_.size())
       ) {
    edm::LogWarning("CustomThreeParticleFilter") << "WARNING: CustomThreeParticleFilter: size of PtMin, EtaMax, motherID, and/or Status does not match ParticleID size!" << std::endl;   
  }
  
  // Fill arrays with defaults if necessary
  while (ptMin_.size() < particleID_.size())
    ptMin_.push_back(defptmin[0]);
  while (etaMax_.size() < particleID_.size())
    etaMax_.push_back(defetamax[0]);
  while (status_.size() < particleID_.size())
    status_.push_back(defstat[0]);
  while (motherID_.size() < particleID_.size())
    motherID_.push_back(defmother[0]);
}

CustomThreeParticleFilter::~CustomThreeParticleFilter()
{
 
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to skim the data  ------------
bool CustomThreeParticleFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
        
        
        edm::Handle<edm::HepMCProduct> evt;
        iEvent.getByToken(token_, evt);
          
        totalEvents_++;
        
        //int nFound = 0;
        bool WhetherThreeMuonPosPass(false);
        bool WhetherThreeMuonNegPass(false);
        bool WhetherTwoMuonandPionPosPass(false);
        bool WhetherTwoMuonandPionNegPass(false);
        
        std::vector<int> posMuons;
        std::vector<int> negMuons;
        std::vector<int> posPionsOrKaons;
        std::vector<int> negPionsOrKaons;
        std::vector<int> OppositeSide;
        
        const HepMC::GenEvent* myGenEvent = evt->GetEvent();
        
        // First loop: categorize particles
        for (HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p) {
          if ((*p)->status() <= 0) continue;
          int pdgId = (*p)->pdg_id();
          float pt = (*p)->momentum().perp();
          float eta = (*p)->momentum().eta();
        
          if (pt < 1.0 || fabs(eta) > 3.1) continue;
        
          if (pdgId == 13 && fabs(eta) < 2.9) negMuons.push_back((*p)->barcode());
          else if (pdgId == -13 && fabs(eta) < 2.9) posMuons.push_back((*p)->barcode());
          else if ((pdgId == 211 || pdgId == 321) && fabs(eta) < 2.9) negPionsOrKaons.push_back((*p)->barcode());
          else if ((pdgId == -211 || pdgId == -321) && fabs(eta) < 2.9) posPionsOrKaons.push_back((*p)->barcode());
          else if (pdgId == 15 || pdgId == -15 || pdgId == 13 || pdgId == -13 || pdgId == 11 || pdgId == -11 ||
                   abs(pdgId) <= 5 || pdgId == 21) {
            OppositeSide.push_back((*p)->barcode());
          }
        }
        
        // Now do matching
        if (posMuons.empty() || negMuons.empty()) return false;
        
        // Helper function to fetch TLorentzVector for a given barcode
        auto getLorentz = [&](int barcode) -> TLorentzVector {
          for (auto p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p) {
            if ((*p)->barcode() == barcode) {
              return TLorentzVector((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
            }
          }
          return TLorentzVector(0,0,0,0); // default if not found (shouldn't happen)
        };
        
        // Triplet search: muon+muon+muon
        for (unsigned int iN = 0; iN < negMuons.size(); ++iN) {
          for (unsigned int iP1 = 0; iP1 < posMuons.size(); ++iP1) {
            for (unsigned int iP2 = 0; iP2 < iP1; ++iP2) {
              TLorentzVector massTriplet = getLorentz(negMuons[iN]) + getLorentz(posMuons[iP1]) + getLorentz(posMuons[iP2]);
              if (massTriplet.M() > invMassMin_ && massTriplet.M() < invMassMax_ && massTriplet.Pt() > 14.5) {
                WhetherThreeMuonPosPass = true;
              }
            }
          }
        }
        
        for (unsigned int iP = 0; iP < posMuons.size(); ++iP) {
          for (unsigned int iN1 = 0; iN1 < negMuons.size(); ++iN1) {
            for (unsigned int iN2 = 0; iN2 < iN1; ++iN2) {
              TLorentzVector massTriplet = getLorentz(posMuons[iP]) + getLorentz(negMuons[iN1]) + getLorentz(negMuons[iN2]);
              if (massTriplet.M() > invMassMin_ && massTriplet.M() < invMassMax_ && massTriplet.Pt() > 14.5) {
                WhetherThreeMuonNegPass = true;
              }
            }
          }
        }
        
        // Triplet search: muon+muon+pion/kaon
        if (!posPionsOrKaons.empty()) {
          for (unsigned int iN = 0; iN < negMuons.size(); ++iN) {
            for (unsigned int iP1 = 0; iP1 < posMuons.size(); ++iP1) {
              for (unsigned int iP2 = 0; iP2 < posPionsOrKaons.size(); ++iP2) {
                TLorentzVector massTriplet = getLorentz(negMuons[iN]) + getLorentz(posMuons[iP1]) + getLorentz(posPionsOrKaons[iP2]);
                if (massTriplet.M() > invMassMin_ && massTriplet.M() < invMassMax_ && massTriplet.Pt() > 14.5) {
                  WhetherTwoMuonandPionPosPass = true;
                }
              }
            }
          }
        }
        
        if (!negPionsOrKaons.empty()) {
          for (unsigned int iP = 0; iP < posMuons.size(); ++iP) {
            for (unsigned int iN1 = 0; iN1 < negMuons.size(); ++iN1) {
              for (unsigned int iN2 = 0; iN2 < negPionsOrKaons.size(); ++iN2) {
                TLorentzVector massTriplet = getLorentz(posMuons[iP]) + getLorentz(negMuons[iN1]) + getLorentz(negPionsOrKaons[iN2]);
                if (massTriplet.M() > invMassMin_ && massTriplet.M() < invMassMax_ && massTriplet.Pt() > 14.5) {
                  WhetherTwoMuonandPionNegPass = true;
                }
              }
            }
          }
        }
        
        // Final decision
        if (WhetherTwoMuonandPionPosPass || WhetherTwoMuonandPionNegPass || WhetherThreeMuonPosPass || WhetherThreeMuonNegPass) {
          passedEvents_++;
          //std::cout << "Something passed." << std::endl;
          //std::cout << "WhetherThreeMuonPosPass: " << WhetherThreeMuonPosPass
          //          << " WhetherThreeMuonNegPass: " << WhetherThreeMuonNegPass
          //          << " WhetherTwoMuonandPionPosPass: " << WhetherTwoMuonandPionPosPass
          //          << " WhetherTwoMuonandPionNegPass: " << WhetherTwoMuonandPionNegPass << std::endl;
          return true;
        }
        else {
          return false;
        }

  
}

// ------------ method called once each job just after ending the event loop  ------------
void CustomThreeParticleFilter::endJob() {
  edm::LogInfo("CustomThreeParticleFilter") << "=== Results of CustomThreeParticleFilter: passed "
                                        << passedEvents_ << "/" << totalEvents_ << " events" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CustomThreeParticleFilter);
