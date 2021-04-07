#include "GeneratorInterface/genFilters/interface/ThreeMuonsSameOriginExtended.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

ThreeMuonsSameOriginExtended::ThreeMuonsSameOriginExtended(const edm::ParameterSet& iConfig) :
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
    edm::LogWarning("ThreeMuonsSameOriginExtended") << "WARNING: ThreeMuonsSameOriginExtended: size of PtMin, EtaMax, motherID, and/or Status does not match ParticleID size!" << std::endl;   
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

ThreeMuonsSameOriginExtended::~ThreeMuonsSameOriginExtended()
{
 
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to skim the data  ------------
bool ThreeMuonsSameOriginExtended::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByToken(token_, evt);
  
  totalEvents_++;
  
  // Problems with this filter: motherID_ is not taken into account; Values of other variables are defined inside the program and not taken from
  // the python file
  
  const HepMC::GenEvent * myGenEvent = evt->GetEvent();
  
  std::vector<float> pts; // Transverse momenta
  std::vector<float> phis;// signed values of phi
  std::vector<float> etas;// signed values of eta
  std::vector<int> charges;// indirectly calculated charges
  
  std::vector<float> es;// muon energy
  std::vector<float> p1;// muon px
  std::vector<float> p2;
  std::vector<float> p3;
  
  std::vector<int> muon_BC; //list of barcodes of the muons
  std::vector<std::vector<int>> ancestor_vtx;// 2D vector of barcodes of all ancestor vertices of each muon
  
  //std::vector<std::unique_ptr<HepMC::GenVertex>> parent_vtx_reference;// stores references to the production vertex of muon. Memory leak?

  //std::cout<<"Test Print out at event number "<<  totalEvents_  <<std::endl;
  

  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); // loop over particles
	p != myGenEvent->particles_end(); ++p ) {
        
     if(abs((*p)->pdg_id())==13 && ((*p)->status()==1)){// Used to collect all the muons (with decay status 1 / not decayed)
     
        pts.push_back((*p)->momentum().perp());
        etas.push_back((*p)->momentum().eta());
        phis.push_back((*p)->momentum().phi());
        charges.push_back(-1*((*p)->pdg_id()/13));
        
        muon_BC.push_back((*p)->barcode());
        
        //std::cout<<"Muon has Px "<<(*p)->momentum().px()<<" Py "<<(*p)->momentum().py()<<" Pz "<<(*p)->momentum().pz()<<" E "<<(*p)->momentum().e()<< std::endl;
        
        es.push_back((*p)->momentum().e());
        p1.push_back((*p)->momentum().px());
        p2.push_back((*p)->momentum().py());
        p3.push_back((*p)->momentum().pz());
        
        std::vector<int> vec_BC;// contains barcode of ancestor vertices of the muon including of prd. vertex
        for (HepMC::GenVertex::vertex_iterator anc=(*p)->production_vertex()->vertices_begin(HepMC::ancestors);anc !=(*p)->production_vertex()->vertices_end(HepMC::ancestors);++anc ){
            vec_BC.push_back((*anc)->barcode());
        }
        std::reverse(vec_BC.begin(),vec_BC.end());// Reverse the vector so that the closest vertex to the muon is at the beginning of the vector
        ancestor_vtx.push_back(vec_BC);
        
        //HepMC::GenVertex* vtx = (*p)->production_vertex();
        //parent_vtx_reference.push_back(std::unique_ptr<HepMC::GenVertex>((*p)->production_vertex())); Memory leak?
         
        //std::cout<<"The particle with barcode "<< (*p)->barcode()  <<" has pdgid "<< (*p)->pdg_id() <<std::endl;
        
        
     }//finish collection of a muon
  }//end gen particle loop
  if(pts.size()>=1){
    //std::cout <<"There are "<<pts.size()<<" muon(s) with max pT, "<<*max_element(std::begin(pts), std::end(pts))<<" ; min eta, "<<*min_element(std::begin(etas), std::end(etas))<<" ; and max energy, "<<*max_element(std::begin(es), std::end(es))<<std::endl;
  }
  if(pts.size()>=40){
    std::cout <<"No. of muons: "<<pts.size()<<std::endl;
  }
  
  int x(0);// stores number of unique muon combinations
  
  std::vector<bool> HLT_pass;//stores whether a particular combination of 3 muons would pass the HLT.
                             //Condition for track not applied. Inv Mass considered separately.
  std::vector<bool> L1_pass;//if the 3 muons would pass L1. Doesn't take into account events with only 2 muons
  std::vector<bool> Charge_pass;//if total charge of muons is +/- 1
  std::vector<bool> InvMass_pass;//Pass if invariant mass is between 1.60-2.02 GeV
  std::vector<bool> CommonVertex_pass;//if all 3 muons have a common vertex (say upto level 10 or higher)
  std::vector<bool> CommonVertex_Low_pass;//if all 3 muons have a common vertex (say upto level 4 or lower)
  std::vector<bool> CommonVertex2Mu_Low_pass;//if at least two of the three muons have a common vertex (level 2 or lower)
  std::vector<bool> dR_pass;// Check if all three pairs of muons are within maxDr_ of each other
  
  if(pts.size()>=3){// reserve space for pass vectors so that the code runs faster
      int reserve_size=((pts.size()*(pts.size()-1)*(pts.size()-2))/6);
      
      HLT_pass.reserve(reserve_size);
      L1_pass.reserve(reserve_size);
      Charge_pass.reserve(reserve_size);
      InvMass_pass.reserve(reserve_size);
      CommonVertex_pass.reserve(reserve_size);
      CommonVertex_Low_pass.reserve(reserve_size);
      CommonVertex2Mu_Low_pass.reserve(reserve_size);
      dR_pass.reserve(reserve_size);
      
  }
  
  IsL1 isl1;
  IsHLT ishlt;
  HasCommonVertex3Mu hascv3mu;
  HasCommonVertex2Mu hascv2mu;
  IsdR isdr;
  
  //This step looks at all unique combinations of muons and checks if they pass certain conditions
  
  if(pts.size()>=3&&pts.size()<=125){// need atleast 3 muons. The number of muons shouldn't be too big, it might crash? I've run with sizes of 120. It should be good as long as you don't print anything inside these loops.
    for (unsigned int A = 0; A < pts.size(); A++) {//all possible combinations of muons (indexed by A, B and C)
        for (unsigned int B = 0; B < pts.size(); B++) {
            for (unsigned int C = 0; C < pts.size(); C++) {
                if(A<B&&B<C){//All unique combinations of muons
                    
                    std::vector<unsigned int> mu_i = {A,B,C};
                    
                    // Check the HLT pTs
                    bool HLT(false);
                    for (int i = 0; i < 3; i++){//sends 3 pairs of muons to the IsHLT class
                        // (b + (a%b)) % b is used to get a positive value for negative modulo
                        if(ishlt(pts.at(mu_i.at((3 + ((i-1)%3)) % 3)),pts.at(mu_i.at((3 + ((i+1)%3)) % 3)))){
                            HLT=true;
                            break;
                        }
                    }
                    HLT_pass.push_back(HLT);
                    
                    
                    // Check the L1 pT and etas
                    bool L1T(false);
                    for (int i = 0; i < 3; i++){//sends 3 pairs of muons to the IsL1 class
                        // (b + (a%b)) % b is used to get a positive value for negative modulo
                        if(isl1(pts.at(mu_i.at((3 + ((i-1)%3)) % 3)),pts.at(mu_i.at((3 + ((i+1)%3)) % 3)),etas.at(mu_i.at((3 + ((i-1)%3)) % 3)),etas.at(mu_i.at((3 + ((i+1)%3)) % 3)))){
                            L1T=true;
                            break;
                        }
                    }
                    L1_pass.push_back(L1T);
                    
                    // Check dR
                    bool dRp(true);
                    for (int i = 0; i < 3; i++){//sends 3 pairs of muons to the IsdR class
                        // (b + (a%b)) % b is used to get a positive value for negative modulo
                        if(isdr(phis.at(mu_i.at((3 + ((i-1)%3)) % 3)),phis.at(mu_i.at((3 + ((i+1)%3)) % 3)),etas.at(mu_i.at((3 + ((i-1)%3)) % 3)),etas.at(mu_i.at((3 + ((i+1)%3)) % 3)),maxDr_)){
                            dRp=false;// if one pair has greater dR than maxdR, it doesn't pass
                            break;
                        }
                    }
                    dR_pass.push_back(dRp);
                    
                    
                    // Check the charges
                    Charge_pass.push_back(abs(charges.at(A)+charges.at(B)+charges.at(C))==1);
                    
                    
                    // Look at their invariant mass
                    TLorentzVector sum_3mu(p1.at(A)+p1.at(B)+p1.at(C),p2.at(A)+p2.at(B)+p2.at(C),p3.at(A)+p3.at(B)+p3.at(C),es.at(A)+es.at(B)+es.at(C));
                    InvMass_pass.push_back((sum_3mu.M()>1.5)&&(sum_3mu.M()<2.1));
                    
                    
                    // Check for common vertices (for 3 muons)
                    std::vector<int> vec_BC_A = ancestor_vtx[A];// contains barcode of ancestor vertices of the muon including of prd. vertex
                    std::vector<int> vec_BC_B = ancestor_vtx[B];
                    std::vector<int> vec_BC_C = ancestor_vtx[C];
                    CommonVertex_pass.push_back(hascv3mu(vec_BC_A,vec_BC_B,vec_BC_C,10));// last argument is the maximum vertex level. 1 equals parent.
                    CommonVertex_Low_pass.push_back(hascv3mu(vec_BC_A,vec_BC_B,vec_BC_C,4));// lower max vertex level (upto great-great grandparent vertex)
                    
                    
                    // Check for common vertices (for 2 of the three muons)
                    bool CV2mu(false);
                    unsigned int max_level(1); // Max number of ancestor vertices that are searched. 1=parent, 2=grand-parent, etc
                    if(hascv2mu(vec_BC_A,vec_BC_B,max_level)) CV2mu=true;
                    if(hascv2mu(vec_BC_A,vec_BC_C,max_level)) CV2mu=true;
                    if(hascv2mu(vec_BC_B,vec_BC_C,max_level)) CV2mu=true;
                    CommonVertex2Mu_Low_pass.push_back(CV2mu);
                    
                    x+=1;
                    /*
                    if(Charge_pass.at(Charge_pass.size()-1)&&InvMass_pass.at(Charge_pass.size()-1)&&L1_pass.at(Charge_pass.size()-1)){
                        std::cout<<"Barcodes of muons which pass: "<<muon_BC.at(A)<<muon_BC.at(B)<<muon_BC.at(C)<<std::endl;
                        std::cout<<"Whether 2mu passed: "<<CommonVertex2Mu_Low_pass.at(Charge_pass.size()-1)<<std::endl;
                    }
                    */
                    
                }// end unique muon loop
            }// C loop
        }// B loop
    }// A loop
  }// end pts.size()>=3
  
  
  bool Final_pass(false);
  // Final_pass checks whether there's atleast one combination of 3 muons that pass all these conditions
  // Poible checks: HLT_pass, L1_pass, Charge_pass, InvMass_pass, CommonVertex_pass, CommonVertex_Low_pass, CommonVertex2Mu_Low_pass, dR_pass
  for (unsigned int comb_idx = 0; comb_idx < Charge_pass.size(); comb_idx++) {//comb_idx indexes all unique combinations of muons
      if(HLT_pass.at(comb_idx)&&L1_pass.at(comb_idx)&&Charge_pass.at(comb_idx)&&InvMass_pass.at(comb_idx)&&CommonVertex_Low_pass.at(comb_idx)){
          Final_pass=true;
          break;
      }
  }
  
  HLT_pass.clear();
  L1_pass.clear();
  Charge_pass.clear();
  InvMass_pass.clear();
  CommonVertex_pass.clear();
  CommonVertex_Low_pass.clear();
  CommonVertex2Mu_Low_pass.clear();
  dR_pass.clear();
  
  if (Final_pass){
    std::cout<<"Three muons found at event no. "<<totalEvents_<<" satisfying all conditions."<<std::endl;
  }
  
  if (Final_pass) {
    
    passedEvents_++;
    return true;
  } else {
    return false;
  }
  
}// End bool filter class

// ------------ method called once each job just after ending the event loop  ------------
void ThreeMuonsSameOriginExtended::endJob() {
  edm::LogInfo("ThreeMuonsSameOriginExtended") << "=== Results of ThreeMuonsSameOriginExtended: passed "
                                        << passedEvents_ << "/" << totalEvents_ << " events" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeMuonsSameOriginExtended);
