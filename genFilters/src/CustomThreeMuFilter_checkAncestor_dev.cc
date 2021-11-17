#include "GeneratorInterface/genFilters/interface/CustomThreeMuFilter_checkAncestor_dev.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

CustomThreeMuFilter_checkAncestor_dev::CustomThreeMuFilter_checkAncestor_dev(const edm::ParameterSet& iConfig) :
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
    edm::LogWarning("CustomThreeMuFilter_checkAncestor_dev") << "WARNING: CustomThreeMuFilter_checkAncestor_dev: size of PtMin, EtaMax, motherID, and/or Status does not match ParticleID size!" << std::endl;   
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

CustomThreeMuFilter_checkAncestor_dev::~CustomThreeMuFilter_checkAncestor_dev()
{
 
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}

bool CustomThreeMuFilter_checkAncestor_dev::isAncestor(HepMC::GenParticle* particle, int IDtoMatch) const {
  for (HepMC::GenVertex::particle_iterator ancestor = particle->production_vertex()->particles_begin(HepMC::ancestors);
       ancestor != particle->production_vertex()->particles_end(HepMC::ancestors);
       ++ancestor) {
     std::cout << __LINE__ << "]\t particle's PDG ID " << particle->pdg_id()
                           << " \t particle's ancestor's PDG ID " << (*ancestor)->pdg_id()
                           << " \t ID to match " << IDtoMatch << std::endl;

    if (abs((*ancestor)->pdg_id()) == abs(IDtoMatch)) {
        std::cout << __LINE__ << "]\t found!" << std::endl;
      return true;
    }
  }

   std::cout << __LINE__ << "]\t nope, no luck" << std::endl;
  return false;
}



// ------------ method called to skim the data  ------------
bool CustomThreeMuFilter_checkAncestor_dev::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByToken(token_, evt);
  
  totalEvents_++;
  int nFound = 0;
  bool accepted(false);
  const HepMC::GenEvent * myGenEvent = evt->GetEvent();
  TLorentzVector sum(0,0,0,0);  

  std::vector<HepMC::GenParticle*>  negMuons;
  std::vector<HepMC::GenParticle*>  posMuons;

  std::vector<float> phis;
  std::vector<float> etas;
  std::vector<float> dR; //dR.push_back(0);


  std::cout<<"  -----   "<< std::endl;
  //  myGenEvent->print();


  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); // loop over particles
	p != myGenEvent->particles_end(); ++p ) {
    
    for (unsigned int i = 0; i < particleID_.size(); ++i) {  // loop over targets
      if ((particleID_[i] == 0 || abs(particleID_[i]) == abs((*p)->pdg_id())) &&
	  (*p)->momentum().perp() > ptMin_[i] &&
	  fabs((*p)->momentum().eta()) < etaMax_[i] &&
	  (status_[i] == 0 || (*p)->status() == status_[i])) 
	{
	

	
	  if((*p)->pdg_id() ==  13) negMuons.push_back(*p);
	  if((*p)->pdg_id() == -13) posMuons.push_back(*p);




	  if(posMuons.size()!=0 && negMuons.size()!=0    &&  ( negMuons.size() + posMuons.size() ) > 2)
	    {

	      for(unsigned int iN = 0; iN < negMuons.size(); iN++){
		float tot_px(0.);
		float tot_py(0.);
		float tot_pz(0.);
		float tot_e(0.);
		float invmass(0.);
		isAncestor(negMuons.at(iN), 431);
		tot_px = negMuons.at(iN)->momentum().px();
		tot_py = negMuons.at(iN)->momentum().py();
		tot_pz = negMuons.at(iN)->momentum().pz();
		tot_e  = negMuons.at(iN)->momentum().e();
		for(unsigned int iP1 = 0; iP1 < posMuons.size()-1; iP1++){
		  isAncestor(posMuons.at(iP1), 431);
		  tot_px += posMuons.at(iP1)->momentum().px();
		  tot_py += posMuons.at(iP1)->momentum().py();
		  tot_pz += posMuons.at(iP1)->momentum().pz();
		  tot_e  += posMuons.at(iP1)->momentum().e();
		  for(unsigned int iP2 = iP1 + 1; iP2 < posMuons.size(); iP2++){
		    isAncestor(posMuons.at(iP2), 431);
		    tot_px += posMuons.at(iP2)->momentum().px();
		    tot_py += posMuons.at(iP2)->momentum().py();
		    tot_pz += posMuons.at(iP2)->momentum().pz();
		    tot_e  += posMuons.at(iP2)->momentum().e();
		    invmass = sqrt(tot_e * tot_e - tot_px * tot_px - tot_py * tot_py - tot_pz * tot_pz);

		    if(invmass > invMassMin_ && invmass < invMassMax_){std::cout<<" P triplet  "<< invmass<< std::endl; accepted = true; }

		  }
		}
	      }


	      //search for Mtriplet
	      for(unsigned int iP = 0; iP < posMuons.size(); iP++){
		float tot_px(0.);
		float tot_py(0.);
		float tot_pz(0.);
		float tot_e(0.);
		float invmass(0.);
                tot_px = posMuons.at(iP)->momentum().px();
                tot_py = posMuons.at(iP)->momentum().py();
                tot_pz = posMuons.at(iP)->momentum().pz();
                tot_e  = posMuons.at(iP)->momentum().e();

		for(unsigned int iN1 = 0; iN1 < negMuons.size()-1; iN1++){
                  tot_px += negMuons.at(iN1)->momentum().px();
                  tot_py += negMuons.at(iN1)->momentum().py();
                  tot_pz += negMuons.at(iN1)->momentum().pz();
                  tot_e  += negMuons.at(iN1)->momentum().e();

		  for(unsigned int iN2 = iN1 + 1; iN2 < negMuons.size(); iN2++){
		    tot_px += negMuons.at(iN2)->momentum().px();
		    tot_py += negMuons.at(iN2)->momentum().py();
		    tot_pz += negMuons.at(iN2)->momentum().pz();
		    tot_e  += negMuons.at(iN2)->momentum().e();
		    

		    invmass = sqrt(tot_e * tot_e - tot_px * tot_px - tot_py * tot_py - tot_pz * tot_pz);
		    if(invmass > invMassMin_ && invmass < invMassMax_){std::cout<<" N triplet  "<< invmass<< std::endl; accepted = true; }

		  }
		}
	      }
	      std::cout<<"  sizes of neg and pos   "<< negMuons.size() << "   "<< posMuons.size() <<  "  accepted   "  << accepted<<std::endl;
	    }


	  
	  if(accepted) break; // find satisfiyng condition once
	  



	
	  /*

	    if(motherID_[i] == 0 ){ // do not check for mother ID if not sepcified
	    TLorentzVector par;par.SetPxPyPzE((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
	    
	    for(unsigned int k=0; k < phis.size(); k++){
	    float dphi = phis.at(k) - (*p)->momentum().phi();
	    if(dphi >= TMath::Pi()) dphi = dphi-2*TMath::Pi();
	    if(dphi <=-TMath::Pi()) dphi = dphi+2*TMath::Pi();
	    dR.push_back(sqrt(pow(dphi,2) + pow(etas.at(k) - (*p)->momentum().eta(),2)));
	    }
	    
	    bool passdr(true);
	    for(auto &l:dR){if(l>=maxDr_){passdr=false;break;}}
	    
	    if(passdr){
	    int size = phis.size();
	    if(size>=numRequired_-1){
	    if( (sum+par).M() > invMassMin_  && (sum+par).M() < invMassMax_ ){
	    phis.push_back((*p)->momentum().phi());
	    etas.push_back((*p)->momentum().eta());
	    sum+=par;
	    nFound++;
	    }
	    }else{
	    phis.push_back((*p)->momentum().phi());
	    etas.push_back((*p)->momentum().eta());
	    sum+=par;
	    nFound++;
	    }
	    
	    
	    break; // only match a given particle once!
	    }
	    }
	    
	    //  ---  Modified in case of one wants to control the Mo particle of muons
	    else{
	    bool hascorrectmother=false;
	    for ( HepMC::GenVertex::particles_in_const_iterator mo = (*p)->production_vertex()->particles_in_const_begin(); mo != (*p)->production_vertex()->particles_in_const_end(); ++mo){
	    if( (*mo)->pdg_id() == motherID_[i]){
	    hascorrectmother = true;
	    break;
	    }
	    }
	    if(hascorrectmother){
	    nFound++;
	    break; // only match a given particle once! 
	    }
	    }
	  */
	}




    } // loop over targets





    //    if (nFound == numRequired_) break;  // stop looking if we don't mind having more
  } // loop over particles


  if(accepted)
    return true;
  return false;


}
/*  if (nFound == numRequired_) {
    passedEvents_++;
    return true;
  } else {
    return false;
  }
*/  


// ------------ method called once each job just after ending the event loop  ------------
void CustomThreeMuFilter_checkAncestor_dev::endJob() {
  edm::LogInfo("CustomThreeMuFilter_checkAncestor_dev") << "=== Results of CustomThreeMuFilter_checkAncestor_dev: passed "
                                        << passedEvents_ << "/" << totalEvents_ << " events" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CustomThreeMuFilter_checkAncestor_dev);
