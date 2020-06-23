#include "GeneratorInterface/genFilters/interface/ThreeMuonsSameOrigin.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

ThreeMuonsSameOrigin::ThreeMuonsSameOrigin(const edm::ParameterSet& iConfig) :
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
    edm::LogWarning("ThreeMuonsSameOrigin") << "WARNING: ThreeMuonsSameOrigin: size of PtMin, EtaMax, motherID, and/or Status does not match ParticleID size!" << std::endl;   
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

ThreeMuonsSameOrigin::~ThreeMuonsSameOrigin()
{
 
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to skim the data  ------------
bool ThreeMuonsSameOrigin::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByToken(token_, evt);
  
  totalEvents_++;
  int nFound = 0;
  
  const HepMC::GenEvent * myGenEvent = evt->GetEvent();
  TLorentzVector sum(0,0,0,0);  
  std::vector<float> phis;
  std::vector<float> etas;
  std::vector<float> dR; //dR.push_back(0);
  std::vector<int> vtxBC;
  bool twomuons(false);
  std::vector<int> vtxBC3;
  bool threemuons(false);

  //std::cout<<"Test Print out at event number "<<  totalEvents_  <<std::endl;


  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); // loop over particles
	p != myGenEvent->particles_end(); ++p ) {
    

    //std::cout<<"What particles we have in the generated pythia event:  "<< (*p)->pdg_id()  << std::endl;
    //std::cout<<"The particle with barcode "<< (*p)->barcode()  <<" has pdgid "<< (*p)->pdg_id() <<std::endl;
    int twomuon_i(0);
    int threemuon_i(0);
    int VertexBC;
    int muoncount1;
    int muoncount2;

    for (unsigned int i = 0; i < particleID_.size(); ++i) {  // loop over targets
      if ((particleID_[i] == 0 || abs(particleID_[i]) == abs((*p)->pdg_id())) &&
	  (*p)->momentum().perp() > ptMin_[i] &&
	  fabs((*p)->momentum().eta()) < etaMax_[i] &&
	  (status_[i] == 0 || (*p)->status() == status_[i])) {
     
     
    //std::cout<<"Barcode of production vertex of muon is:  "<< (*p)->production_vertex()->barcode()  << std::endl;
    std::cout<<"Muon "<<i<< " has Px "<<(*p)->momentum().px()<<" Py "<<(*p)->momentum().py()<<" Pz "<<(*p)->momentum().pz()<<" E "<<(*p)->momentum().e()<< std::endl;
    
    
    if(!twomuons&&twomuon_i==0){
      std::cout<<"Barcode of production vertex of muon2 with barcode "<< (*p)->barcode()  << " is "<< (*p)->production_vertex()->barcode()  << std::endl;
      
      twomuons=std::find(vtxBC.begin(), vtxBC.end(), (*p)->production_vertex()->barcode()) != vtxBC.end();
      vtxBC.push_back((*p)->production_vertex()->barcode());
      twomuon_i++;
      
      //std::cout<<"Pt is:  "<< (*p)->momentum().perp() << " &eta is: "<< fabs((*p)->momentum().eta())  << std::endl;
      if(twomuons){
        std::cout<<"Test Print out at event number "<<  totalEvents_  <<std::endl;
        std::cout<<"Found something at event no. "<<totalEvents_<<std::endl;
      //  std::cout<<"Pt is:  "<< (*p)->momentum().perp() << " &eta is: "<< fabs((*p)->momentum().eta())  << std::endl;
        std::cout<<"Energy of muon is: "<<(*p)->momentum().e()<<std::endl;
      //  break;
      }

    }//end if statement
    
    
    if(!threemuons&&threemuon_i==0){
      for ( HepMC::GenVertex::particles_in_const_iterator otp = (*p)->production_vertex()->particles_in_const_begin(); otp != (*p)->production_vertex()->particles_in_const_end(); ++otp){
        VertexBC=(*otp)->production_vertex()->barcode();
        break;
	    }
      std::cout<<"Barcode of production vertex of muon3 with barcode "<< (*p)->barcode()  << " is "<< (*p)->production_vertex()->barcode()  << " and its grandparent vertex is "<< VertexBC << std::endl;
      
      muoncount1=std::count(vtxBC3.begin(), vtxBC3.end(), (*p)->production_vertex()->barcode());
      muoncount2=std::count(vtxBC3.begin(), vtxBC3.end(), VertexBC);
      std::cout<<"muoncount1 is: "<<muoncount1<<std::endl;
      std::cout<<"muoncount2 is: "<<muoncount2<<std::endl;
      vtxBC3.push_back(VertexBC);
      vtxBC3.push_back((*p)->production_vertex()->barcode());
      threemuon_i++;
      
      //std::cout<<"Pt is:  "<< (*p)->momentum().perp() << " &eta is: "<< fabs((*p)->momentum().eta())  << std::endl;
      if(muoncount1>=2||muoncount2>=2){
        std::cout<<"Test Print out at event number "<<  totalEvents_  <<std::endl;
        threemuons=true;
        std::cout<<"Found 3 muons & same vertex at event no. "<<totalEvents_<<std::endl;
      //  std::cout<<"Pt is:  "<< (*p)->momentum().perp() << " &eta is: "<< fabs((*p)->momentum().eta())  << std::endl;
        std::cout<<"Energy of muon is: "<<(*p)->momentum().e()<<std::endl;
      //  break;
      }

    }//end if statement


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
    


      }
    } // loop over targets


    //if (nFound >= numRequired_&&twomuons){
    //  std::cout<<"Two muons found with same vertex at "<<totalEvents_<<std::endl;
    //}
    
    if (nFound >= numRequired_&&threemuons){
      std::cout<<"Three muons found found at event no. "<<totalEvents_<<std::endl;
      break; // stop looking if we don't mind having more
    }
  } // loop over particles


  
  /*  std::cout<<" ---- "   << std::endl;
  for(unsigned int i=0; i  < phis.size()-1; i++){
    for(unsigned int j=i+1; j  < phis.size(); j++){
      float dphi = phis.at(i) - phis.at(j);
      if(dphi >= TMath::Pi()) dphi = dphi-2*TMath::Pi();
      if(dphi <=-TMath::Pi()) dphi = dphi+2*TMath::Pi();
      std::cout<<"  dphi " << dphi<< " phi1  "<< phis.at(i)<< "  phi2  "<<  phis.at(j)<< std::endl;
      std::cout<<"  deta " << etas.at(i) - etas.at(j) << " eta1  "<< etas.at(i)<< "  eta2  "<<  etas.at(j)<< std::endl;
      std::cout<<"  dR "<< sqrt(pow(dphi,2) + pow(etas.at(i) - etas.at(j),2)) << std::endl;
    }
  }
  */
  //  if (nFound == numRequired_)  std::cout<<" numFound:  "<< nFound<< "  dR size   " << dR.size() <<std::endl;
  //  if (nFound == numRequired_)  for(auto &l:dR){std::cout<<" dR  "<< l <<std::endl;}
  
  if (nFound >= numRequired_&&threemuons) {
    
    //    std::cout<<"sum: "<< sum.M() <<std::endl;
    //    sum.Print();
    passedEvents_++;
    return true;
  } else {
    return false;
  }
  
}

// ------------ method called once each job just after ending the event loop  ------------
void ThreeMuonsSameOrigin::endJob() {
  edm::LogInfo("ThreeMuonsSameOrigin") << "=== Results of ThreeMuonsSameOrigin: passed "
                                        << passedEvents_ << "/" << totalEvents_ << " events" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreeMuonsSameOrigin);
