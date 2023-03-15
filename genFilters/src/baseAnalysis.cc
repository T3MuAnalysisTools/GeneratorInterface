#include <iostream>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/SimpleVector.h"
#include "HepPDT/ParticleData.hh"
#include "CLHEP/Vector/LorentzVector.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

using namespace std;

// ROOT headers
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"

//base analysis header
#include "GeneratorInterface/genFilters/interface/baseAnalysis.h"

//**********************

/**
Constructor
*/
baseAnalysis::baseAnalysis()
{ 
  m_Jetfinder_enabled=true;
}

/**
Destructor: delete all the Pointer
*/
baseAnalysis::~baseAnalysis()
{
  for (std::vector<TH1D*>::const_iterator i(m_histVector.begin()); i != m_histVector.end(); ++i) {
    delete (*i);
  }
  m_histVector.resize(0);

  //probably not needed anymore, due to averagedHistograms function
  //   while (!m_histVectorVariableBin.empty()){
  //     delete m_histVectorVariableBin.back();
  //     m_histVectorVariableBin.pop_back();
  //   }
}

/**
InitJetFinder: Initialisation of JetFinder
*/
int baseAnalysis::InitJetFinder(double coneRadius, double overlapThreshold, double jet_ptmin, double lepton_ptmin, double DeltaR_lepton_track)
{

  //I removed the delete
  m_jetDef=0;
  m_plugin=0;
    
 // initialise fastjet
  m_coneRadius = coneRadius;
  m_overlapThreshold = overlapThreshold;
  m_jet_ptmin = jet_ptmin;
  m_lepton_ptmin = lepton_ptmin;
  m_DeltaR_lepton_track = DeltaR_lepton_track;

  m_plugin = new fastjet::SISConePlugin(m_coneRadius, m_overlapThreshold);
  m_jetDef = new fastjet::JetDefinition(m_plugin);

  if(!m_jetDef) return false;

  m_Jetfinder_enabled=true;
  return true;

}

/**
FindJetableParticles
*/
int baseAnalysis::FindJetableParticles(const HepMC::GenEvent* hepmcevt)
{

  //delete m_clust_seq;
  //m_clust_seq=0;
  if(!m_Jetfinder_enabled || !m_jetDef){
    std::cout << "baseAnalysis: FastJet Algorithm not initialized!" << std::endl;
    return false;
  }
  CLHEP::HepLorentzVector lv_leadingJet(0,0,0,0);

  for ( HepMC::GenEvent::particle_const_iterator p =
          hepmcevt->particles_begin(); p != hepmcevt->particles_end(); ++p ){

    // use the pdg id to identify the particles
    int pid = (*p)->pdg_id();


    
    // first of all we need stable particle
    if (!(((*p)->status()%1000) == 1)) continue;
    if ((*p)->end_vertex()) continue;
    
    // cut out the neutral particles (charge is not stored in HepMC)
    // just remove the neutrinos and Neutron from the final state
    if(IsNeutrino(pid) || IsNeutron(pid)) continue;

    // cut out all isolated leptons
    // variables for DeltaR calculation
    double deltaPhi = 0;
    double deltaEta = 0;
    double DeltaR = 9999;

    // cut variables
    double DeltaRCut = m_DeltaR_lepton_track;
    double LeptonpTCut = m_lepton_ptmin;

    if(IsElectron(pid) || IsMuon(pid)) {
      int lepton_barcode = (*p)->barcode();
      HepMC::GenParticle *lepton = hepmcevt->barcode_to_particle(lepton_barcode);
      
      // set pt cut for lepton
      if(fabs((*lepton).momentum().perp()) > LeptonpTCut) continue;
      
      // check if it is an isolated lepton
      for ( HepMC::GenEvent::particle_const_iterator q = hepmcevt->particles_begin(); 
            q != hepmcevt->particles_end(); ++q ){
        if(lepton_barcode==(*q)->barcode()) continue;

        // set DeltaR cut
        deltaPhi = (*p)->momentum().phi()-(*q)->momentum().phi();

        if( deltaPhi > TMath::Pi() ||  deltaPhi == TMath::Pi() )
          {
            deltaPhi = deltaPhi - 2 * TMath::Pi() ;
          }
        
        if( deltaPhi < - TMath::Pi() ||  deltaPhi == - TMath::Pi() )
          {
            deltaPhi = deltaPhi + 2 * TMath::Pi();
          }

        deltaEta = (*p)->momentum().eta()-(*q)->momentum().eta();
        double DeltaRnew = sqrt( pow(deltaPhi,2) + pow(deltaEta,2));

        if(DeltaRnew < DeltaR) DeltaR = DeltaRnew;
      }
      if(DeltaR > DeltaRCut) continue;
    }

    // or ... remove all charged particles 
    //if(chargedParticle(pid)==NOTCHARGED) continue;
    
    // They need to have a production Vertex --> we will not select incoming Protons from the Beam
    if (!(*p)->production_vertex()) continue;
    // And since they are stable they should not have a decay vertex.
    if ((*p)->end_vertex()) continue;
    
//     // check wether the track directly comes from a Z-Boson
//     if( trackfromPID(23,(*p)) ) continue;
    
//     // check wether the track directly comes from a W-Boson
//     if( trackfromPID(24,(*p)) ) continue;
//     if( trackfromPID(-24,(*p)) ) continue;
    
    //set the eta cut for charged particles 
    double eta = (*p)->momentum().eta(); 
    if(std::abs(eta) > m_max_eta )  continue;
    
    //set the pt cut for charged particles 
    double pt = (*p)->momentum().perp(); 
    if(pt < m_min_pt) continue;
    
    double px = (*p)->momentum().x();
    double py = (*p)->momentum().y();
    double pz = (*p)->momentum().z();
    double pe = (*p)->momentum().e();
    
    
    // jet finding in stable particles record for charged particles
    fastjet::PseudoJet tempVec;
    tempVec = fastjet::PseudoJet(px,py,pz,pe);
    m_input_particles.push_back(tempVec);
  }
  return true;
}

/**
FindJet: run JetFinder
*/
int baseAnalysis::FindJet(const HepMC::GenEvent* hepmcevt) {

  if(!m_Jetfinder_enabled || !m_jetDef){
    std::cout << "baseAnalysis: FastJet Algorithm not initialized!" << std::endl;
    return false;
  }
  
  // run the jet finding
  if(m_input_particles.empty()) {
    int input = FindJetableParticles(hepmcevt); 
    if(!input) 
      return false;
  }
  
  if(m_input_particles.size()) 
    m_clust_seq = new fastjet::ClusterSequence(m_input_particles,*m_jetDef);
  else
    m_clust_seq = 0;
  
  if(m_input_particles.empty() || m_clust_seq==0) 
    return false;
  
  // select jets above a threshold, e.g. 0.5 GeV                                                                              
  m_inclusive_jets= sorted_by_pt(m_clust_seq->inclusive_jets(m_jet_ptmin));
  
  return true;
}

vector<fastjet::PseudoJet> baseAnalysis::GetJet(const HepMC::GenEvent* hepmcevt) {

  if(m_inclusive_jets.empty()) {
    FindJet(hepmcevt);
    //int inputJet = FindJet(hepmcevt);
    //if(inputJet==FAILURE) std::cout<<"baseAnalysis::no jets found"<<std::endl;
  }
  return m_inclusive_jets;
  //return SUCCESS;
}

/**
DeleteJetObject: delete all the jet objects
*/
int baseAnalysis::DeleteJetObject(const HepMC::GenEvent* hepmcevt) {
  //if(!m_input_particles.size()&&!m_inclusive_jets.size()&&!m_clust_seq) return FAILURE;
  //if(!(m_input_particles.size()&&m_inclusive_jets.size()&&m_clust_seq)) return FAILURE;

  m_input_particles.clear();
  m_inclusive_jets.clear();
  delete m_clust_seq;
  m_clust_seq=0;

  return true;
}

/**
FindMissingEt: calculate the missing transversal energy of each event
*/
int baseAnalysis::FindMissingEt(const HepMC::GenEvent* hepmcevt) {

  int properStatus = 1; //stable particles
  
  exMissTruth = 0.0;
  eyMissTruth = 0.0;
  etMissTruth = 0.0;
  etsumMissTruth = 0.0;
  
  // loop over the particles and select pdgid and pt
  for ( HepMC::GenEvent::particle_const_iterator p =  hepmcevt->particles_begin(); p != hepmcevt->particles_end(); ++p ){

    int pid = (*p)->pdg_id();

    // skip all not stable particles
    if ( (*p)->status()%1000 != properStatus ) continue;
    if ((*p)->end_vertex()) continue;

    // fiducial range eta cut on 2.5
    if(fabs((*p)->momentum().eta()) > 2.5) continue;
    // minimum pt for tracks on 0.5 GeV
    if((*p)->momentum().perp() < 0.5) continue;

    double px  = (*p)->momentum().px();
    double py  = (*p)->momentum().py();
    double pt  = (*p)->momentum().perp();

    if ( IsNeutrino(pid) || abs(pid)==1000039 || abs(pid) ==1000022 || abs(pid) ==39 || abs(pid) ==5100022 ) {
      exMissTruth += px;    
      eyMissTruth += py;
      etsumMissTruth += pt;
    }
  }
  
  etMissTruth = sqrt(exMissTruth*exMissTruth + eyMissTruth*eyMissTruth); 
  
  return true;
}

/**
ClearMissingEt: setting the variables to calculate missing energy to zero
*/
int baseAnalysis::ClearMissingEt(const HepMC::GenEvent* hepmcevt) {
  exMissTruth = 0.0;
  eyMissTruth = 0.0;
  etMissTruth = 0.0;
  etsumMissTruth = 0.0;

  return true;
}

/** 
Check if final state particle: returns 0 if it is not a final state particle, 1 otherwise
*/
int baseAnalysis::IsFinalStateParticle(HepMC::GenParticle *p) {

  // first of all we need stable particle
  if (!((*p).status()%1000 == 1)) return false;
  // And since they are stable they should not have a decay vertex. 
  if ((*p).end_vertex()) return false;
    
    // fiducial range eta cut on 2.5
    if(fabs((*p).momentum().eta()) > 2.5) return false;
    // minimum pt for tracks on 0.5
    if((*p).momentum().perp() < 0.5) return false;

    return true;
}


/**
ClearEvent: delete and clear some variables from the event
*/
int baseAnalysis::ClearEvent(const HepMC::GenEvent* hepmcevt) {

  baseAnalysis::DeleteJetObject(hepmcevt);
  baseAnalysis::ClearMissingEt(hepmcevt);

  return true;
}

/**
In the final step all the histogramms are stored in a rootfile.
The name of the rootfile can be set with the function setOutpuFileName(const char* filename). 
 */
int baseAnalysis::Finalize() 
{
  //write the output in a root file
  TFile f(m_outputFileName.c_str(),"RECREATE");
  Finalize(&f);
  f.Close();
  
  return true;
}

/**
In the final step all the histogramms are stored in a rootfile.
The name of the rootfile can be set with the function setOutpuFileName(const char* filename). 
 */
int baseAnalysis::Finalize(TFile* output) 
{
  TDirectory* dir = output->GetDirectory(m_outputRootDir.c_str());
  if (dir) {
    dir->cd();
  } else {
    std::cout << "The directory " << m_outputRootDir << " will be created" << std::endl;
    dir = output->mkdir(m_outputRootDir.c_str());
    dir->cd();
  }

  // loop over the standard vector
  for ( vector< TH1D * >::const_iterator iter = m_histVector.begin(); iter!=m_histVector.end(); iter++ ) 
    (*iter)->Write(); 
  
  // probably not needed anymore, due to averagedHistograms function
  //   for ( vector< TH1D * >::const_iterator iter = m_histVectorVariableBin.begin(); iter!=m_histVectorVariableBin.end(); iter++ ) 
  //     (*iter)->Write(); 
  
  return true;
}

/**
return the last entry  of the histogram vector, needed for ATHENA (ATLAS software) implementation
*/
//TH1D* baseAnalysis::popHisto(std::vector<TH1D*> m_histVector) {
TH1D* baseAnalysis::popHisto() {

  // return NULL pointer if no histogram
  if (m_histVector.empty()) return 0;

  // get last element from vector (don't change vector)
  TH1D* temp = m_histVector.back();
  std::cout<<"m_histVector.back()"<< temp;
  // delete last entry from vector
  m_histVector.pop_back();
  return temp;
}

/// helper function: convert integer into string
const std::string intToString(const int i) {
  std::ostringstream s;
  s << i;
  return s.str();
}


/**
To get a better handling of the several Histogramms. This class initializes the histograms (name, title, binning, x-axis label) 
and collects them in the std::vector m_histVector. By doing this looping over all the 
histograms becomes much more convenient.
*/
TH1D* baseAnalysis::initHist(string name, string title, string xlabel, int nrBins, double xmin, double xmax)
{
  TDirectory* dir = gDirectory->GetDirectory(("/"+m_outputRootDir).c_str());
  if (dir) {
    dir->cd();
  } else {
    gDirectory->cd("/");
    dir = gDirectory->mkdir(m_outputRootDir.c_str());
    dir->cd();
  }

  // TH1D * histo=new TH1D((m_outputRootDir + "_" + intToString(m_histCounter[m_outputRootDir]) + "_" + name).c_str(),title.c_str(),nrBins,xmin,xmax);
  TH1D* histo=new TH1D((m_outputRootDir+"_"+intToString(m_histCounter[m_outputRootDir]) + "_" + name).c_str(),title.c_str(),nrBins,xmin,xmax);
  histo->GetXaxis()->SetTitle(xlabel.c_str());
  histo->GetYaxis()->SetTitle("Count");
  baseAnalysis::m_histVector.push_back(histo);

  ++m_histCounter[m_outputRootDir];

  return histo;
}

/**
histograms with variable bin size
*/
TH1D* baseAnalysis::initHistVariableBin(string name, string title, string xlabel, int nbin, Double_t nbinRange[])
{
  TDirectory* dir = gDirectory->GetDirectory(("/"+m_outputRootDir).c_str());
  if (dir) {
    dir->cd();
  } else {
    gDirectory->cd("/");
    dir = gDirectory->mkdir(m_outputRootDir.c_str());
    dir->cd();
  }

  TH1D* histoVariableBin=new TH1D(name.c_str(),title.c_str(),nbin,nbinRange);
  histoVariableBin->GetXaxis()->SetTitle(xlabel.c_str());
  histoVariableBin->GetYaxis()->SetTitle("Count");
  baseAnalysis::m_histVector.push_back(histoVariableBin);

  return histoVariableBin;
}

/**
This function returns true if the particle daughter is the decay product of the particle mother. 
The function scans the mother particle up to a certain level,
which can be specified with maxGenerations (negativ number --> all levels).
*/
bool baseAnalysis::checkDaughter(HepMC::GenParticle* mother, HepMC::GenParticle* daughter, int maxGenerations)
{
  bool found=false;
  int maxLoops=maxGenerations;
  int loop=1;
  
  // initialize the maximal number of mothers (Levels) to be checked
  if(maxGenerations<0) maxLoops=1000;
  
  // initialize current Production Vertex and mother particles
  // There should always be only one mother particle for decaying particles
  HepMC::GenVertex* current_vertex=daughter->production_vertex();
  HepMC::GenVertex::particle_iterator current_mother;
  
  // iterate through the mother particles and compare with mother
  // If found match, than stop the loop, otherwise try until there is no production vertex
  // or until you have reached the max number of loops specified
  while((current_vertex) && !(found) && (loop<maxLoops))
    {
      current_mother = current_vertex->particles_begin(HepMC::parents);
      
      if((*current_mother)==mother) 
        found=true;
      else
        current_vertex=(*current_mother)->production_vertex();
      
      loop++;
    }
  return found;
}


/** check if the track comes from a specific particle e.g. pid(W-Boson)=23 */
bool baseAnalysis::trackfromPID(int pid, HepMC::GenParticle* track, int maxGenerations) 
{
  bool found=false;
  int maxLoops=maxGenerations;
  int loop=1;
  
  // initialize the maximal number of mothers (Levels) to be checked
  if(maxGenerations<0) maxLoops=1000;
  
  // initialize current Production Vertex and mother particles
  // There should always be only one mother particle for decaying particles
  HepMC::GenVertex* current_vertex=track->production_vertex();
  HepMC::GenVertex::particle_iterator current_mother;
  
  // iterate through the mother particles and compare with mother
  // If found match, than stop the loop, otherwise try until there is no production vertex
  // or until you have reached the max number of loops specified
  while((current_vertex) && !(found) && (loop<maxLoops))
    {
      current_mother = current_vertex->particles_begin(HepMC::parents);
      
     //if no mother anymore, jump out the loop
      if(*current_mother==0)break;     

      if((*current_mother)->pdg_id()==pid) 
        found=true;
      else
        current_vertex=(*current_mother)->production_vertex();
      
      loop++;
    }
  return found;
}