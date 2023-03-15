#ifndef baseAnalysis_H
#define baseAnalysis_H

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "HepMC/GenEvent.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

// forward declarations
namespace fastjet {
  class ClusterSequence;
}
namespace HepPDT {
  class ParticleDataTable;
}
class TFile;
class TH1D;

//sj//enum RETURN{FAILURE=-1, SUCCESS=1, CHARGED=2, NOTCHARGED=3, FINALSTATEPARTICLE=4};

/**
@class baseAnalysis.h
@brief Base class for all  HepMCanalysis classes to implement common 
       functionality and interfaces.

@author Cano Ay Dec 2008 */

class baseAnalysis
{
 public:

  baseAnalysis();
  virtual ~baseAnalysis();
  
  virtual int Finalize();
  virtual int Finalize(TFile* output);
  //virtual TH1D* popHisto(std::vector<TH1D*>);
  virtual TH1D* popHisto();

  virtual int Init(double maxeta=2.5, double minpt=0.5) { std::cout << "baseAnalysis: WARNING: You have called the dummy function Init()." << std::endl; return 0;};
  virtual int Process(const HepMC::GenEvent* hepmcevt) { std::cout << "baseAnalysis: WARNING: You have called the dummy function Process()." << std::endl; return 0;} ;
  //virtual std::vector<TH1D*> averagedHistograms() {TH1D* temp; std::vector<TH1D*> temp_vec; temp_vec.clear(); temp_vec.push_back(temp); return temp_vec;} ;
  virtual std::vector<TH1D*> averagedHistograms() {return m_histVector;} ;

  /** Initialize FastJet Algorithm*/
  int InitJetFinder(double coneRadius, double overlapThreshold, double jet_ptmin, double lepton_ptmin, double DeltaR_lepton_track);
  int FindJetableParticles(const HepMC::GenEvent* hepmcevt);
  int FindJet(const HepMC::GenEvent* hepmcevt);
  int DeleteJetObject(const HepMC::GenEvent* hepmcevt);
  std::vector<fastjet::PseudoJet> GetJet(const HepMC::GenEvent* hepmcevt);
  virtual void SetJet(std::vector<fastjet::PseudoJet>* inclusive_jets) {m_inclusive_jets=*inclusive_jets;};

  /** Calculate Missing Et*/
  int FindMissingEt(const HepMC::GenEvent* hepmcevt);
  int ClearMissingEt(const HepMC::GenEvent* hepmcevt);

  /** Clear some values from the event */
  int ClearEvent(const HepMC::GenEvent* hepmcevt);

  /** Set the Output filename*/  
  inline int setOutpuFileName(const char* filename){ m_outputFileName=filename; return 0;};
  
  /** Set the directory name in Output root file*/  
  inline int setOutpuRootDir(const char* dirname){ m_outputRootDir=dirname; return 0;};

  /** Check some special neutral Particles */
  inline bool IsNeutrino(int pid) { if( abs(pid)==12 || abs(pid)==14 || abs(pid)==16) return true; else return false; };
  inline bool IsGamma(int pid) { if( abs(pid)==22 ) return true; else return false; };
  inline bool IsNeutron(int pid) { if( abs(pid)==2112 ) return true; else return false; };
  inline bool IsK0(int pid) { if( abs(pid)==130 || abs(pid)==310 || abs(pid)==311 ) return true; else return false; };
  inline bool IsPi0(int pid) { if( abs(pid)==111 || abs(pid)==113 || abs(pid)==221 ) return true; else return false; };
  inline bool IsElectron(int pid) { if( abs(pid)==11 ) return true; else return false; };
  inline bool IsMuon(int pid) { if( abs(pid)==13 ) return true; else return false; };
  inline bool IsGluon(int pid) { if( abs(pid)==21 ) return true; else return false; };


  /** Check if neutral particle*/  
  inline int chargedParticle(int pid) { 
    if( IsNeutrino(pid) || IsGamma(pid) || IsNeutron(pid) || IsK0(pid) || IsPi0(pid) || IsGluon(pid) || abs(pid)== 94 || abs(pid)==1000039 || abs(pid) ==1000022 || abs(pid) ==39 || abs(pid) ==5100022) //above neutral particles and some "special" neutral particles
      {return false;} 
    else 
      {return true;}
    return 0;};
  
  /** Check if final state particle*/
  virtual int IsFinalStateParticle(HepMC::GenParticle *p);

  /** Set the maximum allowed eta range*/  
  inline int setMaxEtaCut(double maxeta){ m_max_eta=maxeta; return 0;};
  
  /** Set maximum pt of tracks */  
  inline int setMinPtCut(double minpt){ m_min_pt=minpt; return 0;};
  
  /** Initialization of histograms  */ 
  TH1D* initHist(std::string name, std::string title, std::string xlabel, int nrBins=100, double xmin=0., double xmax=100.) ;
  TH1D* initHistVariableBin(std::string name, std::string title, std::string xlabel, int nbin, double nbinRange[]) ;
  
  /** check if mother decayed into daugther */ 
  bool checkDaughter(HepMC::GenParticle* mother, HepMC::GenParticle* daughter, int maxGenerations=-1) ;
  
  /** check if the track comes from a specific particle e.g. pid(W-Boson)=24 */ 
  bool trackfromPID(int pid,HepMC::GenParticle* track, int maxGenerations=-1) ;

  /** calculate the rapidity of a particle */
  inline double getRapidity(HepMC::GenEvent::particle_const_iterator p) {
    double e = (*p)->momentum().e();
    double pz = (*p)->momentum().pz();
    double rapidity = 0.5 * log((e + pz) / (e - pz));
    return rapidity;
  };

 protected:
  
  /** jet finding */
  fastjet::JetDefinition::Plugin*   m_plugin;
  fastjet::JetDefinition*           m_jetDef;

  std::vector<fastjet::PseudoJet>  m_input_particles;
  std::vector<fastjet::PseudoJet> m_inclusive_jets;
  fastjet::ClusterSequence* m_clust_seq;

  std::map< std::string, int > m_histCounter;
  std::vector<TH1D*> m_histVector;
  std::vector<TH1D*> m_histVectorVariableBin;  
  std::string m_outputFileName;
  std::string m_outputRootDir;

  //for jet finder
  double m_coneRadius;
  double m_overlapThreshold;
  double m_jet_ptmin;
  double m_lepton_ptmin;
  double m_DeltaR_lepton_track;
  bool m_Jetfinder_enabled;
  
  // specify the maximum eta of tracks
  double m_max_eta;
  double m_min_pt;

  // for calculation of the missing energy
  double exMissTruth;
  double eyMissTruth;
  double etMissTruth;
  double etsumMissTruth;
  
  const HepPDT::ParticleDataTable* m_particleTable;
};

#endif