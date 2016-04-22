#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <math.h>

#include "TLorentzVector.h"
#include "TParticle.h"
#include <memory>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
//~ #include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//~ #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//~ #include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
//~ #include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
//~ #ifndef TUPLEREADER_h
//~ #define TUPLEREADER_h

class TupleReader{
   public:
      TupleReader( std::vector<std::string> fileList );
      TupleReader( );
      ~TupleReader( );
      bool nextEvent();
      void readTreeFromFile();
      bool setFileList( std::vector<std::string> & fileList );
      bool openNextFile();
      TLorentzVector *m_muonN_p4, *m_muonP_p4, *m_dimuon_p4;
      ULong64_t m_event;
      UInt_t m_run, m_lumiblock, m_trigger;
      
   private:
      TBranch *m_branch_muonN_p4, *m_branch_muonP_p4, *m_branch_dimuon_p4;
      int m_eventsInTree, m_currentEvent;
      
      TFile * m_currentFile;
      std::vector<std::string> m_fileList;
      std::vector<std::string>::iterator m_file_it;
      TTree *m_tree;
      
};
//~ #endif
