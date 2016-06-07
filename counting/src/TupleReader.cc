#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <math.h>

#include "TLorentzVector.h"
#include "TParticle.h"
#include <memory>

#include "jpsiCounting/counting/interface/TupleReader.h"

TupleReader::TupleReader( std::vector<std::string> fileList ) {
   m_fileList = fileList;
   std::cout << "TupleReader: Initialized with a file list of length " << m_fileList.size() << "." << std::endl;
   m_file_it = m_fileList.begin();
   m_currentFile = 0;
   openNextFile();

   m_muonN_p4 = new TLorentzVector( 0., 0., 0., 0.);
   m_muonP_p4 = new TLorentzVector( 0., 0., 0., 0.);
   m_dimuon_p4 = new TLorentzVector( 0., 0., 0., 0.);
   return;
}
TupleReader::TupleReader(  ) {
   //~ return;
}
TupleReader::~TupleReader(){
   if( m_currentFile ) m_currentFile->Close();
}
bool TupleReader::openNextFile() {
   if( m_currentFile != 0 ) m_currentFile->Close();
   std::cout << "TupleReader: Reading file '" << *m_file_it << "'." << std::endl;
   m_currentFile = new TFile( (*m_file_it).c_str() );
   m_tree = (TTree*) m_currentFile->Get("rootuple/oniaTree");

   m_tree ->SetBranchAddress( "run",         &m_run );
   m_tree ->SetBranchAddress( "lumiblock",   &m_lumiblock );
   m_tree ->SetBranchAddress( "event",       &m_event );


   m_tree ->SetBranchAddress( "muonN_p4", &m_muonN_p4 );
   m_tree ->SetBranchAddress( "muonP_p4", &m_muonP_p4 );
   m_tree ->SetBranchAddress( "dimuon_p4", &m_dimuon_p4 );
   m_tree ->SetBranchAddress( "vProb", &m_dimuon_vProb );

   m_tree -> SetBranchAddress( "trigger", &m_trigger );
   m_tree -> SetBranchAddress( "numPrimaryVertices", &m_nvtx );
   m_eventsInTree = m_tree -> GetEntries();
   m_currentEvent = 0;

   m_file_it++;

   return true;
}

bool TupleReader::nextEvent() {
   if( m_currentEvent == ( m_eventsInTree - 1) ) {
      // We have reached the end of this file, so load next one
      if( m_file_it != m_fileList.end() ) {
         return openNextFile();
      } else {
         // No more files to load
         return false;
      }
   } else {
      // Load event from tree
      m_tree->GetEvent( m_currentEvent++ );
      return true;
   }

   return false;
}

bool TupleReader::setFileList(std::vector<std::string> & fileList) {
   m_fileList = fileList;
   return true;
}


