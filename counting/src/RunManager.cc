#include "jpsiCounting/counting/interface/RunManager.h"


std::vector<run> RunManager::getRuns() const{
   return runs;
}


bool RunManager::registerCandidate( unsigned long runId, unsigned long lumiId, TLorentzVector * p4, int nvertex, TupleReader const & reader ) {
   // Try to find the right run
   bool runFound = false;
   for( run & thisrun : runs ) {
      if( thisrun.id == runId ) {
         thisrun.registerCandidate( lumiId, p4, nvertex, reader );
         runFound = true;
         break;
      }
   }

   // If it does not exist yet, create it
   if( not runFound ) {
      run newrun( runId );
      newrun.registerCandidate( lumiId, p4, nvertex,reader );
      runs.push_back( newrun );
   }

   return 1;
}
bool RunManager::setOutputFile( TFile * outfile ) {
   for( auto &thisrun:runs ) {
      thisrun.h_mass -> SetDirectory( outfile );
      thisrun.h_pt -> SetDirectory( outfile );
      thisrun.h_eta -> SetDirectory( outfile );
      thisrun.h_nvtx -> SetDirectory( outfile );
      thisrun.h_vProb -> SetDirectory( outfile );
   }
   return true;
}




RunManager::RunManager() {
}
RunManager::~RunManager() {
}
