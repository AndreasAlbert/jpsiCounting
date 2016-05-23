#include "jpsiCounting/counting/interface/RunManager.h"


std::vector<run> RunManager::getRuns() const{
   return runs;
}


bool RunManager::registerCandidate( unsigned long runId, unsigned long lumiId, TLorentzVector * p4, int nvertex ) {
   // Try to find the right run
   bool runFound = false;
   for( run & thisrun : runs ) {
      if( thisrun.id == runId ) {
         thisrun.registerCandidate( lumiId, p4, nvertex );
         runFound = true;
         break;
      }
   }

   // If it does not exist yet, create it
   if( not runFound ) {
      run newrun( runId );
      newrun.registerCandidate( lumiId, p4, nvertex );
      runs.push_back( newrun );
   }

   return 1;
}
bool RunManager::setOutputFile( TFile * outfile ) {
   for( auto &thisrun:runs ) {
      thisrun.h_mass -> SetDirectory( outfile );
      thisrun.h_nvtx -> SetDirectory( outfile );
   }
   return true;
}




RunManager::RunManager() {
}
RunManager::~RunManager() {
}
