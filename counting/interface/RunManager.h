#include "TH1F.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <set>
#include <stdio.h>
#include <iostream>
#include "Util.h"
#include "math.h"
#include "jpsiCounting/counting/interface/TupleReader.h"

struct timebinning;
struct run;
struct lumisection;

//// Template struct defining members common to runs and lumisections
struct timebinning {
   unsigned long id;
   TH1F *h_mass, *h_pt, *h_eta, *h_nvtx,*h_vProb, *h_effweight;
   timebinning ( unsigned long inputId ) {
      id = inputId;
      h_mass = new TH1F( ( std::to_string(inputId)+"_mass" ).c_str(),( std::to_string(inputId)+"_mass" ).c_str(),100,2.9,3.3);
      h_pt = new TH1F( ( std::to_string(inputId)+"_pt" ).c_str(),( std::to_string(inputId)+"_pt" ).c_str(),100,0.,50.);
      h_eta = new TH1F( ( std::to_string(inputId)+"_eta" ).c_str(),( std::to_string(inputId)+"_eta" ).c_str(),100,-2.5,2.5);
      h_nvtx = new TH1F(  ( std::to_string(inputId)+"_nvtx" ).c_str(),( std::to_string(inputId)+"_nvtx" ).c_str(),31,-0.5,30.5);
      h_vProb = new TH1F(  ( std::to_string(inputId)+"_vProb" ).c_str(),( std::to_string(inputId)+"_vProb" ).c_str(),100,-0.1,0.9);
      h_mass->SetDirectory(0);
      h_pt->SetDirectory(0);
      h_eta->SetDirectory(0);
      h_nvtx->SetDirectory(0);
      h_vProb->SetDirectory(0);
   }
   bool operator<(const timebinning other) const{
      return id<other.id;
   }
   int getN() {
      return h_mass->GetEntries();
   }
   
   double getN_weighted(TH1F* h_weight) {
      double N = 0;
      for( int i=1; i <= h_nvtx->GetNbinsX(); i++ ) {
         int njpsi = h_nvtx->GetBinContent(i);
         int nvtx =  h_nvtx->GetBinCenter(i);
         double eff = getValuefrom1DHist(nvtx,h_weight);
         N += njpsi / eff;
      }
      return N;
   }
   double getDN_weighted(TH1F* h_weight) {
      double DNSQR = 0;
      for( int i=1; i <= h_nvtx->GetNbinsX(); i++ ) {
         int nvtx =  h_nvtx->GetBinCenter( i );                                 // Number of vertices for this bin
         double eff = getValuefrom1DHist(nvtx, h_weight);                              // Get the efficiency, function of nvtx
         double eff_err = getValuefrom1DHist(nvtx, h_weight);                           // Get the efficiency uncertainty
         DNSQR += pow( h_nvtx->GetBinError(i) / eff, 2 );                       // Error propagation
         DNSQR += pow( h_nvtx->GetBinContent(i) * eff_err / pow( eff, 2 ), 2 ); // Error propagation
      }
      return sqrt( DNSQR );
   }

   
   bool fillHistos( TLorentzVector * p4, int nvertex,TupleReader const & reader ) {
      h_mass->Fill( p4->M() );
      h_pt->Fill( p4->Pt() );
      h_eta->Fill( p4->Eta() );
      h_nvtx->Fill(nvertex);
      h_vProb->Fill(reader.m_dimuon_vProb);
      return true;
   }
};


//// Struct with lumi section information
struct lumisection:timebinning {
   lumisection(unsigned long inputId):timebinning(inputId){};
   bool registerCandidate( TLorentzVector * p4, int nvertex, TupleReader const & reader ) {
      fillHistos(p4, nvertex, reader);
      return 1;
   }
   
};
//// Struct with run information
struct run:timebinning {
   run(unsigned long inputId):timebinning(inputId){};
   std::vector<lumisection> lumisections;

   bool registerCandidate( unsigned long lumiId, TLorentzVector * p4, int nvertex,TupleReader const & reader ) {

      // Try to Find the right lumi section
      bool lsFound = false;
      for( lumisection & thisls : lumisections ) {
         if( thisls.id == lumiId ) {
            thisls.registerCandidate( p4, nvertex, reader );
            lsFound = true;
            break;
         }
      }
      // If it does not exist yet, create it
      if( not lsFound ) {
         lumisection newls( lumiId );
         newls.registerCandidate( p4, nvertex, reader );
         lumisections.push_back( newls );
      }

      fillHistos( p4, nvertex, reader );

      return 1;
   }


};

class RunManager {
   public:
      RunManager();
      ~RunManager();
      std::vector<run> getRuns() const;
      bool registerCandidate( unsigned long runId, unsigned long lumiId, TLorentzVector * p4, int nvertex, TupleReader const & reader );
      bool setOutputFile( TFile * outfile );
   private:
      std::vector<run> runs;
};
