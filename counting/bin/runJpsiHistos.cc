// system include files
#include <memory>

// user include files
//~ #include "FWCore/Framework/interface/Frameworkfwd.h"
//~ #include "FWCore/Framework/interface/EDAnalyzer.h"
//~ #include "FWCore/Framework/interface/Event.h"
//~ #include "FWCore/Framework/interface/MakerMacros.h"
//~ #include "FWCore/ParameterSet/interface/ParameterSet.h"
//~ #include "FWCore/ServiceRegistry/interface/Service.h"
//~
//~ #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//~ #include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
//~ #include "DataFormats/PatCandidates/interface/Muon.h"
//~ #include "DataFormats/Candidate/interface/Candidate.h"
//~ #include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//~ #include "DataFormats/VertexReco/interface/VertexFwd.h"
//~
//~ #include "DataFormats/Common/interface/TriggerResults.h"
//~ #include "FWCore/Common/interface/TriggerNames.h"
//~
//~ #include "CommonTools/UtilAlgos/interface/TFileService.h"
//~ #include "TLorentzVector.h"
//~ #include "TTree.h"



#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <math.h>

#include "TLorentzVector.h"
#include "TParticle.h"
#include <memory>
#include "jpsiCounting/counting/interface/TupleReader.h"
#include <bitset>
struct runs_t {
    // Outer Map Key:   Run Identifier
    // Inner Map Key:   Lumi Section Identifier
    // Inner Value:     Number of J/Psi candidates
    // example: cand[runIdentifier][lumiSectionIdentifier] = numberOfCandidates;
    std::map< unsigned int, std::map< unsigned int, unsigned int > >cand;

    // Check if a certain run has been registered
    bool hasRun( unsigned int run ) { return cand.count( run ); };

    // Check if a certain run and lumi section have been registered
    bool hasRunAndLumi( unsigned int run, unsigned int ls ) { return cand.count( run ) && cand[run].count( ls ); }

    // Register a J/Psi candidate with given run and lumi section
    // If no candidates have been registered before, the counter is initialized.
    // If candidates have been registered before, the counter is incremented.
    void registerCandidate( unsigned int run, unsigned int ls ) {
        if( hasRunAndLumi( run, ls ) ) {
            cand[run][ls]++;
        } else if ( hasRun( run ) ) {
            cand[run][ls] = 1;
        } else {
            cand[run] = std::map< unsigned int, unsigned int >();
        }
    }

    int getTotalCount() {
        int njpsi(0);
        for( auto &perRun:cand ) {
            for( auto &perLumiSection:perRun.second ) {
                njpsi+=perLumiSection.second;
            }
        }
        return njpsi;
    }
};

// Read Luminosity From File
// The file should contain two columns separated by white space
// First column:    run number
// Second columns:  luminosity for this run
std::map< unsigned int, double > readLumiFile( std::string path ) {
    // Map to store luminosity information
    // Key:     Run number
    // Value:   Luminosity for the run
    std::map< unsigned int, double > lumimap;

    unsigned int run;
    double lumi;
    std::ifstream f( path );

    if( not f ) {
        std::cout << "readLumiFile: File not found: " << path << std::endl;
        throw;
    }
    while(!f.eof()) {
        f >> run >> lumi;
        lumimap[run] = lumi;
    }
    f.close();
    return lumimap;
}

// Decode trigger information
// The convention corresponds to that in the Muonia Tupler
std::vector<std::string> getTriggerTags( int trigger ) {
    std::bitset<16> bits( trigger );
    std::vector<std::string> tags;
    if( bits.test( 0 ) ) { tags.push_back( "_dimu16" ); }
    if( bits.test( 3 ) ) { tags.push_back( "_dimu10" ); }
    if( bits.test( 6 ) ) { tags.push_back( "_dimu20" ); }
    if( bits.any() ) { tags.push_back( "_all" ); };
    return tags;
}

void makeStabilityPlot( runs_t runs, std::string triggerTag, std::string normTag, TFile * file ) {
    std::cout << "Making the stability plot." << std::endl;

    // Initiate the histograms
    TH1D * candidatesPerRun = new TH1D( ("candidatesPerRun" + normTag + triggerTag ).c_str(),
                                         "J/Psi candidates detected in a given run;Run Number;J/Psi Count",
                                         runs.cand.size(), 0, runs.cand.size() );
    TH1D * stabilityPerRun  = new TH1D( ("stabilityPerRun"  + normTag + triggerTag ).c_str(),
                                         "Stability of Luminosity (Normalized to mean of all runs);Run Number;J/Psi Count / Luminosity",
                                         runs.cand.size(), 0, runs.cand.size() );
    TH1D * luminosityPerRun = new TH1D( ("luminosityPerRun" + normTag + triggerTag ).c_str(),
                                         "Luminosity in run via BrilCalc;Run Number; Luminosity",
                                         runs.cand.size(), 0, runs.cand.size() );
    TH1D * lumiSectionsPerRun = new TH1D( ("lumiSectionsPerRun" + normTag + triggerTag ).c_str(),
                                         "Luminosity sections run;Run Number; Number of Lumi Sections",
                                         runs.cand.size(), 0, runs.cand.size() );
    // Read the luminosity information for each run from file
    std::map< unsigned int, double > lumimap = readLumiFile( "norm/lumiPerRun"+ normTag + triggerTag + ".txt" );


    // Loop over runs
    int bin(1), nls(0);
    for( auto &lsmap:runs.cand ) {
        // Label the bin with the run number
        candidatesPerRun ->GetXaxis()->SetBinLabel( bin, std::to_string( lsmap.first ).c_str() );
        stabilityPerRun  ->GetXaxis()->SetBinLabel( bin, std::to_string( lsmap.first ).c_str() );
        luminosityPerRun ->GetXaxis()->SetBinLabel( bin, std::to_string( lsmap.first ).c_str() );

        // Add up the J/Psi counts for all Lumi Sections in this run
        int njpsi(0);
        for( auto &ls:lsmap.second ) {
            njpsi += ls.second;
            nls++;
        }
        candidatesPerRun->SetBinContent( bin, njpsi );
        candidatesPerRun->SetBinError( bin, sqrt(njpsi) );

        stabilityPerRun->SetBinContent( bin, njpsi / lumimap[lsmap.first] );
        stabilityPerRun->SetBinError( bin, sqrt(njpsi) / lumimap[lsmap.first] );

        luminosityPerRun->SetBinContent( bin, lumimap[lsmap.first] );
        luminosityPerRun->SetBinError( bin, 0 );

        lumiSectionsPerRun->SetBinContent( bin, nls );
        lumiSectionsPerRun->SetBinError( bin, 0 );
        bin++;
    }

    // Normalize the stability plot to its average
    stabilityPerRun -> Scale( luminosityPerRun->Integral() / candidatesPerRun->Integral() );

    // Attach histograms to files
    candidatesPerRun                -> SetDirectory( file );
    stabilityPerRun                 -> SetDirectory( file );
    luminosityPerRun                -> SetDirectory( file );
    lumiSectionsPerRun              -> SetDirectory( file );
}
std::map<std::string, TH1D*> initHistoMap(std::vector<std::string> triggerTags) {
    std::map<std::string,TH1D*> histoMap;

    for( auto &triggerTag:triggerTags ) {
        // Single Mu
        histoMap["rec_eta_mu1"  +triggerTag] = new TH1D( ("rec_eta_mu1" +triggerTag).c_str(), "Rec. #eta of leading muon", 50, -2.5, 2.5 );
        histoMap["rec_eta_mu2"  +triggerTag] = new TH1D( ("rec_eta_mu2" +triggerTag).c_str(), "Rec. #eta of trailing muon", 50, -2.5, 2.5 );
        histoMap["rec_pt_mu1"   +triggerTag] = new TH1D( ("rec_pt_mu1"  +triggerTag).c_str(), "Rec. p_{T} of leading muon", 100, 0., 100. );
        histoMap["rec_pt_mu2"   +triggerTag] = new TH1D( ("rec_pt_mu2"  +triggerTag).c_str(), "Rec. p_{T} of trailing muon", 100, 0., 100. );
        // DiMu
        histoMap["rec_eta_dimu" +triggerTag] = new TH1D( ("rec_eta_dimu"  +triggerTag).c_str(), "Rec. #eta of the di-muon system", 50, -2.5, 2.5 );
        histoMap["rec_pt_dimu"  +triggerTag] = new TH1D( ("rec_pt_dimu"   +triggerTag).c_str(), "Rec. p_{T} of the di-muon system", 100, 0., 100. );
        histoMap["rec_mass_dimu"+triggerTag] = new TH1D( ("rec_mass_dimu" +triggerTag).c_str(), "Rec. mass of the di-muon system", 1000, 0., 10. );
    }

    return histoMap;
}
std::map<std::string, runs_t> initRunsMap(std::vector<std::string> triggerTags ) {
    std::map<std::string,runs_t> runsMap;
    for( auto &triggerTag:triggerTags ) {
            runsMap[ triggerTag ] = runs_t();
    }
    return runsMap;
}

int main(int argc, char* argv[])
{
    std::cout << "runJpsiHistos: Preparations." << std::endl;
    // Initialize a bunch of stuff
    std::vector<std::string> triggerTags = {"_dimu16","_dimu10","_dimu20", "_all"};
    std::vector<std::string> normTags = {"_moriond","_dtv2"};
    TFile * outputFile = TFile::Open("jpsi_histos.root", "RECREATE");
    std::map<std::string,TH1D*> histoMap = initHistoMap( triggerTags );
    std::map<std::string,runs_t> runsMap = initRunsMap( triggerTags );



    // Input is a list of paths to files
    // Save them into a vector
    std::vector<std::string> inputFiles;
    for( int i = 1; i < argc; i++ ) {
        std::string mystring = *(argv+i);
        inputFiles.push_back( mystring );
    }


    std::cout << "runJpsiHistos: Working on a total of '" << inputFiles.size() << "' files." << std::endl;

    // Initiate the reader and start event loop
    TupleReader reader( inputFiles );

    std::cout << "runJpsiHistos: Starting Event Loop." << std::endl;
    while( reader.nextEvent() ) {
        std::vector<std::string> triggerTags = getTriggerTags( reader.m_trigger );

        // Cut on mass
        if( reader.m_dimuon_p4->M() > 3.3 or reader.m_dimuon_p4->M() < 2.9 ) continue;

        for( auto &triggerTag:triggerTags ) {
                histoMap["rec_pt_mu1"   + triggerTag] -> Fill( reader.m_muonN_p4->Pt() );
                histoMap["rec_pt_mu2"   + triggerTag] -> Fill( reader.m_muonP_p4->Pt() );
                histoMap["rec_pt_dimu"  + triggerTag] -> Fill( reader.m_dimuon_p4->Pt() );
                histoMap["rec_eta_mu1"  + triggerTag] -> Fill( reader.m_muonN_p4->Eta() );
                histoMap["rec_eta_mu2"  + triggerTag] -> Fill( reader.m_muonP_p4->Eta() );
                histoMap["rec_eta_dimu" + triggerTag] -> Fill( reader.m_dimuon_p4->Eta() );
                histoMap["rec_mass_dimu"+ triggerTag] -> Fill( reader.m_dimuon_p4->M() );

                runsMap[triggerTag].registerCandidate( reader.m_run, reader.m_lumiblock );
        }
    }
    std::cout << "runJpsiHistos: Event loop finished." << std::endl;

    std::cout << "runJpsiHistos: Make JPsi Count/Lumi/Stability Histograms." << std::endl;
    for( auto &normTag : normTags ){
        for( auto &runs_it : runsMap ) {
            makeStabilityPlot( runs_it.second, runs_it.first, normTag, outputFile );
        }
    }

    std::cout << "runJpsiHistos: Finishing up." << std::endl;
    outputFile->Write();
    outputFile->Close();

    // Print the number of runs, lumi sections and J/Psi candidates we have registered
    int nrun(0), nls(0), njpsi(0);
    for( auto &lsmap:runsMap["_all"].cand ) {
        nrun++;
        for( auto &ls:lsmap.second ) {
            nls++;
            njpsi += ls.second;
        }
    }


    std::cout << "runJpsiHistos: Observed " << nrun << " runs, " << nls << " lumi sections and " << njpsi << " J/Psi candidates." << std::endl;
    std::cout << "runJpsiHistos: Exiting." << std::endl;

    return 1;
}
