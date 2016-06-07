// system include files
#include <memory>
#include <math.h>
#include <bitset>

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TParticle.h"

#include "jpsiCounting/counting/interface/TupleReader.h"
#include "jpsiCounting/counting/interface/RunManager.h"



// Read Luminosity From File
// The file should contain two columns separated by white space
// First column:    run number
// Second columns:  luminosity for this run
std::map< unsigned int, double > readLumiFile( TString path ) {
    // Map to store luminosity information
    // Key:     Run number
    // Value:   Luminosity for the run
    std::map< unsigned int, double > lumimap;

    unsigned int run;
    double lumi;
    gSystem->ExpandPathName(path);
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
    //~ if( bits.test( 3 ) ) { tags.push_back( "_dimu10" ); }
    //~ if( bits.test( 6 ) ) { tags.push_back( "_dimu20" ); }
    //~ if( bits.any() ) { tags.push_back( "_all" ); };
    return tags;
}
TH1F* getVertexEfficiencyHisto() {
    TFile* f = TFile::Open("${CMSSW_BASE}/src/jpsiCounting/counting/data/TnP_Vertexing__data_all__vtx_ptpair.root", "r");
    TH1F* h_temp = (TH1F*)f->Get("tpTreeOnePair/Dimuon16_Jpsi_wrt_Dimuon6_Jpsi_NoVertexing");
    TH1F* h_weight = (TH1F*) h_temp->Clone();
    h_weight->SetDirectory(0);
    f->Close();
    return h_weight;
}
void makeStabilityPlot( std::vector<run> runs, std::string triggerTag, std::string normTag, TFile * file ) {
    std::cout << "Making the xs plot." << std::endl;

    // Initiate the histograms
    TH1D * candidatesPerRun = new TH1D( ("candidatesPerRun" + normTag + triggerTag ).c_str(),
                                         "J/Psi candidates detected in a given run;Run Number;J/Psi Count",
                                         runs.size(), 0, runs.size() );
    //~ TH1D * candidatesPerRun_corr = new TH1D( ("candidatesPerRun_corr" + normTag + triggerTag ).c_str(),
                                         //~ "J/Psi candidates detected in a given run, corrected for vertexing eff.;Run Number;J/Psi Count",
                                         //~ runs.size(), 0, runs.size() );
    TH1D * xsPerRun  = new TH1D( ("xsPerRun"  + normTag + triggerTag ).c_str(),
                                         "Cross section per run divided by average;Run Number;#sigma / #sigma_{average}",
                                         runs.size(), 0, runs.size() );
    //~ TH1D * xsPerRun_corr  = new TH1D( ("xsPerRun_corr"  + normTag + triggerTag ).c_str(),
                                         //~ "Stability of Luminosity (Normalized to mean of all runs), corrected for vertexing eff.;Run Number;J/Psi Count / Luminosity",
                                         //~ runs.size(), 0, runs.size() );
    TH1D * luminosityPerRun = new TH1D( ("luminosityPerRun" + normTag + triggerTag ).c_str(),
                                         "Luminosity in run via BrilCalc;Run Number; Luminosity",
                                         runs.size(), 0, runs.size() );
    TH1D * lumiSectionsPerRun = new TH1D( ("lumiSectionsPerRun" + normTag + triggerTag ).c_str(),
                                         "Luminosity sections run;Run Number; Number of Lumi Sections",
                                         runs.size(), 0, runs.size() );
    // Read the luminosity information for each run from file
    std::map< unsigned int, double > lumimap = readLumiFile( "$CMSSW_BASE/src/jpsiCounting/counting/data/lumiPerRun"+ normTag + triggerTag + ".txt" );

    // Read the efficiency information
    //~ TFile* f_vtx = TFile::Open("${CMSSW_BASE}/src/jpsiCounting/counting/data/TnP_Vertexing__data_all__vtx_ptpair.root", "r");
    //~ if( f_vtx == 0 or f_vtx->IsZombie() ) {
        //~ std::cout << "Vertexing file not found!" << std::endl;
        //~ return;
    //~ }
    //~ TH1F* h_vertexEfficiency = (TH1F*)f_vtx->Get("tpTreeOnePair/Dimuon16_Jpsi_wrt_Dimuon6_Jpsi_NoVertexing/fit_eff_plots/tag_nVertices_PLOT/");
    //~ if( h_vertexEfficiency == NULL ) {
        //~ std::cout << "Vertexing histogram not read!" << std::endl;
        //~ return;
    //~ }
    // Loop over runs
    int bin(1), nls(0);
    unsigned long runId(0);
    double lumi(0.);
    for( auto thisrun:runs ) {
        runId=thisrun.id;
        lumi=lumimap[runId];
        if( lumi == 0 ) {
            std::cout << "Found zero lumi in run: " << runId << std::endl;
            continue;
        }

        // Label the bin with the run number
        candidatesPerRun ->GetXaxis()->SetBinLabel( bin, std::to_string( runId ).c_str() );
        xsPerRun  ->GetXaxis()->SetBinLabel( bin, std::to_string( runId ).c_str() );
        luminosityPerRun ->GetXaxis()->SetBinLabel( bin, std::to_string( runId ).c_str() );

        // Add up the J/Psi counts for all Lumi Sections in this run
        int njpsi = thisrun.getN();
        //~ int njpsi_corr = thisrun.getN_weighted(h_vertexEfficiency);
        //~ int njpsi_corr_err = thisrun.getDN_weighted(h_vertexEfficiency);

        int njpsiFromLs=0;
        for( auto ls:thisrun.lumisections ) {
            njpsiFromLs += ls.getN();
        }
        if( njpsi != njpsiFromLs ) { std::cout << "Mismatch: " << njpsi << " vs " << njpsiFromLs << std::endl; }


        nls = thisrun.lumisections.size();
        candidatesPerRun->SetBinContent( bin, njpsi );
        candidatesPerRun->SetBinError(   bin, sqrt(njpsi) );

        //~ candidatesPerRun_corr->SetBinContent( bin, njpsi_corr );
        //~ candidatesPerRun_corr->SetBinError(   bin, njpsi_corr_err );

        xsPerRun->SetBinContent( bin, njpsi / lumi );
        xsPerRun->SetBinError(   bin, sqrt(njpsi) / lumi );

        //~ xsPerRun_corr->SetBinContent( bin, njpsi_corr / lumi );
        //~ xsPerRun_corr->SetBinError(   bin, njpsi_corr_err );

        luminosityPerRun->SetBinContent( bin, lumi );
        luminosityPerRun->SetBinError(   bin, 0 );

        lumiSectionsPerRun->SetBinContent( bin, nls );
        lumiSectionsPerRun->SetBinError(   bin, 0 );
        bin++;
    }

    // Normalize the xs plot to its average
    xsPerRun -> Scale( luminosityPerRun->Integral() / candidatesPerRun->Integral() );

    // Attach histograms to files
    candidatesPerRun                -> SetDirectory( file );
    xsPerRun                 -> SetDirectory( file );
    luminosityPerRun                -> SetDirectory( file );
    lumiSectionsPerRun              -> SetDirectory( file );

    //~ f_vtx->Close();

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

int main(int argc, char* argv[])
{
    std::cout << "runJpsiHistos: Preparations." << std::endl;
    // Initialize a bunch of stuff
    //~ std::vector<std::string> triggerTags = {"_dimu16","_dimu10","_dimu20", "_all"};
    std::vector<std::string> triggerTags = {"_dimu16"};
    std::vector<std::string> normTags = {"_moriond","_dtv2"};
    TFile * outputFile = TFile::Open("jpsi_histos.root", "RECREATE");
    std::map<std::string,TH1D*> histoMap = initHistoMap( triggerTags );


    // Input is a list of paths to files
    // Save them into a vector
    std::vector<std::string> inputFiles;
    for( int i = 1; i < argc; i++ ) {
        std::string mystring = *(argv+i);
        inputFiles.push_back( mystring );
    }

    RunManager manager;

    std::cout << "runJpsiHistos: Working on a total of '" << inputFiles.size() << "' files." << std::endl;

    // Initiate the reader and start event loop
    TupleReader reader( inputFiles );

    std::cout << "runJpsiHistos: Starting Event Loop." << std::endl;
    while( reader.nextEvent() ) {
        std::vector<std::string> triggerTags = getTriggerTags( reader.m_trigger );

        // Cut on mass
        if( reader.m_dimuon_p4->M() > 3.097+0.05 or reader.m_dimuon_p4->M() < 3.097-0.05 ) continue;
        // Cut on Muon Pt
        if( std::max(reader.m_muonP_p4->Pt(),reader.m_muonP_p4->Pt()) < 10 ) continue;
        if( std::min(reader.m_muonP_p4->Pt(),reader.m_muonP_p4->Pt()) < 10 ) continue;
        // Cut on Dimuon Pt
        if( ( reader.m_dimuon_p4->Pt() < 20 ) ) continue;
        if( abs( reader.m_dimuon_p4->Eta() ) > 1.5 ) continue;

        for( auto &triggerTag:triggerTags ) {
                histoMap["rec_pt_mu1"   + triggerTag] -> Fill( reader.m_muonN_p4->Pt() );
                histoMap["rec_pt_mu2"   + triggerTag] -> Fill( reader.m_muonP_p4->Pt() );
                histoMap["rec_pt_dimu"  + triggerTag] -> Fill( reader.m_dimuon_p4->Pt() );
                histoMap["rec_eta_mu1"  + triggerTag] -> Fill( reader.m_muonN_p4->Eta() );
                histoMap["rec_eta_mu2"  + triggerTag] -> Fill( reader.m_muonP_p4->Eta() );
                histoMap["rec_eta_dimu" + triggerTag] -> Fill( reader.m_dimuon_p4->Eta() );
                histoMap["rec_mass_dimu"+ triggerTag] -> Fill( reader.m_dimuon_p4->M() );
        }
        if( triggerTags.size() ) {
            manager.registerCandidate( reader.m_run, reader.m_lumiblock, reader.m_dimuon_p4, reader.m_nvtx );
        }

    }
    std::cout << "runJpsiHistos: Event loop finished." << std::endl;

    std::string normtag="_moriond";
    makeStabilityPlot( manager.getRuns(), "_dimu16",normtag, outputFile);

    manager.setOutputFile(outputFile);
    outputFile->Write();
    outputFile->Close();
    return 1;
}
