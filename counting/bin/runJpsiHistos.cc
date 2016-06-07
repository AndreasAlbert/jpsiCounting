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
void makeStabilityPlot( std::vector<run> runs, std::string triggerTag, std::string normTag, TFile * file, bool doEfficiencies ) {
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
    /*
    //// Read the efficiency information
    // ID
    TFile * f_id = TFile::Open("/user/albert/lumi/CMSSW_7_6_3_patch2/src/MuonAnalysis/TagAndProbe/crab/TnP_MuonID__data_all__Soft2012_vtx.root","R");
    TCanvas * c_id  = (TCanvas*) f_id->Get("tpTree/Soft2012_vtx_Mu7p5_Track2_Jpsi/fit_eff_plots/tag_nVertices_PLOT_tag_Mu7p5_Track2_Jpsi_MU_pass");
    TGraphAsymmErrors* eff_id = (TGraphAsymmErrors*) c_id->GetPrimitive("hxy_fit_eff");

    // Tracking
    TFile * f_trk = TFile::Open("/user/albert/lumi/CMSSW_7_6_3_patch2/src/MuonAnalysis/TagAndProbe/crab/TnP_Tracking__data_all__vtx_trk.root","R");
    TCanvas * c_trk = (TCanvas*) f_trk->Get("tpTreeSta/eff_DR1p0_DEta0p4/fit_eff_plots/tag_nVertices_PLOT_Mu7p5_L2Mu2_Jpsi_L2_pass_&_tag_Mu7p5_L2Mu2_Jpsi_MU_pass");
    TGraphAsymmErrors* eff_trk = (TGraphAsymmErrors*) c_trk->GetPrimitive("hxy_fit_eff");

    // Vertexing
    TFile * f_vtx = TFile::Open("/user/albert/lumi/CMSSW_7_6_3_patch2/src/MuonAnalysis/TagAndProbe/crab/TnP_Vertexing__data_all__vtx_ptpair.root","R");
    TCanvas * c_vtx = (TCanvas*) f_vtx->Get("tpTreeOnePair/Dimuon16_Jpsi_wrt_Dimuon6_Jpsi_NoVertexing/fit_eff_plots/tag_nVertices_PLOT_Dimuon6_Jpsi_NoVertexing_pass_&_tag_Dimuon6_Jpsi_NoVertexing_pass");
    TGraphAsymmErrors* eff_vtx = (TGraphAsymmErrors*) c_vtx->GetPrimitive("hxy_fit_eff");
    */

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
        TH1F* h = thisrun.h_nvtx;
        double eff1(1.), eff2(1.), eff3(1.);
        double deff1(0.), deff2(0.), deff3(0.);
        double njpsi(0.),dnjpsisqr(0.);

        for( int i = 1; i < h -> GetNbinsX(); i++ ){
            /*
            getValuefromTGraph( h->GetBinCenter(i), eff_id, eff1, deff1 );
            getValuefromTGraph( h->GetBinCenter(i), eff_trk, eff2, deff2 );
            getValuefromTGraph( h->GetBinCenter(i), eff_vtx, eff3, deff3 );
            */
            double this_njpsi(0.), this_dnjpsisqr(0.);
            if( doEfficiencies ) {
                this_njpsi = h->GetBinContent(i) / (eff1 * eff2 * eff3 );
                this_dnjpsisqr = pow(this_njpsi,2) * (1 / this_njpsi + pow( deff1 / eff1, 2 )+ pow( deff2 / eff2, 2 )+ pow( deff3 / eff3, 2 ));
            } else {
                this_njpsi = h->GetBinContent(i);
                this_dnjpsisqr = this_njpsi;
            }
            //~ std::cout << this_njpsi << " " << sqrt(this_dnjpsisqr) << std::endl;
            if( this_njpsi > 0 ) {
                njpsi += this_njpsi;
                dnjpsisqr += this_dnjpsisqr;
            }

        }

        nls = thisrun.lumisections.size();
        candidatesPerRun->SetBinContent( bin, njpsi );
        candidatesPerRun->SetBinError(   bin, sqrt(dnjpsisqr) );

        xsPerRun->SetBinContent( bin, njpsi / lumi );
        xsPerRun->SetBinError(   bin, sqrt(dnjpsisqr) / lumi );

        luminosityPerRun->SetBinContent( bin, lumi );
        luminosityPerRun->SetBinError(   bin, 0 );

        lumiSectionsPerRun->SetBinContent( bin, nls );
        lumiSectionsPerRun->SetBinError(   bin, 0 );
        bin++;
    }

    // Normalize the xs plot to its average
    std::vector<float> cross_sections;
    for(int i=1; i<luminosityPerRun->GetNbinsX(); i++) {
        cross_sections.push_back( candidatesPerRun->GetBinContent (i) / luminosityPerRun->GetBinContent(i) );
    }
    std::nth_element(cross_sections.begin(), cross_sections.begin() + cross_sections.size()/2, cross_sections.end());
    xsPerRun -> Scale( 1 / cross_sections[cross_sections.size()/2] );

    // Attach histograms to files
    candidatesPerRun    -> SetDirectory( file );
    xsPerRun            -> SetDirectory( file );
    luminosityPerRun    -> SetDirectory( file );
    lumiSectionsPerRun  -> SetDirectory( file );
/*
    f_id->Close();
    f_vtx->Close();
    f_trk->Close();
*/
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
        // Check for Triggers
        std::vector<std::string> triggerTags = getTriggerTags( reader.m_trigger );
        if( not triggerTags.size() ) continue;


        // Order muons by pT
        TLorentzVector *mu1, *mu2;
        if ( reader.m_muonP_p4->Pt() > reader.m_muonN_p4->Pt() ) {
            mu1 = reader.m_muonP_p4;
            mu2 = reader.m_muonN_p4;
        } else {
            mu1 = reader.m_muonN_p4;
            mu2 = reader.m_muonP_p4;
        }

        // Cut on mass
        if( reader.m_dimuon_p4->M() > (3.097+0.05) or reader.m_dimuon_p4->M() < (3.097-0.05) ) continue;
        // Cut on Muon Pt
        if( mu1->Pt() < 11 ) continue;
        if( mu2->Pt() < 4.5 ) continue;

        // Cut on Muon Eta
        if( fabs(mu1->Eta()) > 2.1 ) continue;
        if( fabs(mu2->Eta()) > 2.1 ) continue;

        // Cut on Dimuon Pt
        //~ if( ( reader.m_dimuon_p4->Pt() < 20 ) ) continue;
        //~ if( ( reader.m_dimuon_vProb > 0.05 ) ) continue;
        //~ if( fabs( reader.m_dimuon_p4->Eta() ) > 1.5 ) continue;

        if ( reader.m_dimuon_vProb < 0.01 ) continue;
        for( auto &triggerTag:triggerTags ) {
            histoMap["rec_pt_mu1"   + triggerTag] -> Fill( mu1->Pt() );
            histoMap["rec_pt_mu2"   + triggerTag] -> Fill( mu2->Pt() );
            histoMap["rec_pt_dimu"  + triggerTag] -> Fill( reader.m_dimuon_p4->Pt() );

            histoMap["rec_eta_mu1"  + triggerTag] -> Fill( mu1->Eta() );
            histoMap["rec_eta_mu2"  + triggerTag] -> Fill( mu2->Eta() );
            histoMap["rec_eta_dimu" + triggerTag] -> Fill( reader.m_dimuon_p4->Eta() );

            histoMap["rec_mass_dimu"+ triggerTag] -> Fill( reader.m_dimuon_p4->M() );
        }
        manager.registerCandidate( reader.m_run, reader.m_lumiblock, reader.m_dimuon_p4, reader.m_nvtx, reader );
    }
    std::cout << "runJpsiHistos: Event loop finished." << std::endl;

    std::string normtag="_moriond";
//~
    //~ TFile * outputFile_noEff = TFile::Open("jpsi_histos_noEff.root", "RECREATE");
    //~ makeStabilityPlot( manager.getRuns(), "_dimu16",normtag, outputFile_noEff, false);
    //~ outputFile_noEff->Write();
    //~ outputFile_noEff->Close();

    makeStabilityPlot( manager.getRuns(), "_dimu10",normtag, outputFile, false);




    manager.setOutputFile(outputFile);
    outputFile->Write();
    outputFile->Close();

    return 1;
}
