#include "TH1F.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <set>
#include <stdio.h>
#include <iostream>

float getValuefrom1DHist(double xval, TH1F* h_)
{
    if(h_==NULL) {
        std::cout << "[getValuefrom1DHist]: empty hist! " << std::endl;
        return 1;
    }
    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    float sf_ = h_->GetBinContent(binx);

    if(sf_==0.) return 1.;
    else return sf_;
}
float getErrorfrom1DHist(double xval, TH1F* h_)
{
    if(h_==NULL) {
        std::cout << "[getValuefrom1DHist]: empty hist! " << std::endl;
        return 1;
    }
    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    float sf_ = h_->GetBinError(binx);

    if(sf_==0.) return 1.;
    else return sf_;
}

void printEffTable(TString dir, TString meas, TString vars="pt:eta", TString method="fit") {
    RooDataSet *ds = (RooDataSet *) gFile->Get(dir+"/"+meas+"/"+method+"_eff");
    if (ds == 0) {
        std::cerr << "NOT FOUND: " << (dir+"/"+meas+"/"+method+"_eff") << std::endl;
        if (gFile->GetDirectory(dir) != 0) {
            gFile->GetDirectory(dir)->ls();
        } else {
            gFile->ls();
        }
    } else {
        ds->tree()->Scan("efficiency:efficiency_aerr_lo:efficiency_aerr_hi:"+vars);
    }
}
