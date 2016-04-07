#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>

double GetEfficiency(TTree* tree);

int main(int argc, char *argv[]) {
    
    // argc should be =2 for correct execution
    if ( argc != 2 ) {
        // We print argv[0] assuming it is the program name
        std::cout << "usage: " << argv[0] << " <efficiency.root>" << std::endl;
        return 0;
    }
    
    TString fileName = argv[1];
    
    // open file
    TFile* File = TFile::Open(fileName);
    
    if (!File) {
        std::cout << "##### ERROR: could not open " << fileName << std::endl;
        return 0;
    }

    // get information on run
    TTree* RunTree = (TTree*) File->Get("RunTree");
    
    double m_WIMP, delta, delta_low, delta_high, mu_WIMP, mu_WIMP_low, mu_WIMP_high;
    int nSteps_delta, nSteps_mu;
    
    RunTree->SetBranchAddress("m_WIMP",&m_WIMP);
    RunTree->SetBranchAddress("nSteps_delta",&nSteps_delta);
    RunTree->SetBranchAddress("nSteps_mu",&nSteps_mu);
    RunTree->SetBranchAddress("delta_low",&delta_low);
    RunTree->SetBranchAddress("delta_high",&delta_high);
    RunTree->SetBranchAddress("mu_WIMP_low",&mu_WIMP_low);
    RunTree->SetBranchAddress("mu_WIMP_high",&mu_WIMP_high);
    
    RunTree->GetEntry(0);
    
    // calculate histogram bins
    double bin_delta_low = delta_low-(delta_high-delta_low)/(nSteps_delta-1)/2.;
    double bin_delta_high = delta_high+(delta_high-delta_low)/(nSteps_delta-1)/2.;
    double bin_mu_low = mu_WIMP_low-(mu_WIMP_high-mu_WIMP_low)/(nSteps_mu-1)/2.;
    double bin_mu_high = mu_WIMP_high+(mu_WIMP_high-mu_WIMP_low)/(nSteps_mu-1)/2.;
    
    TH2D* hist = new TH2D("hist",";#delta (keV);#mu_{#chi} (#times10^{-3}#mu_{N})",nSteps_delta,bin_delta_low,bin_delta_high,nSteps_mu,bin_mu_low,bin_mu_high);
    double efficiency;
    
    int Ntot = nSteps_delta*nSteps_mu;
    int N = 0;
    
    // loop over all trees
    for (int i=0; i<nSteps_delta; ++i) {
        
        delta=delta_low+i*(delta_high-delta_low)/(nSteps_delta-1);
        
        for (int j=0; j<nSteps_mu; ++j) {
            N++;
            std::cout << "#### analyzing simulation " << N << " of " << Ntot << std::endl;

            mu_WIMP=mu_WIMP_low+j*(mu_WIMP_high-mu_WIMP_low)/(nSteps_mu-1);
        
            // find tree
            TString tree_name;
            tree_name = TString::Format("%1.f_%1.f_%1.1f",m_WIMP,delta,mu_WIMP);
            
            if (!File->GetListOfKeys()->FindObject(tree_name)) {
                std::cout << "##### ERROR: found no tree named " << tree_name << std::endl;
                std::cout << "setting efficiency to zero!" << std::endl;
                efficiency = 0.;
            }
            else {
                TTree* dataTree = (TTree*) File->Get(tree_name);
                // get efficiency
                efficiency = GetEfficiency(dataTree);
                std::cout << "efficiency: " << efficiency << std::endl;
            }
            
            // fill histogram
            int ibin = hist->FindBin(delta,mu_WIMP);
            hist->SetBinContent(ibin,efficiency);
            
            //std::cout << delta << "\t" << mu_WIMP << "\t" << efficiency << std::endl;

        }
    }
    
    // draw results
    TCanvas* c = new TCanvas("c");
    gStyle->SetOptStat(0);
    hist->Draw("COLZ");
    c->SaveAs(fileName+"_results.pdf");
    c->SaveAs(fileName+"_results.root");
    
    File->Close();

    return 1;
}


double GetEfficiency(TTree* dataTree) {
    
    
    
    // minimum distance between S1s is about 0.5us according to https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=xenon:xenon100:analysis:mayra:run14:rn220:bipo
    
    // minimum distance between S2s is 3 mm/(1.73 m/s)=1.73us
    
    // set alias
    //
    // WIMP position at scattering
    dataTree->SetAlias("x0","x_WIMP0.X()");
    dataTree->SetAlias("y0","x_WIMP0.Y()");
    dataTree->SetAlias("z0","x_WIMP0.Z()");
    dataTree->SetAlias("r0","sqrt(x_WIMP0.X()**2+x_WIMP0.Y()**2)");
    
    // WIMP position at de-excitation
    dataTree->SetAlias("x1","x_WIMP1.X()");
    dataTree->SetAlias("y1","x_WIMP1.Y()");
    dataTree->SetAlias("z1","x_WIMP1.Z()");
    dataTree->SetAlias("r1","sqrt(x_WIMP1.X()**2+x_WIMP1.Y()**2)");
    
    // time of S1s and S2s
    dataTree->SetAlias("S1_2","tau");
    dataTree->SetAlias("S2_1","(0.153-z0)/1.73e3");
    dataTree->SetAlias("S2_2","tau+(0.153-z1)/1.73e3");

    
    // ---------------------------------------------------
    // cuts
    // ---------------------------------------------------

    // nuclear recoil in 34kg
    dataTree->SetAlias("XNR34kg","(TMath::Power(TMath::Abs(-z0+0.001)/0.1268,2.7)+TMath::Power((r0*r0)/0.0175,2.7)<1)"); // X34kg2
    
    // de-excitation in 34kg
    dataTree->SetAlias("XG34kg","(TMath::Power(TMath::Abs(-z0+0.001)/0.1268,2.7)+TMath::Power((r0*r0)/0.0175,2.7)<1)"); // X34kg2
    
    // NR in 48kg
    dataTree->SetAlias("XNR48kg","((pow(TMath::Abs((z0+0.003)/0.146),4)+pow(TMath::Abs(r0*r0/0.02),4))<1)"); // X48kg0
    
    // de-excitation in 48kg
    dataTree->SetAlias("XG48kg","((pow(TMath::Abs((z1+0.003)/0.146),4)+pow(TMath::Abs(r1*r1/0.02),4))<1)"); // X48kg0
    
    // second S1 before both S2s
    dataTree->SetAlias("XS1S2S2","(S1_2<S2_1) && (S1_2<S2_2)");

    // difference between S1s >=0.5us
    dataTree->SetAlias("XdtS1","S1_2>=0.5e-6");
    
    // difference between S2s >=1.73us
    dataTree->SetAlias("XdtS2","TMath::Abs(S2_2-S2_1)>=1.73e-6");


    
    // events with NR in 48kg fiducial volume
    double N_all = dataTree->Draw("","XNR48kg","goff");
    
    // events passing all cuts
    double N_cuts = dataTree->Draw("","XNR48kg && XG48kg && XS1S2S2 && XdtS1 && XdtS2","goff");
    
    // detection efficiency
    double efficiency = N_cuts/N_all;
    
    return efficiency;

}

