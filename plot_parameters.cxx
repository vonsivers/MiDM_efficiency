#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>

double GetEfficiency(TTree* tree);
TH1D* GetdtS1(TTree* tree);
TH1D* GetdtS2(TTree* tree);
TH1D* GetdtS1S2(TTree* tree);
TH1D* Getdx(TTree* tree);

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
    
    // create histos
    TH2D* hist = new TH2D("hist",";#delta (keV);#mu_{#chi} (#times10^{-3}#mu_{N})",nSteps_delta,bin_delta_low,bin_delta_high,nSteps_mu,bin_mu_low,bin_mu_high);
    hist->SetTitle("Efficiency");
    double efficiency;
    
    TH2D* hist2D_S1 = new TH2D("hist2D_S1",";#delta (keV);#mu_{#chi} (#times10^{-3}#mu_{N})",nSteps_delta,bin_delta_low,bin_delta_high,nSteps_mu,bin_mu_low,bin_mu_high);
    hist2D_S1->SetTitle("Time Difference S1_{2}-S1_{1}");
    double dS1;
    
    TH2D* hist2D_S2 = new TH2D("hist2D_S2",";#delta (keV);#mu_{#chi} (#times10^{-3}#mu_{N})",nSteps_delta,bin_delta_low,bin_delta_high,nSteps_mu,bin_mu_low,bin_mu_high);
    hist2D_S2->SetTitle("Time Difference S2_{2}-S2_{1}");
    double dS2;
    
    TH2D* hist2D_S1S2 = new TH2D("hist2D_S1S2",";#delta (keV);#mu_{#chi} (#times10^{-3}#mu_{N})",nSteps_delta,bin_delta_low,bin_delta_high,nSteps_mu,bin_mu_low,bin_mu_high);
    hist2D_S1S2->SetTitle("Time Difference S2_{1}-S1_{2}");
    double dS1S2;
    
    TH2D* hist2D_dx = new TH2D("hist2D_dx",";#delta (keV);#mu_{#chi} (#times10^{-3}#mu_{N})",nSteps_delta,bin_delta_low,bin_delta_high,nSteps_mu,bin_mu_low,bin_mu_high);
    hist2D_dx->SetTitle("Distance x_{2}-x_{1}");
    double dx;
    
    TH1D* hist_S1 = new TH1D("hist_S1",";Time Difference S1_{2}-S1_{1} (s);Counts (a.u.)",100,0.,10.e-6);
    TH1D* hist_S2 = new TH1D("hist_S2",";Time Difference S2_{2}-S2_{1} (s);Counts (a.u.)",1000,-500.e-6,500.e-6);
    TH1D* hist_S1S2 = new TH1D("hist_S1S2",";Time Difference S2_{1}-S1_{2} (s);Counts (a.u.)",1000,-500.e-6,500.e-6);
    TH1D* hist_dx = new TH1D("hist_dx",";Distance x_{2}-x_{1} (m);Counts (a.u.)",1000,0.,1.);
    
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
                dx = 0;
                dS1 = 0;
                dS2 = 0;
                dS1S2 = 0;
            }
            else {
                TTree* dataTree = (TTree*) File->Get(tree_name);
                // get efficiency
                efficiency = GetEfficiency(dataTree);
                std::cout << "efficiency: " << efficiency << std::endl;
                
                dS1 = GetdtS1(dataTree)->GetMean();
                dS2 = GetdtS2(dataTree)->GetMean();
                dS1S2 = GetdtS1S2(dataTree)->GetMean();
                dx = Getdx(dataTree)->GetMean();
                
                hist_dx->Add(Getdx(dataTree));
                hist_S1->Add(GetdtS1(dataTree));
                hist_S2->Add(GetdtS2(dataTree));
                hist_S1S2->Add(GetdtS1S2(dataTree));
                
            }
            
            // fill histogram
            int ibin = hist->FindBin(delta,mu_WIMP);
            hist->SetBinContent(ibin,efficiency);
            hist2D_S1->SetBinContent(ibin,dS1);
            hist2D_S2->SetBinContent(ibin,dS2);
            hist2D_S1S2->SetBinContent(ibin,dS1S2);
            hist2D_dx->SetBinContent(ibin,dx);
 
            //std::cout << delta << "\t" << mu_WIMP << "\t" << efficiency << std::endl;

        }
    }
    
    // draw results
    TCanvas* c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    hist->Draw("COLZ");
    c1->SaveAs(fileName+"_results.pdf");
    c1->SaveAs(fileName+"_results.root");
    
    // draw results
    TCanvas* c2 = new TCanvas("c2");
    gStyle->SetOptStat(0);
    hist_dx->Draw("");
    c2->SaveAs(fileName+"_dx.pdf");
    c2->SaveAs(fileName+"_dx.root");
    
    // draw results
    TCanvas* c3 = new TCanvas("c3");
    gStyle->SetOptStat(0);
    hist_S1->Draw("");
    c3->SaveAs(fileName+"_dS1.pdf");
    c3->SaveAs(fileName+"_dS1.root");
    
    // draw results
    TCanvas* c4 = new TCanvas("c4");
    gStyle->SetOptStat(0);
    hist_S2->Draw("");
    c4->SaveAs(fileName+"_dS2.pdf");
    c4->SaveAs(fileName+"_dS2.root");
    
    // draw results
    TCanvas* c5 = new TCanvas("c5");
    gStyle->SetOptStat(0);
    hist_S1S2->Draw("");
    c5->SaveAs(fileName+"_dS1S2.pdf");
    c5->SaveAs(fileName+"_dS1S2.root");
    
    // draw results
    TCanvas* c6 = new TCanvas("c6");
    gStyle->SetOptStat(0);
    hist2D_S1->Draw("COLZ");
    c6->SaveAs(fileName+"_2D_dS1.pdf");
    c6->SaveAs(fileName+"_2D_dS1.root");
    
    // draw results
    TCanvas* c7 = new TCanvas("c7");
    gStyle->SetOptStat(0);
    hist2D_S2->Draw("COLZ");
    c7->SaveAs(fileName+"_2D_dS2.pdf");
    c7->SaveAs(fileName+"_2D_dS2.root");
    
    // draw results
    TCanvas* c8 = new TCanvas("c8");
    gStyle->SetOptStat(0);
    hist2D_S1S2->Draw("COLZ");
    c8->SaveAs(fileName+"_2D_dS1S2.pdf");
    c8->SaveAs(fileName+"_2D_dS1S2.root");
    
    // draw results
    TCanvas* c9 = new TCanvas("c9");
    gStyle->SetOptStat(0);
    hist2D_dx->Draw("COLZ");
    c9->SaveAs(fileName+"_2D_dx.pdf");
    c9->SaveAs(fileName+"_2D_dx.root");

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
    
    // time difference between signals
    dataTree->SetAlias("dS1","S1_2");
    dataTree->SetAlias("dS2","S2_2-S2_1");
    
    // position difference
    dataTree->SetAlias("dx","Abs(sqrt(x1*x1+y1*y1+z1*z1)-sqrt(x0*x0+y0*y0+z0*z0))");

    
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

TH1D* Getdx(TTree* dataTree) {
    
    if(TH1D* hdx = (TH1D*)gROOT->FindObject("hdx")) {
        delete hdx;
    }
    
    TH1D* hdx = new TH1D("hdx",";#Delta x (m);Counts (a.u.)",1000,0.,1.);
    
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
    
    // position difference
    dataTree->SetAlias("dx","sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))");
    
    
    // ---------------------------------------------------
    // cuts
    // ---------------------------------------------------
    
    // NR in 48kg
    dataTree->SetAlias("XNR48kg","((pow(TMath::Abs((z0+0.003)/0.146),4)+pow(TMath::Abs(r0*r0/0.02),4))<1)"); // X48kg0
    
    // de-excitation in 48kg
    dataTree->SetAlias("XG48kg","((pow(TMath::Abs((z1+0.003)/0.146),4)+pow(TMath::Abs(r1*r1/0.02),4))<1)"); // X48kg0
    
    dataTree->Draw("dx>>hdx","XNR48kg && XG48kg","goff");
    

    return hdx;
    
}

TH1D* GetdtS1(TTree* dataTree) {
    
    if(TH1D* hS1 = (TH1D*)gROOT->FindObject("hS1")) {
        delete hS1;
    }

    
    TH1D* hS1 = new TH1D("hS1",";#Delta t_{S1} (s);Counts (a.u.)",100,0.,10.e-6);
    
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
    
    // time difference between signals
    dataTree->SetAlias("dS1","S1_2");
    dataTree->SetAlias("dS2","S2_2-S2_1");
    
    
    // ---------------------------------------------------
    // cuts
    // ---------------------------------------------------
    
    // NR in 48kg
    dataTree->SetAlias("XNR48kg","((pow(TMath::Abs((z0+0.003)/0.146),4)+pow(TMath::Abs(r0*r0/0.02),4))<1)"); // X48kg0
    
    // de-excitation in 48kg
    dataTree->SetAlias("XG48kg","((pow(TMath::Abs((z1+0.003)/0.146),4)+pow(TMath::Abs(r1*r1/0.02),4))<1)"); // X48kg0
    
    dataTree->Draw("dS1>>hS1","XNR48kg && XG48kg","goff");
    
    return hS1;
    
}

TH1D* GetdtS2(TTree* dataTree) {
    
    if(TH1D* hS2 = (TH1D*)gROOT->FindObject("hS2")) {
        delete hS2;
    }

    
    TH1D* hS2 = new TH1D("hS2",";#Delta t_{S2} (s);Counts (a.u.)",1000,-500.e-6,500.e-6);
    
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
    
    // time difference between signals
    dataTree->SetAlias("dS1","S1_2");
    dataTree->SetAlias("dS2","S2_2-S2_1");
    
    
    // ---------------------------------------------------
    // cuts
    // ---------------------------------------------------
    
    // NR in 48kg
    dataTree->SetAlias("XNR48kg","((pow(TMath::Abs((z0+0.003)/0.146),4)+pow(TMath::Abs(r0*r0/0.02),4))<1)"); // X48kg0
    
    // de-excitation in 48kg
    dataTree->SetAlias("XG48kg","((pow(TMath::Abs((z1+0.003)/0.146),4)+pow(TMath::Abs(r1*r1/0.02),4))<1)"); // X48kg0
    
    dataTree->Draw("dS2>>hS2","XNR48kg && XG48kg","goff");
    
    return hS2;
    
}

TH1D* GetdtS1S2(TTree* dataTree) {
    
    if(TH1D* hS1S2 = (TH1D*)gROOT->FindObject("hS1S2")) {
        delete hS1S2;
    }
    
    TH1D* hS1S2 = new TH1D("hS1S2",";#Delta t_{S1S2} (s);Counts (a.u.)",1000,-500.e-6,500.e-6);
    
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
    
    // time difference between signals
    dataTree->SetAlias("dS1S2","S2_1-S1_2");
    
    
    // ---------------------------------------------------
    // cuts
    // ---------------------------------------------------
    
    // NR in 48kg
    dataTree->SetAlias("XNR48kg","((pow(TMath::Abs((z0+0.003)/0.146),4)+pow(TMath::Abs(r0*r0/0.02),4))<1)"); // X48kg0
    
    // de-excitation in 48kg
    dataTree->SetAlias("XG48kg","((pow(TMath::Abs((z1+0.003)/0.146),4)+pow(TMath::Abs(r1*r1/0.02),4))<1)"); // X48kg0
    
    dataTree->Draw("dS1S2>>hS1S2","XNR48kg && XG48kg","goff");
    
    return hS1S2;
    
}

