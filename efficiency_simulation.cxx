#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TF3.h>
#include <TFile.h>
#include <TUnuran.h>
#include <TUnuranMultiContDist.h>
#include <TRotation.h>

#include <fstream>
#include <time.h>

using namespace std;

// initialize random number generator
int seed = 41302;
TRandom3* fRandom = new TRandom3(seed);

// functions
TVector3 GetUniPos();
TVector3 GetV_scatt(TVector3 v_init);
void GetV_start(TF3* samplePDF, TVector3* Vin, double* E_R, double* cosTheta);
TVector3 GetV_scatt(TVector3 Vin, double E_R);
TVector3 RotateEarth(TVector3 v, double hour);
TVector3 SimpleRotation(TVector3 v, double phi, double theta);
double GetLifetime();
int ReadParameters(TString FileName);
int InitParameters(int j, int k);

// global variables
TF1* fGST;
TF1* fh_WIMP;
TF1* fA_WIMP;
TF1* fa_WIMP;
TF1* fphi_WIMP;
TF1* ftheta_WIMP;

// input parameters read from file
int fNParticles = -1;
double fm_WIMP = -1;
double fmu_WIMP_low = -1;
double fmu_WIMP_high = -1;
double fdelta_low = -1;
double fdelta_high = -1;
int fNSteps_delta = -1;
int fNSteps_mu = -1;

// parameters set by init function
double fmu_WIMP;
double fdelta;
double fv_min;
double fmu;

// common physics parameters
//
// nucleus mass (eV/c^2)
double fm_Nucleus = 122.2e9/(2.998e8*2.998e8);

// size of cylindrical detector (m)
double fR_det = 0.153;
double fh_det = 0.306;

// velocity of sun in galaxy (m/s)
double fv_sun = 232.e3;

// galactic escape velocity
double fv_esc = 550.e3;

// mean velocity of WIMPs in halo
double fv_mean = 220.e3;


// velocity distribution, x=v_WIMP, y=cos(Theta)
TF3* velocityPDF = new TF3("velocityPDF","(x/[0])**3*exp(-(x/[0])**2)*exp(-2.*x*[1]/([0]*[0])*y)+0*z");

// cross section, x=v_WIMP, z=E_R
TF3* cross_section = new TF3("cross_section","1/z*(1-z/(x*x)*(1/(2.*[0])+1./[1])-[2]/(x*x)*(([0]+[1])/([0]*[1])+[2]/(2.*[0]*z)))+0*y");

// form factor F(E_R)^2 (arXiv:1202.6073v2), z=E_R
TF3* FormFactor2 = new TF3("FormFactor2","(3.*exp(-[0]*z/2.)*(sin([1]*sqrt(z))-[1]*sqrt(z)*cos([1]*sqrt(z)))/([1]*sqrt(z))**3)**2+0*x+0*y");



// ---------------------------------------------------------
// main program
int main(int argc, char *argv[]) {
    
    // argc should be =2 for correct execution
    if ( argc != 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <efficiency_simulation_input.txt>" << std::endl;
        
        return 1;
    }
    
    TString FileName = argv[1];

    
    clock_t t1,t2;
    t1=clock();
    
    // read parameters from input file
    if (ReadParameters(FileName)) {
        cout << "######## ERROR: Could not read input parameters" << endl;
        return 1;
    }
    
    double m_WIMP = fm_WIMP*2.998e8*2.998e8/1.e9;
    double delta_low = fdelta_low/1.e3;
    double delta_high = fdelta_high/1.e3;
    double mu_WIMP_low = fmu_WIMP_low/1.e-3;
    double mu_WIMP_high = fmu_WIMP_high/1.e-3;
    
    // root file to save results
    TString results_name;
    results_name = TString::Format("%1.f_%1.f-%1.f_%1.1f-%1.1f",m_WIMP, delta_low,delta_high, mu_WIMP_low,mu_WIMP_high);
    TFile* file = new TFile("efficiency_"+results_name+".root","RECREATE");
    
    TTree* RunTree = new TTree("RunTree","RunTree");
    
    RunTree->Branch("m_WIMP",&m_WIMP);
    RunTree->Branch("nSteps_delta",&fNSteps_delta);
    RunTree->Branch("nSteps_mu",&fNSteps_mu);
    RunTree->Branch("delta_low",&delta_low);
    RunTree->Branch("delta_high",&delta_high);
    RunTree->Branch("mu_WIMP_low",&mu_WIMP_low);
    RunTree->Branch("mu_WIMP_high",&mu_WIMP_high);
    
    RunTree->Fill();
    RunTree->Write();
    
    int N=1;
    int Ntot=fNSteps_delta*fNSteps_mu;
    
    // loop over whole parameter space
    //for (int i=0; i<fNSteps_m; ++i) {     // fix wimp mass
        
        for (int j=0; j<fNSteps_delta; ++j) {

            for (int k=0; k<fNSteps_mu; ++k) {
                
                std::cout << "#### starting simulation nr. " << N << " of " << Ntot << std::endl;
                N++;

                // initialize parameters
                if (InitParameters(j,k)) {
                    continue;
                }
                
                // make TTree
                TString tree_name;
                tree_name = TString::Format("%1.f_%1.f_%1.1f",fm_WIMP*2.998e8*2.998e8/1.e9,fdelta/1.e3, fmu_WIMP/1.e-3);

                TTree* dataTree = new TTree(tree_name,tree_name);
                
                TVector3 x_WIMP0, v_WIMP0, x_WIMP1, v_WIMP1;
                
                double dt, E_R;
                
                dataTree->Branch("x_WIMP0", &x_WIMP0);
                dataTree->Branch("v_WIMP0", &v_WIMP0);
                dataTree->Branch("x_WIMP1", &x_WIMP1);
                dataTree->Branch("v_WIMP1", &v_WIMP1);
                dataTree->Branch("E_R", &E_R);
                dataTree->Branch("tau", &dt);

                
                double cosTheta;
                TVector3 v_WIMPin(1,1,1);
                
               
                // set sample PDF parameters
                velocityPDF->SetParameter(0,fv_mean);
                velocityPDF->SetParameter(1,fv_sun);
                
                cross_section->SetParameter(0,fm_Nucleus);
                cross_section->SetParameter(1,fm_WIMP);
                cross_section->SetParameter(2,fdelta);
                
                double mN = fm_Nucleus*(2.998e8*2.998e8);
                double s=1./0.197*1.e-9;
                double R=1.2*pow(131.2,1./3.)*s;
                
                FormFactor2->SetParameter(0,2.*mN*s*s);
                FormFactor2->SetParameter(1,sqrt(2.*mN)*sqrt(R*R-5.*s*s));
                
                // sample PDF, sigma*f(v)*F(E_R)
                TF3* samplePDF = new TF3("samplePDF","TMath::Max(0,velocityPDF*cross_section*FormFactor2)",fv_min,fv_esc+fv_sun,-1.,1.,0.,500.e3);
                
                // need to set a reasonable number of points in TF1 to get acceptable quality from GetRandom
                int np = 200;
                samplePDF->SetNpx(np);
                samplePDF->SetNpy(np);
                samplePDF->SetNpz(np);
                
                // loop over all particles
                for (int i=0; i<fNParticles; ++i) {
                    
                    if (i%(fNParticles/100)==0) {
                        cout << "\r" << "progress " << i*100./fNParticles << "%" << std::flush;
                    }
                    
                    // start position of WIMP (position of scattering)
                    x_WIMP0 = GetUniPos();
                    
                    //cout << "start position: " << x_WIMP0.X() << ", " << x_WIMP0.Y() << ", " << x_WIMP0.Z() << endl;
                    
                    
                    // velocity of WIMP before scattering
                    //
                    GetV_start(samplePDF, &v_WIMPin, &E_R, &cosTheta);
                    
                    //cout << "start velocity: " << v_WIMP0.X() << ", " << v_WIMP0.Y() << ", " << v_WIMP0.Z() << endl;
                    
                    // velocity of WIMP after scattering
                    TVector3 v_WIMPout = GetV_scatt(v_WIMPin, E_R);
                    
                    v_WIMP0 = RotateEarth(v_WIMPin,(i % 24));
                    v_WIMP1 = RotateEarth(v_WIMPout,(i % 24));
                    
                    //cout << "new velocity: " << v_WIMP1.X() << ", " << v_WIMP1.Y() << ", " << v_WIMP1.Z() << endl;
                    
                    // lifetime of WIMP
                    dt = GetLifetime();
                    
                    //cout << "lifetime: " << dt << endl;
                    
                    // position of de-excitation
                    x_WIMP1 = x_WIMP0+v_WIMP1*dt;
                    
                    //cout << "end position: " << x_WIMP1.X() << ", " << x_WIMP1.Y() << ", " << x_WIMP1.Z() << endl;
                    
                    
                    dataTree->Fill();
                    
                }
                
                cout << "\n" << endl;
                
                file->cd();
                dataTree->Write();
                
                // clean up
                delete dataTree;
                //delete velocityPDF;
                //delete cross_section;
                //delete FormFactor2;
                delete samplePDF;
                
            }
        }
    //}
    
    file->Close();
   
    t2=clock();
    
    cout << "run time: " << (t2-t1) / CLOCKS_PER_SEC << " sec." << endl;
    
    return 0;
    
}


// ---------------------------------------------------------
// get correctly distributed start values
void GetV_start(TF3* samplePDF, TVector3* Vin, double* E_R, double* cosTheta) {
    
    /*
    TUnuran unr;
    
    // Multi- dimensional case from a TF1 (TF2 or TF3) objects
    TUnuranMultiContDist dist(samplePDF);
    
    // the recommended method for multi-dimensional function is "hitro"
    unr.Init(dist,"method=hitro");
    */
    
    double x[3];
    
    // azimuth angle for incoming WIMP
    double phi = (fRandom->Uniform())*2.*TMath::Pi();
    
    // polar angel for incoming WIMP
    double theta;
    
    // velocity of incoming WIMP
    double v_WIMP;
    
    // recoil energy
    double ER;
    
    // kinematic constraints
    double ERmin, ERmax, cosTheta_min, cosTheta_max;
    
    // sample distribution
    do {
        
        //unr.SampleMulti(x);
        samplePDF->GetRandom3(x[0],x[1],x[2]);
        
        // velocity of incoming WIMP
        v_WIMP=x[0];
        
        // polar angel for incoming WIMP
        theta = acos(x[1]);
        
        // recoil energy
        ER=x[2];
        
        // check kinematic constraints
        ERmin=fmu*fmu*v_WIMP*v_WIMP/fm_Nucleus*(1-fdelta/(fmu*v_WIMP*v_WIMP)-sqrt(1-2.*fdelta/(fmu*v_WIMP*v_WIMP)));
        ERmax=fmu*fmu*v_WIMP*v_WIMP/fm_Nucleus*(1-fdelta/(fmu*v_WIMP*v_WIMP)+sqrt(1-2.*fdelta/(fmu*v_WIMP*v_WIMP)));
        
        cosTheta_min = -1;
        cosTheta_max = (fv_esc*fv_esc-fv_sun*fv_sun-v_WIMP*v_WIMP)/(2.*fv_sun*v_WIMP);
    }
    while(ER<ERmin||ER>ERmax||cos(theta)<cosTheta_min||cos(theta)>cosTheta_max);
    
        
    Vin->SetPhi(phi);
    Vin->SetTheta(theta);
    Vin->SetMag(v_WIMP);
    
    *cosTheta = x[1];
    *E_R = x[2];
    
}


// ---------------------------------------------------------
// get position uniformly distributed in cylinder
TVector3 GetUniPos() {
    
    double phi = (fRandom->Uniform())*2.*TMath::Pi();
    double rho = sqrt(fRandom->Uniform())*fR_det;
    
    double xPos = rho*cos(phi);
    double yPos = rho*sin(phi);
    double zPos = (fRandom->Uniform()-0.5)*fh_det;
    
    TVector3 Pos(xPos,yPos,zPos);
    
    return Pos;
    
}

// ---------------------------------------------------------
// calculate WIMP velocity vector after scattering
TVector3 GetV_scatt(TVector3 Vin, double E_R) {
    
    // velocity of incoming WIMP
    double v = Vin.Mag();
    
    // azimuth angle for incoming WIMP
    double phi = Vin.Phi();
    
    // polar angle for incoming WIMP
    double theta = Vin.Theta();
    
    // azimuth angle of outgoing WIMP in COM frame
    double phi_COM = (fRandom->Uniform())*2.*TMath::Pi();
    
    // polar angle of outgoing WIMP in COM frame
    double PSfactor = sqrt(1-2.*fdelta/(fmu*v*v));
    double theta_COM = acos(1./PSfactor*((fm_Nucleus*E_R)/(fmu*fmu*v*v)-1+fdelta/(fmu*v*v)));
    
    // velocity of nuclear recoil in COM frame
    //TVector3 v_nuc_COM(1,1,1);
    //v_nuc_COM.SetMag(fmu*v/fm_Nucleus*PSfactor);
    //v_nuc_COM.SetPhi(phi_COM);
    //v_nuc_COM.SetTheta(theta_COM);
    
    // outgoing velocity of WIMP in COM frame
    TVector3 v_WIMP_COM(1,1,1);
    v_WIMP_COM.SetMag(fmu*v/fm_WIMP*PSfactor);
    v_WIMP_COM.SetPhi(phi_COM);
    v_WIMP_COM.SetTheta(theta_COM);
    v_WIMP_COM = -v_WIMP_COM;
    
    // boost back to the earth frame
    TVector3 v_boost(0,0,fmu/fm_Nucleus*v);
    //Tvector3 v_nuc = v_nuc_COM + v_boost;
    TVector3 v_WIMPout = v_WIMP_COM + v_boost;
    
    // rotate to the right frame
    v_WIMPout = SimpleRotation(v_WIMPout,-phi,-theta);
    
    return v_WIMPout;
    
}

// ---------------------------------------------------------
// rotate vector
TVector3 SimpleRotation(TVector3 v, double phi, double theta) {
    
    TRotation r;
    r.SetYEulerAngles(phi,theta,phi);
    r.Invert();
    
    v*=r;
    
    return v;
    
    
}

// ---------------------------------------------------------
// rotate vectors to the lab frame
TVector3 RotateEarth(TVector3 v, double hour) {
    
    double tlab = 2.*TMath::Pi()/360.*(15.*(hour+13.7/15.));
    double llab = 42.45*2.*TMath::Pi()/360.;
    
    TVectorD vd(3);
    
    // at the moment the average velocity is in the z-direction
    // to get galactic coordinates cycle {x,y,z} → {y,z,x}
    vd(0)=v(1);
    vd(1)=v(2);
    vd(2)=v(0);

    TMatrixD r1(3,3);
    TMatrixD r2(3,3);
    
    r1(0,0)=-0.06699;
    r1(0,1)=0.4927;
    r1(0,2)=-0.8676;
    r1(1,0)=-0.8728;
    r1(1,1)=-0.4503;
    r1(1,2)=-0.1884;
    r1(2,0)=-0.4835;
    r1(2,1)=0.7446;
    r1(2,2)=0.4602;
    
    r2(0,0)=-sin(llab)*cos(tlab);
    r2(0,1)=-sin(llab)*sin(tlab);
    r2(0,2)=cos(llab);
    r2(1,0)=sin(tlab);
    r2(1,1)=-cos(tlab);
    r2(1,2)=0;
    r2(2,0)=cos(llab)*cos(tlab);
    r2(2,1)=cos(llab)*sin(tlab);
    r2(2,2)=sin(llab);
    
    vd*=r1;
    vd*=r2;
    
    TVector3 v1(vd(0),vd(1),vd(2));
    
    return v1;
    
    
}

// ---------------------------------------------------------
// calculate lifetime (sec.) of excited WIMP
double GetLifetime() {
    
    double tau = TMath::Pi()*pow(6.626e-34/(2.*TMath::Pi()),3)*pow(5.29e-19,2)*pow(2.998e8,4)/(pow(fdelta*1.6022e-19,3)*pow(fmu_WIMP*5.0508e-27,2));
    double lifetime = fRandom->Exp(tau);
    
    return lifetime;
    
    //cout << "tau: " << tau << endl;
    //cout << "lifetime: " << lifetime << endl;
    
}

// ---------------------------------------------------------
// read input parameters from file
int ReadParameters(TString FileName) {
    
    ifstream File;
    File.open(FileName);
    
    if(!File.is_open()) {
        cout << "ERROR: Could not open " << FileName << endl;
        return 0;
    }
    
    string headerline;
    double m_WIMP, delta_low, delta_high, mu_WIMP_low, mu_WIMP_high;
    int NParticles, nSteps_delta, nSteps_mu;
    
    getline(File, headerline);
    File >> NParticles;
    getline(File, headerline);
    getline(File, headerline);
    
    File >> m_WIMP;
    getline(File, headerline);
    getline(File, headerline);

    File >> delta_low;
    getline(File, headerline);
    getline(File, headerline);
    
    File >> delta_high;
    getline(File, headerline);
    getline(File, headerline);
    
    File >> nSteps_delta;
    getline(File, headerline);
    getline(File, headerline);

    File >> mu_WIMP_low;
    getline(File, headerline);
    getline(File, headerline);
    
    File >> mu_WIMP_high;
    getline(File, headerline);
    getline(File, headerline);
    
    File >> nSteps_mu;
    getline(File, headerline);
    getline(File, headerline);
    
    File.close();
    
    
    fNParticles=NParticles;
    fm_WIMP=m_WIMP/(2.998e8*2.998e8);
    fdelta_low=delta_low;
    fdelta_high=delta_high;
    fNSteps_delta=nSteps_delta;
    fmu_WIMP_low=mu_WIMP_low;
    fmu_WIMP_high=mu_WIMP_high;
    fNSteps_mu=nSteps_mu;
    
    if (fNParticles!=-1&&fm_WIMP!=-1&&fdelta_low!=-1&&fdelta_high!=-1&&fmu_WIMP_low!=-1&&fmu_WIMP_high!=-1&&fNSteps_delta!=-1&&fNSteps_mu!=-1) {
        return 0;
    }
    
    else {
        return 1;
    }
    
}

// ---------------------------------------------------------
// initialize parameters
int InitParameters(int step_delta, int step_mu) {
    
    //fm_WIMP=fm_WIMP_low+step_m*(fm_WIMP_high-fm_WIMP_low)/(fNSteps_m-1);

    fdelta=fdelta_low+step_delta*(fdelta_high-fdelta_low)/(fNSteps_delta-1);
    
    fmu_WIMP=fmu_WIMP_low+step_mu*(fmu_WIMP_high-fmu_WIMP_low)/(fNSteps_mu-1);

    // reduced mass
    fmu=fm_WIMP*fm_Nucleus/(fm_WIMP+fm_Nucleus);
    
    // minimum velocity due to inelastic scattering
    fv_min = sqrt(2.*fdelta/fmu);
    
    if  (fv_esc+fv_sun<fv_min) {
        cout << "ERROR: Inelastic scattering is kinematically not possible! Choose a different combination of m_WIMP and delta!" << endl;
        return 1;
    }
    else {
        return 0;
    }
    
}



