//
//  main.cpp
//  lab1
//
//  Created by Алексей on 16.02.2024.
//

#define Pi 3.1415926535
#define mPi 135.0
#define Epi 6.0e3
#define zpos 243.0

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <time.h>
#include <random>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom.h"

//#include "classes.cpp"

using namespace std;

const double Vpi = sqrt(1 - pow(mPi/Epi, 2));

class Hep3Vector{
public:
    Hep3Vector();
    Hep3Vector(double x, double y, double z);
    Hep3Vector(const Hep3Vector &v);
    double x();
    double y();
    double z();
    double phi();
    double cosTheta();
    double mag();
    void symmetry();
    void lorentz();
    void revLorentz();
private:
    double r, cos_Theta, Phi;
};

Hep3Vector::Hep3Vector(double x, double y, double z){
    r = x;
    cos_Theta = y;
    Phi = z;
}

Hep3Vector::Hep3Vector(const Hep3Vector &v) {
    r = v.r;
    Phi = v.Phi;
    cos_Theta = v.cos_Theta;
}

Hep3Vector::Hep3Vector(){
}

double Hep3Vector::phi(){
    return Phi;
}

double Hep3Vector::cosTheta(){
    return cos_Theta;
}

double Hep3Vector::mag(){
    return r;
}

double Hep3Vector::x(){
    return sin(acos(cos_Theta)) * cos(Phi) * r;
}

double Hep3Vector::y(){
    return sin(acos(cos_Theta)) * sin(Phi) * r;
}

double Hep3Vector::z(){
    return cos_Theta * r;
}

void Hep3Vector::symmetry(){
    cos_Theta *= -1;
    if (Phi > Pi)
        Phi -= Pi;
    else
        Phi += Pi;
}

void Hep3Vector::lorentz(){
    r = r / mPi * (Epi + sqrt(Epi*Epi - mPi*mPi) * cos_Theta);
    cos_Theta = (cos_Theta + Vpi) / (1 + Vpi * cos_Theta);
}

void Hep3Vector::revLorentz(){
    r = r / mPi * (Epi - sqrt(Epi*Epi - mPi*mPi) * cos_Theta);
    cos_Theta = (cos_Theta - Vpi) / (1 - Vpi * cos_Theta);
}

int decayPi0(){
    //file creaion
    TFile *oGrph = new TFile("decayPi0.root", "RECREATE");
    
    //MCS hists definition
    TH1* TE1 = new TH1D("E1", "E1", 500, 50, 80);
    TH1* TE2 = new TH1D("E2", "E2", 500, 50, 80);
    TH1* TE1E2 = new TH1D("E1-E2", "E1-E2", 500, -100, 100);
    TH1* TPhi1 = new TH1D("Phi1", "Phi1", 500, 0, 2*Pi+0.1);
    TH1* TPhi2 = new TH1D("Phi2", "Phi2", 500, 0, 2*Pi+0.1);
    TH1* TTheta1 = new TH1D("Theta1", "Theta1", 500, 0, Pi+0.1);
    TH1* TTheta2 = new TH1D("Theta2", "Theta2", 500, 0, Pi+0.1);
    TH1* TPx = new TH1D("SumPx", "Px1+Px2", 500, -100, 100);
    TH1* TPy = new TH1D("SumPy", "Py1+Py2", 500, -100, 100);
    TH1* TPz = new TH1D("SumPz", "Pz1+Pz2", 500, -100, 100);

    //Lab hists 
    TH1* TlE1 = new TH1D("lorentz E1", "lorentz E1", 500, -1000, 6500);
    TH1* TlE2 = new TH1D("lorentz E2", "lorentz E2", 500, -1000, 6500);
    TH1* TlE1E2 = new TH1D("lorentz total Energy", "lorentz total Energy", 500, -100, 100);
    TH1* TlPhi1 = new TH1D("loerntz Phi1", "loerntz Phi1", 500, -1, 2*Pi+1);
    TH1* TlPhi2 = new TH1D("loerntz Phi2", "loerntz Phi2", 500, -1, 2*Pi+1);
    TH1* TlTheta1 = new TH1D("lorentz Theta1", "lorentz Theta1", 500, -1, Pi+0.1);
    TH1* TlTheta2 = new TH1D("lorentz Theta2", "lorentz Theta1", 500, -1, Pi+0.1);
    TH1* TlPx = new TH1D("lorentz SumPx", "lorentz Px1+Px2", 500, -100, 100);
    TH1* TlPy = new TH1D("lorentz SumPy", "lorentz Py1+Py2", 500, -100, 100);
    TH1* TlPz = new TH1D("lorentz SumPz", "lorentz Pz1+Pz2", 500, -100, 100);

    TH2* Strike1 = new TH2D("strike1", "Particle 1 Strike", 1000, -3000, 3000, 1000, -3000, 3000);
    TH2* Strike2 = new TH2D("strike2", "Particle 2 Strike", 1000, -3000, 3000, 1000, -3000, 3000);

    TTree* det = new TTree("det_stat", "Detector statistics");
    double x1, y1, x2, y2, E1, E2;
    det -> Branch("x1", &x1, "x1/D");
    det -> Branch("y1", &y1, "y1/D");
    det -> Branch("E1", &E1, "E1/D");
    det -> Branch("x2", &x2, "x2/D");
    det -> Branch("y2", &y2, "y2/D");
    det -> Branch("E2", &E2, "E2/D");
    
    int n = 1E6;
    vector<Hep3Vector> gamma1, gamma2;
    TRandom *r1 = new TRandom();
    double p = 0, cos_Theta1= 0, cos_Theta2= 0, Phi1 = 0, Phi2 = 0;
    for (int i = 0; i < n; i++){
        p = mPi / 2;
        TE1->Fill(p);
        TE2->Fill(p);
        TE1E2->Fill(p-p);
        
        cos_Theta1 = r1 -> Uniform(-1,1);
        TTheta1->Fill(acos(cos_Theta1));
        TTheta2->Fill(Pi - acos(cos_Theta1));
        
        Phi1 = r1 -> Uniform(0, 2*Pi);
        TPhi1->Fill(Phi1);
        TPhi2->Fill(2*Pi - Phi1);
        
        Hep3Vector draft(p, cos_Theta1, Phi1);
        gamma1.push_back(draft);
        draft.symmetry();
        gamma2.push_back(draft);
        
        TPx -> Fill(gamma1[i].x() + gamma2[i].x());
        TPy -> Fill(gamma1[i].y() + gamma2[i].y());
        TPz -> Fill(gamma1[i].z() + gamma2[i].z());
        
        gamma1[i].lorentz();
        gamma2[i].lorentz();

        TlE1->Fill(gamma1[i].mag());
        TlE2->Fill(gamma2[i].mag());
        TlE1E2->Fill(gamma1[i].mag()+gamma2[i].mag()-Epi);
        
        TlTheta1->Fill(acos(gamma1[i].cosTheta()));
        TlTheta2->Fill(acos(gamma2[i].cosTheta()));
        
        TlPhi1->Fill(gamma1[i].phi());
        TlPhi2->Fill(gamma2[i].phi());
        
        TlPx -> Fill(gamma1[i].x() + gamma2[i].x());
        TlPy -> Fill(gamma1[i].y() + gamma2[i].y());
        TlPz -> Fill(gamma1[i].z() + gamma2[i].z() - sqrt(Epi*Epi - mPi*mPi));

        double t = zpos / gamma1[i].z();
        x1 = gamma1[i].x() * t;
        y1 = gamma1[i].y() * t;
        gamma1[i].z() > 0 ? E1 = gamma1[i].mag() : E1 = 0;

        t = zpos / gamma2[i].z();
        x2 = gamma2[i].x() * t;
        y2 = gamma2[i].y() * t;
        gamma2[i].z() > 0 ? E2 = gamma2[i].mag() : E2 = 0;

        if (E1)
            Strike1 -> Fill(x1, y1);
        if (E2)
            Strike2 -> Fill(x2, y2);

        det -> Fill();  
    } 

    
    TE1 -> Write();
    TE2 -> Write();
    TE1E2 -> Write();
    TPhi1 -> Write();
    TPhi2 -> Write();
    TTheta1 -> Write();
    TTheta2 -> Write();
    TPx -> Write();
    TPy -> Write();
    TPz -> Write();
    
    TlE1 -> Write();
    TlE2 -> Write();
    TlE1E2 -> Write();
    TlPhi1 -> Write();
    TlPhi2 -> Write();
    TlTheta1 -> Write();
    TlTheta2 -> Write();
    TlPx -> Write();
    TlPy -> Write();
    TlPz -> Write();

    Strike1 -> Write();
    Strike2 -> Write();

    det -> Write();

    oGrph->Close();
    
    
    return 0;
}

using namespace std;

