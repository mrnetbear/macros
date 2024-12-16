#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <random>
#include <thread>
#include <chrono>
#include <math.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TDecompLU.h"
#include "TMultiGraph.h"

//define step size
const double stepSize = 1e-5;

//define accurate integral value
const double accInt = 1;

//define function 
double f(double x){
        return TMath::Exp(-x);
}

//define quasi-uniform grid
double qGrid(double h){
        return h/pow(1-h,3);
}

//define central rectangles integration function
void centralIntegrator(double *result, double stepSize){
    double integral = 0;
    for(double xi = stepSize; xi <= 1; xi += stepSize){
        double h = qGrid(xi) - qGrid(xi-stepSize);
        integral += f(qGrid(xi)-h/2)*h;
    }
    *result = integral;
        
}

//efine trapezoid integration function
void trapezoidIntegrator(double *result, double stepSize){
    double integral = 0;
    for(double xi = stepSize; xi <= 1; xi += stepSize){
        double h = qGrid(xi) - qGrid(xi-stepSize);
        integral += (f(qGrid(xi))+f(qGrid(xi-stepSize)))/2 * h;
    }
    *result = integral;
        
}

void GKIntegrator2(double *C, double *adX1, double stepSize){
    double  C0 = 0,
            C1 = 0;
    for(double xi = stepSize; xi <= 1; xi += stepSize){
        double h = qGrid(xi) - qGrid(xi-stepSize);
        C0 += (qGrid(xi-stepSize) + 0.5*h - adX1[1])/(adX1[0] - adX1[1])*f(qGrid(xi-stepSize)+0.5*h)*h;
        C1 += (qGrid(xi-stepSize) + 0.5*h - adX1[0])/(adX1[1] - adX1[0])*f(qGrid(xi-stepSize)+0.5*h)*h;
    }
    C[0] = C0;
    C[1] = C1;
}

void GKIntegrator3(double *C, double *adX2, double stepSize){
    double  C0 = 0,
            C1 = 0,
            C2 = 0;
    for(double xi = stepSize; xi <= 1; xi += stepSize){
        double h = qGrid(xi) - qGrid(xi-stepSize);
        C0 += (qGrid(xi-stepSize) + 0.5*h - adX2[1]) * (qGrid(xi-stepSize) + 0.5*h - adX2[2]) / (adX2[0] - adX2[1]) / (adX2[0] - adX2[2])*f(qGrid(xi-stepSize)+0.5*h)*h;
        C1 += (qGrid(xi-stepSize) + 0.5*h - adX2[0]) * (qGrid(xi-stepSize) + 0.5*h - adX2[2]) / (adX2[1] - adX2[0]) / (adX2[1] - adX2[2])*f(qGrid(xi-stepSize)+0.5*h)*h;
        C2 += (qGrid(xi-stepSize) + 0.5*h - adX2[0]) * (qGrid(xi-stepSize) + 0.5*h - adX2[1]) / (adX2[2] - adX2[0]) / (adX2[2] - adX2[1])*f(qGrid(xi-stepSize)+0.5*h)*h;
    }
    C[0] = C0;
    C[1] = C1;
    C[2] = C2;
}

void numIntMore(){

    //define step size
    const double fStepSize = 1e-5;

    //define accurate integral value
    const double dAccInt = 1;

    double dCtrInt, dTrpInt;
    std::thread th1(centralIntegrator, &dCtrInt, fStepSize);
    std::thread th2(trapezoidIntegrator, &dTrpInt, fStepSize);
    th1.join();
    th2.join();

    double adCtrErr[10000], adTrpErr[10000];
    double adH1[10000];
    size_t i = 0;
    for (double h = fStepSize; h < 1e-1 || i < 10000; h += fStepSize){
        adH1[i] = h;
        double ctr, trp;
        std::thread th1(centralIntegrator, &ctr, h);
        std::thread th2(trapezoidIntegrator, &trp, h);
        th1.join();
        th2.join();
        adCtrErr[i] = abs(dAccInt - ctr);
        adTrpErr[i] = abs(dAccInt - trp);
        i++;
    }

    //define Gauss-Kristoffel algorithm
    const int   iNum0 = 1,
                iNum1 = 2,
                iNum2 = 3;

    double  adX0[iNum0] = {1},
            adX1[iNum1] = {2 - sqrt(2), 2 + sqrt(2)},
            adX2[iNum2] = {0.4158, 2.2942, 6.2899};
    
    double  adC0[iNum0] = {dCtrInt},
            adC1[iNum1],
            adC2[iNum2] = {0, 0, 0};

    GKIntegrator2(adC1, adX1, fStepSize);
    GKIntegrator3(adC2, adX2, fStepSize);

    double  dGKInt2 = adC1[0] + adC1[1],
            dGKInt3 = adC2[0] + adC2[1] + adC2[2];
    
    std::cout << "Central Rectangles Integration Result: " << dCtrInt << std::endl;
    std::cout << "Trapezoidal Rule Integration Result: " << dTrpInt << std::endl;
    std::cout << "Gauss-Kristoffel 0 Integration Result: " << adC0[0] << std::endl;
    std::cout << "Gauss-Kristoffel 1 Integration Result: " << dGKInt2 << std::endl;
    std::cout << "Gauss-Kristoffel 2 Integration Result: " << dGKInt3 << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->SetGrid();
    c1->SetLogx();
    c1->SetLogy();

    TGraph *gr1 = new TGraph(i, adH1, adCtrErr);
    gr1->SetTitle("Central Rectangles Error vs Step Size");
    gr1->SetLineStyle(5);
    gr1->SetLineWidth(6);
    gr1->SetLineColor(kRed);

    TGraph *gr2 = new TGraph(i, adH1, adTrpErr);
    gr2->SetTitle("Trapezoidal Rule Error vs Step Size");
    gr2->SetLineStyle(2);
    gr2->SetLineWidth(6);
    gr2->SetLineColor(kBlue);

    TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg->AddEntry(gr1, "Central Rectangles", "L");
    leg->AddEntry(gr2, "Trapezoidal Rule", "L");

    TMultiGraph *mgr = new TMultiGraph();
    mgr->Add(gr1);
    mgr->Add(gr2);
    mgr->SetTitle("Error vs Step Size");

    mgr->Draw("ALP");
    leg->Draw("SAME");
}