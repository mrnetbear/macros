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

void numIntMore(){

    //define step size
    const double stepSize = 1e-5;

    //define accurate integral value
    const double accInt = 1;

    double ctrInt, trpInt;
    std::thread th1(centralIntegrator, &ctrInt, stepSize);
    std::thread th2(trapezoidIntegrator, &trpInt, stepSize);
    th1.join();
    th2.join();

    double ctrErr[10000], trpErr[10000];
    double h1[10000];
    size_t i = 0;
    for (double h = stepSize; h < 1e-1 || i < 10000; h += stepSize){
        h1[i] = h;
        double ctr, trp;
        std::thread th1(centralIntegrator, &ctr, h);
        std::thread th2(trapezoidIntegrator, &trp, h);
        th1.join();
        th2.join();
        ctrErr[i] = abs(accInt - ctr);
        trpErr[i] = abs(accInt - trp);
        i++;
    }

    std::cout << "Central Rectangles Integration Result: " << ctrInt << std::endl;
    std::cout << "Trapezoidal Rule Integration Result: " << trpInt << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->SetGrid();
    c1->SetLogx();
    c1->SetLogy();

    TGraph *gr1 = new TGraph(i, h1, ctrErr);
    gr1->SetTitle("Central Rectangles Error vs Step Size");
    gr1->SetLineStyle(5);
    gr1->SetLineWidth(6);
    gr1->SetLineColor(kRed);

    TGraph *gr2 = new TGraph(i, h1, trpErr);
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