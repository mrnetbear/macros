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

double maxValue(std::pair<double, double> arr){
    if (arr.first > arr.second) return arr.first;
    else return arr.second;
}
double maxValue(TMatrixD &arr){
    double max = -INFINITY;
    for (int i = 0; i < arr.GetNrows(); i++) {
        for (int j = 0; j < arr.GetNcols(); j++) {
            if (abs(arr(i, j)) > max) max = abs(arr(i, j));
        }
    }
    return max;
}
void simpleIterations(){
    // Define accuracy
    const double dEpsilon = 1e-6;

    // Define initial guess for x
    std::pair<double, double> pdX; // nitial guess for x
    std::vector<std::pair<double, double>> vpdSolutions; // every guess for solution
    pdX.first = 0;
    pdX.second = 0;
    vpdSolutions.push_back(pdX);
    int  i = 0;
    // Iterate until accuracy is reached or maximum number of iterations is reached (1000)
    while (i < 1000){
        i++;
        std::pair <double,  double> pdXTemp = {exp(pdX.first * pdX.first - 1), pdX.first * pdX.first - 1}, // Temporary solution
                                    pdDelta = {abs(pdX.first - pdXTemp.first), abs(pdX.second - pdXTemp.second)}; // Delta
        vpdSolutions.push_back(pdXTemp);
        pdX.first = pdXTemp.first;
        pdX.second = pdXTemp.second;
        if (maxValue(pdDelta) < dEpsilon) break; // Check for epsilon
    }

    // Print results
    std::cout << "========Simple Iterations========" << std::endl;
    std::cout << "Final values:\tx = " << pdX.first << ", y = " << pdX.second << std::endl;
    std::cout << "Iterations:\t" << i << std::endl;
    std::cout << "=================================" << std::endl;

    // Plot results
    TCanvas *c1 = new TCanvas("c1", "Simple Iterations", 800, 600);
    c1->SetGrid();
    
    TGraph *sol = new TGraph(vpdSolutions.size(), &vpdSolutions[0].first, &vpdSolutions[0].second);
    sol->SetMarkerStyle(20);
    sol->SetMarkerColor(kBlue);
    sol->SetLineWidth(0);

    TF1 *f1 = new TF1("f1", "x^2-1", -2, 2);
    TF1 *f2 = new TF1("f2", "log(x)", 0, 2);
    f1->Draw("");
    f2->Draw("SAME");
    //TMultiGraph *g = new TMultiGraph();
    //g->Add(sol);
    //g->Add(f1);
    //g->Add(f2);

    sol->Draw("P");
}

void NewtonIterations(){
    // Define accuracy
    const double dEpsilon = 1e-6;

    // Define initial guess for x
    TMatrixD pdX (2, 1);

    pdX(0, 0) = 0.1;
    pdX(1, 0) = -1;

    //std::vector<TMatrixD> vpdSolutions;
    //vpdSolutions.push_back(pdX);
    std::vector <std::pair <double, double>> vpdSolutions;
    vpdSolutions.push_back(std::make_pair(pdX(0, 0), pdX(1, 0)));

    int  i = 0;
    while (i < 1000){
        i++;
        TMatrixD pdF (2, 1);
        pdF(0, 0) = pdX(0, 0) * pdX(0, 0) - 1 - pdX(1, 0);
        pdF(1, 0) = log(pdX(0, 0)) - pdX(1, 0);

        // Calculation of Jacobian
        TMatrixD pdJacobian (2, 2);
        pdJacobian(0, 0) = 2 * pdX(0, 0);
        pdJacobian(0, 1) = -1;
        pdJacobian(1, 0) = 1 / pdX(0, 0);
        pdJacobian(1, 1) = -1;

        // Calculation of delta
        TMatrixD pdDelta (2, 1);
        pdDelta = pdJacobian.Invert() * pdF;

        pdX -= pdDelta;
        vpdSolutions.push_back(std::make_pair(pdX(0, 0), pdX(1, 0)));
        
        if (maxValue(pdDelta) < dEpsilon) break;
    }

    TCanvas *c2 = new TCanvas("c2", "Newton Iterations", 800, 600);
    c2->SetGrid();
    
    TGraph *sol = new TGraph(vpdSolutions.size(), &vpdSolutions[0].first, &vpdSolutions[0].second);
    sol->SetMarkerStyle(20);
    sol->SetMarkerColor(kBlue);
    sol->SetLineWidth(0);

    TF1 *f1 = new TF1("f1", "x^2-1", -2, 2);
    TF1 *f2 = new TF1("f2", "log(x)", 0, 2);
    f1->Draw("");
    f2->Draw("SAME");
    sol->Draw("P");

    // Print results
    std::cout << "========Newton Iterations==========" << std::endl;
    std::cout << "Final values:\tx = " << pdX(0, 0) << ", y = " << pdX(1, 0) << std::endl;
    std::cout << "Iterations:\t" << i << std::endl;
    std::cout << "=================================" << std::endl;
}