#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
//#include <math.h>
#include <cmath>
#include <time.h>
#include <random>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TDecompLU.h"

int interPoly(){

    //Define the coefficients of the polynomial
    const int n = 7;
    TMatrixD A(n,n);
    TMatrixD f(n,1);
    TMatrixD aj(n,1);
    TMatrixD x(n,1);

    //Set the coefficients

    const double left = 0, right = 3;

    for (int i=0; i<n; i++){
        A(i,0) = 1;
        x(i,0) = right/(double)n * i;

        f(i,0) = sin(x(i,0))*exp(-x(i,0));
        for (int j=1; j<n; j++){
            A(i,j) = pow(x(i,0), j);
        }
    }
    //Print all coefficients
    std::cout << "================================================================" << std::endl;
    for (int i=0; i<n; i++){
        for (int j = 0; j < n; j++)
        {
            std::cout << A(i,j) << " ";
        }
        std::cout << "| " << x(i,0) << " | " << f(i,0) << std::endl;
        
    }
    std::cout << "================================================================" << std::endl;
    //Perform the transformation
    double det = A.Determinant();
    if (!det) return 1;

    aj = A.Invert(&det)*f;

    //Print inverted matrix
    std::cout << "================Inverted matrix================" << std::endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
        {
            std::cout << A(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "===============================================" << std::endl;

    //Print found coefficients
    std::cout << "================Original coefficients================" << std::endl;
    for (int i = 0; i < n; i++){
        std::cout << "aj(" << i << ") = " << aj(i,0) << std::endl;
    }
    std::cout << "=====================================================" << std::endl;

    //Generate Lagrange coefficients
    double lagrangeCoef[n];
    for (int i = 0; i < n; i++) {
        lagrangeCoef[i] = 1.0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                lagrangeCoef[i] *= (x(i,0) - x(j,0));
            }
        }
    }

    //Generate original function
    double x1 [100];
    double fx1 [100];
    for (int i = 0; i < 100; i++){
        x1[i] = 3.0/100.00 * i;
        fx1[i] = sin(x1[i])*exp(-x1[i]);
    }

    

    //Generate approximated based Polynomial function 
    double x2 [100];
    double fx2 [100];
    for (int i = 0; i < 100; i++){
        x2[i] = 3.0/100.00 * i;
        for (int j = 0; j < n; j++){
            fx2[i] += aj(j,0)*pow(x2[i], j);
        }
    }

    //Calculate Lagrange polynomial
    double x3 [100];
    double fx3 [100];
    for (int i = 0; i < 100; i++){
        x3[i] = 3.0/100.00 * i;
        fx3[i] = 0.0;
        for (int j = 0; j < n; j++){
            double numerator = 1.0;
            for (int k = 0; k < n; k++){
                if (k != j) {
                    numerator *= (x3[i] - x(k,0));
                }
            }
            fx3[i] += (numerator / lagrangeCoef[j]) * f(j,0);
        }
    }

    //Draw original function
    TGraph* origin = new TGraph (100, x1, fx1);
    origin->SetMarkerStyle(20);
    origin->SetMarkerColor(kRed);
    origin->SetLineColor(kRed);
    origin->SetLineWidth(2);

    //Draw based approximation
    TGraph* approxPol = new TGraph (100, x2, fx2);
    approxPol->SetMarkerStyle(30);
    approxPol->SetMarkerColor(kBlue);
    approxPol->SetLineColor(kBlue);
    approxPol->SetLineWidth(4);

    //Draw Lagrange polynomial
    TGraph* approxLag = new TGraph (100, x3, fx3);
    approxLag->SetMarkerStyle(40);
    approxLag->SetMarkerColor(kGreen);
    approxLag->SetLineColor(kGreen);
    approxLag->SetLineWidth(3);


    origin->Draw();
    approxPol->Draw("same");
    approxLag->Draw("same");
    
    return 0;
}