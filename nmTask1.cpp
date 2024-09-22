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

#define NUM_OF_NODS 7
#define LEFT_NOD 0
#define RIGHT_NOD 3
#define NUM_OF_POINTS 100


void basedPoly(TMatrixD &x, TMatrixD &f, double *x1, double *fx2/*const double (&x1)[NUM_OF_POINTS], double (&fx2)[NUM_OF_POINTS]*/){
    //Define the coefficients
    TMatrixD A(NUM_OF_NODS,NUM_OF_NODS);
    TMatrixD aj(NUM_OF_NODS,1);

    for (int i=0; i<NUM_OF_NODS; i++){
        A(i,0) = 1;
        for (int j=1; j<NUM_OF_NODS; j++){
            A(i,j) = pow(x(i,0), j);
        }
    }

    //Calculate the coefficients
    double det = A.Determinant();
    //if (!det) exit;
    aj = A.Invert(&det)*f;
    
    for (int i = 0; i < NUM_OF_POINTS; i++){
        for (int j = 0; j < NUM_OF_NODS; j++){
            fx2[i] += aj(j,0)*pow(x1[i], j);
        }
    }
}

void lagrangePoly(TMatrixD &x, TMatrixD &f,  double *x1, double *fx3/*const double (&x1)[NUM_OF_POINTS], double (&fx3)[NUM_OF_POINTS]*/){
    //Calculate Lagrange coefficients
    double lagrangeCoef[NUM_OF_NODS];
    for (int i = 0; i < NUM_OF_NODS; i++) {
        lagrangeCoef[i] = 1.0;
        for (int j = 0; j < NUM_OF_NODS; j++) {
            if (j != i) {
                lagrangeCoef[i] *= (x(i,0) - x(j,0));
            }
        }
    }
    //Calculate Lagrange polynom
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fx3[i] = 0.0;
        for (int j = 0; j < NUM_OF_NODS; j++){
            double numerator = 1.0;
            for (int k = 0; k < NUM_OF_NODS; k++){
                if (k != j) {
                    numerator *= (x1[i] - x(k,0));
                }
            }
            fx3[i] += (numerator / lagrangeCoef[j]) * f(j,0);
        }
    }
}

int interPoly(){

    //Define the nods of the polynom
    TMatrixD f(NUM_OF_NODS,1);
    TMatrixD x(NUM_OF_NODS,1);


    for (int i=0; i<NUM_OF_NODS; i++){
        x(i,0) = (double)RIGHT_NOD/(double)NUM_OF_NODS * i;
        f(i,0) = sin(x(i,0))*exp(-x(i,0));
    }

    //Calculate divided differences
    double dividedDifferences[NUM_OF_NODS][NUM_OF_NODS];
    for (int i = 0; i < NUM_OF_NODS; i++) {
        dividedDifferences[i][0] = f(i, 0);
    }
    for (int i = 0; i < NUM_OF_NODS; i++){
        for (int j =1; j <= i; j++){
            dividedDifferences[i][j] = (dividedDifferences[i][j-1] - dividedDifferences[i-1][j-1])/(x(i,0)-x(i-j,0));
            std::cout << " F[" << i << "," << j << "] = " << dividedDifferences[i][j];
        }
    }

    std::cout << std::endl;
    std::cout << "=====================================================" << std::endl;


    //Generate original function
    double x1 [NUM_OF_POINTS];
    double fx1 [NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        x1[i] = (double)RIGHT_NOD/(double)NUM_OF_POINTS * i;
        fx1[i] = sin(x1[i])*exp(-x1[i]);
    }

    

    //Generate approximated based Polynomial function 
    double fx2 [NUM_OF_POINTS];
    basedPoly(x, f, x1, fx2);

    //Generate Lagrange polynomial function
    double fx3 [NUM_OF_POINTS];
    lagrangePoly(x, f, x1, fx3);
    
    //Calculate Newton polynomial

    double x4 [NUM_OF_POINTS];
    double fx4 [NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        x4[i] = (double)RIGHT_NOD/(double)NUM_OF_POINTS * i;
        for (int j = 0; j < NUM_OF_NODS; j++){
            double numerator = 1.0;
            for (int k = 0; k < j; k++){
                numerator *= (x4[i] - x(k,0));
            }
            fx4[i] += numerator * dividedDifferences[j][j];  
        }
        std::cout << "f4[" << i << "] = " << fx4[i] << "; x4[" << i << "] = " << x4[i] << std::endl;
    }

    //Calculate difference between original and newton values
    double fdiff[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiff[i] = abs(fx4[i] - fx1[i]);
    }

    //Draw original function
    TGraph* origin = new TGraph (NUM_OF_POINTS, x1, fx1);
    origin->SetMarkerStyle(20);
    origin->SetMarkerColor(kRed);
    origin->SetLineColor(kRed);
    origin->SetLineWidth(14);

    //Draw based approximation
    TGraph* approxPol = new TGraph (NUM_OF_POINTS, x1, fx2);
    approxPol->SetMarkerStyle(30);
    approxPol->SetMarkerColor(kBlue);
    approxPol->SetLineColor(kBlue);
    approxPol->SetLineWidth(10);

    //Draw Lagrange polynomial
    TGraph* approxLag = new TGraph (NUM_OF_POINTS, x1, fx3);
    approxLag->SetMarkerStyle(40);
    approxLag->SetMarkerColor(kGreen);
    approxLag->SetLineColor(kGreen);
    approxLag->SetLineWidth(6);

    //Draw Newton Polynomial
    TGraph* approxNewton = new TGraph (NUM_OF_POINTS, x1, fx4);
    approxNewton->SetMarkerStyle(50);
    approxNewton->SetMarkerColor(kOrange);
    approxNewton->SetLineColor(kOrange);
    approxNewton->SetLineWidth(3);

    //Draw difference between original and newton values
    TGraph* diff = new TGraph (NUM_OF_POINTS, x1, fdiff);
    diff->SetMarkerStyle(60);
    diff->SetMarkerColor(kBlack);
    diff->SetLineColor(kBlack);
    diff->SetLineWidth(2);

    TCanvas* canvas1 = new TCanvas("canvas1", "Graph and approximations", 600, 900);
    canvas1->Divide(1, 2);
    canvas1->cd(1);
    gPad->SetGrid();

    TLegend *legend = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend->SetHeader("Interpolation Functions","C");
    legend->AddEntry(origin, "Original Function", "lp");
    legend->AddEntry(approxPol, "Polynomial approx", "lp");
    legend->AddEntry(approxLag, "Lagrange Polynomial", "lp");
    legend->AddEntry(approxNewton, "Newton Polynomial", "lp");

    origin->Draw();
    approxPol->Draw("same");
    approxLag->Draw("same");
    approxNewton->Draw("same");
    legend->Draw();

    canvas1->cd(2);
    gPad->SetGrid();
    diff->Draw("AL");
    return 0;
}