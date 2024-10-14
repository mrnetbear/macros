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
#define NUM_OF_NODS_ACC 8
#define LEFT_NOD 0
#define RIGHT_NOD 3
#define NUM_OF_POINTS 100


int factorial(int n){
    int res = 1;
    for (int i=1; i<=n; i++){
        res *= i;
    }
    return res;
}

void basedPoly(TMatrixD &x, TMatrixD &f, double *x1, double *fx2){
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
        fx2[i] = 0.0;
        for (int j = 0; j < NUM_OF_NODS; j++){
            fx2[i] += aj(j,0)*pow(x1[i], j);
        }
    }
}

void lagrangePoly(TMatrixD &x, TMatrixD &f,  double *x1, double *fx3){
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

void newtonPoly(TMatrixD &x, TMatrixD &f,  double *x1, double *fx4, int num){
    //Calculate divided differences
    double dividedDifferences[num][num];
    for (int i = 0; i < num; i++) {
        dividedDifferences[i][0] = f(i, 0);
    }
    for (int i = 0; i < num; i++){
        for (int j =1; j <= i; j++){
            dividedDifferences[i][j] = (dividedDifferences[i][j-1] - dividedDifferences[i-1][j-1])/(x(i,0)-x(i-j,0));
        }
        std::cout << dividedDifferences[i][i]  << std::endl;
    }

    //Calculate Newton polynomial
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fx4[i] = 0.0;
        for (int j = 0; j < num; j++){
            double numerator = 1.0;
            for (int k = 0; k < j; k++){
                numerator *= (x1[i] - x(k,0));
            }
            fx4[i] += numerator * dividedDifferences[j][j];  
        }
    }
}

void errorEvaluation (TMatrixD &x, TMatrixD &f, double *x1, double *R, int num){
    //Calculate coefficient
    double coefR = 16.0 * sin(M_PI/4.0)*exp(-M_PI/4.0) / factorial(num+1);
    
    //Calculate R
    for (int i = 0; i < NUM_OF_POINTS; i++){
        R[i] = 0.0;
        for (int j = 0; j < num; j++){
            double numerator = 1.0;
            for (int k = 0; k < num; k++){
                numerator *= (x1[i] - x(k,0));
            }
            R[i] += std::abs(numerator) * coefR;  
        }
    }

}

void splineInterpolation(TMatrixD &x, TMatrixD &f, double *x1, double *S, int num){

    //Calculate the function differences
    double h[num-1];
    for (int i = 0; i < num-1; i++){
        h[i] = f(i+1,0) - f(i,0);
    }
    TMatrixD spline(num, num);
}

void fillSplineMatrix(TMatrixD &spline, double *h, int num){
    //Fill the spline matrix
    spline(0,0) = 1.0;
    for (int i = i; i < num; i++)
        spline(0,i) = 0;
    for (int i=1; i<num-1; i++){
        for (int j=0; j<num; j++){
            if (j==i-1) spline(i,j) = h[i];
            else if (j==i) spline(i,j) = 2 * (h[i] + h[i+1]);
            else if (j==i+1) spline(i,j) = h[i+1];
            else spline(i,j) = 0.0;
        }
    }
    spline(num-1,num-1) = 1.0;
}

int interPoly(){

    //Define the nodes of the polynom

    //Number of nodes is 7
    TMatrixD f(NUM_OF_NODS,1);
    TMatrixD x(NUM_OF_NODS,1);


    for (int i=0; i<NUM_OF_NODS; i++){
        x(i,0) = (double)RIGHT_NOD/(double)NUM_OF_NODS * i;
        f(i,0) = sin(x(i,0))*exp(-x(i,0));
    }

    //Additional nodes for Newton polynomial
    TMatrixD fmore(NUM_OF_NODS_ACC, 1);
    TMatrixD xmore(NUM_OF_NODS_ACC, 1);


    for (int i=0; i<NUM_OF_NODS_ACC; i++){
        xmore(i,0) = (double)RIGHT_NOD/(double)NUM_OF_NODS_ACC * i;
        fmore(i,0) = sin(xmore(i,0))*exp(-xmore(i,0));
    }

    // Additional nodes for Chebyshev polynom

    TMatrixD fchebyshev(NUM_OF_NODS_ACC, 1);
    TMatrixD xchebyshev(NUM_OF_NODS_ACC, 1);

    for (int i=0; i<NUM_OF_NODS_ACC; i++){
        xchebyshev(i,0) = 1.5 * cos(M_PI*(2*i+1)/(2*NUM_OF_NODS_ACC)) + 1.5;
        fchebyshev(i,0) = sin(xchebyshev(i,0))*exp(-xchebyshev(i,0));
    }

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

    //For 7 nodes
    double fx4 [NUM_OF_POINTS];
    newtonPoly(x, f, x1, fx4, NUM_OF_NODS);

    //For 8 nodes
    double fx5 [NUM_OF_POINTS];
    newtonPoly(xmore, fmore, x1, fx5, NUM_OF_NODS_ACC);

    //For 8 nodes with Chebyshev distribution
    double fx6 [NUM_OF_POINTS];
    newtonPoly(xchebyshev, fchebyshev, x1, fx6, NUM_OF_NODS_ACC);

    //Calculate difference between original and newton values with 7 nodes
    double fdiff[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiff[i] = abs(fx4[i] - fx1[i]);
    }

    //Calculate difference between original and newton values with 8 nodes

    //For 8 nodes with based distribution
    double fdiffmore[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffmore[i] = abs(fx5[i] - fx1[i]);
    }

    //For 8 nodes with Chebyshev distribution
    double fdiffchebyshev[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffchebyshev[i] = abs(fx6[i] - fx1[i]);
    }

    //Calculate R
    double R7[NUM_OF_POINTS];
    double R8[NUM_OF_POINTS];
    double R9[NUM_OF_POINTS];
    errorEvaluation(x, f, x1, R7, NUM_OF_NODS);
    errorEvaluation(xmore, fmore, x1, R8, NUM_OF_NODS_ACC);
    errorEvaluation(xchebyshev, fchebyshev, x1, R9, NUM_OF_NODS_ACC);

    //Draw nodes

    //Draw 7 nodes
    double fnodes7[NUM_OF_NODS], xnodes7[NUM_OF_NODS];
    for (int i = 0; i < NUM_OF_NODS; i++)
    {
        xnodes7[i] = x(i, 0);
        fnodes7[i] = f(i, 0);
    }
    
    TGraph* nodes7 = new TGraph (NUM_OF_NODS, xnodes7, fnodes7);
    nodes7->SetMarkerStyle(21);
    nodes7->SetMarkerColor(7);

    //Draw 8 nodes

    //Based distribution
    double fnodes8[NUM_OF_NODS_ACC], xnodes8[NUM_OF_NODS_ACC];
    for (int i = 0; i < NUM_OF_NODS_ACC; i++)
    {
        xnodes8[i] = xmore(i, 0);
        fnodes8[i] = fmore(i, 0);
    }
    
    TGraph* nodes8 = new TGraph (NUM_OF_NODS_ACC, xnodes8, fnodes8);
    nodes8->SetMarkerStyle(22);
    nodes8->SetMarkerColor(8);

    //Chebyshev distribution
    double fnodes9[NUM_OF_NODS_ACC], xnodes9[NUM_OF_NODS_ACC];
    for (int i = 0; i < NUM_OF_NODS_ACC; i++)
    {
        xnodes9[i] = xchebyshev(i, 0);
        fnodes9[i] = fchebyshev(i, 0);
    }
    
    TGraph* nodes9 = new TGraph (NUM_OF_NODS_ACC, xnodes9, fnodes9);
    nodes9->SetMarkerStyle(23);
    nodes9->SetMarkerColor(9);

    //Draw original function
    TGraph* origin = new TGraph (NUM_OF_POINTS, x1, fx1);
    origin->SetMarkerStyle(20);
    origin->SetMarkerSize(0.1);
    origin->SetMarkerColor(kRed);
    origin->SetLineColor(kRed);
    origin->SetLineStyle(10);
    origin->SetLineWidth(2);

    //Draw based approximation
    TGraph* approxPol = new TGraph (NUM_OF_POINTS, x1, fx2);
    approxPol->SetMarkerStyle(30);
    approxPol->SetMarkerSize(0.1);
    approxPol->SetMarkerColor(kBlue);
    approxPol->SetLineColor(kBlue);
    approxPol->SetLineStyle(8);
    approxPol->SetLineWidth(2);

    //Draw Lagrange polynomial
    TGraph* approxLag = new TGraph (NUM_OF_POINTS, x1, fx3);
    approxLag->SetMarkerStyle(40);
    approxLag->SetMarkerSize(0.1);
    approxLag->SetMarkerColor(kGreen);
    approxLag->SetLineColor(kGreen);
    approxLag->SetLineStyle(7);
    approxLag->SetLineWidth(2);

    //Draw Newton Polynomial with 7 nodes
    TGraph* approxNewton = new TGraph (NUM_OF_POINTS, x1, fx4);
    approxNewton->SetMarkerStyle(50);
    approxNewton->SetMarkerSize(0.1);
    approxNewton->SetMarkerColor(kBlack);
    approxNewton->SetLineColor(kBlack);
    approxNewton->SetLineStyle(9);
    approxNewton->SetLineWidth(2);

    //Draw Newton Polynomial with 8 nodes
    TGraph* approxNewtonMore = new TGraph (NUM_OF_POINTS, x1, fx5);
    approxNewtonMore->SetMarkerStyle(50);
    approxNewtonMore->SetMarkerSize(0.1);
    approxNewtonMore->SetMarkerColor(kGray);
    approxNewtonMore->SetLineColor(kGray);
    approxNewtonMore->SetLineStyle(6);
    approxNewtonMore->SetLineWidth(2);

    //Draw Newton Polynomial with 8 nodes with Chebyshev distribution
    TGraph* approxNewtonChebyshev = new TGraph (NUM_OF_POINTS, x1, fx6);
    approxNewtonChebyshev->SetMarkerStyle(50);
    approxNewtonChebyshev->SetMarkerSize(0.1);
    approxNewtonChebyshev->SetMarkerColor(kMagenta);
    approxNewtonChebyshev->SetLineColor(kMagenta);
    approxNewtonChebyshev->SetLineStyle(5);
    approxNewtonChebyshev->SetLineWidth(2);

    //Draw difference between original and newton values
    TGraph* diff = new TGraph (NUM_OF_POINTS, x1, fdiff);
    diff->SetMarkerStyle(60);
    diff->SetMarkerColor(kBlack);
    diff->SetLineColor(kBlack);
    diff->SetLineStyle(9);
    diff->SetLineWidth(2);

    //Draw difference between original and more newton values
    TGraph* diffmore = new TGraph (NUM_OF_POINTS, x1, fdiffmore);
    diffmore->SetMarkerStyle(60);
    diffmore->SetMarkerColor(kGray);
    diffmore->SetLineColor(kGray);
    diffmore->SetLineStyle(9);
    diffmore->SetLineWidth(2);

    //Draw difference between original and Chebyshev newton values
    TGraph* diffchebyshev = new TGraph (NUM_OF_POINTS, x1, fdiffchebyshev);
    diffchebyshev->SetMarkerStyle(60);
    diffchebyshev->SetMarkerColor(kMagenta);
    diffchebyshev->SetLineColor(kMagenta);
    diffchebyshev->SetLineStyle(9);
    diffchebyshev->SetLineWidth(2);

    //Draw R7
    TGraph* error7 = new TGraph (NUM_OF_POINTS, x1, R7);
    error7->SetMarkerStyle(60);
    error7->SetMarkerColor(kBlack);
    error7->SetLineColor(kBlack);
    error7->SetLineWidth(2);

    //Draw R8
    TGraph* error = new TGraph (NUM_OF_POINTS, x1, R8);
    error->SetMarkerStyle(60);
    error->SetMarkerColor(kGray);
    error->SetLineColor(kGray);
    error->SetLineWidth(2);

    //Draw R9
    TGraph* errorChebyshev = new TGraph (NUM_OF_POINTS, x1, R9);
    errorChebyshev->SetMarkerStyle(60);
    errorChebyshev->SetMarkerColor(kMagenta);
    errorChebyshev->SetLineColor(kMagenta);
    errorChebyshev->SetLineWidth(2);

    TCanvas* canvas1 = new TCanvas("canvas1", "Graph and approximations", 900, 600);
    canvas1->SetGrid();

    TLegend *legend = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend->SetHeader("Interpolation Functions","C");
    legend->AddEntry(origin, "Original Function", "lp");
    legend->AddEntry(approxPol, "Polynomial approx", "lp");
    legend->AddEntry(approxLag, "Lagrange Polynomial", "lp");
    legend->AddEntry(approxNewton, "Newton Polynomial", "lp");

    origin->Draw();
    nodes7->Draw("Psame");
    origin->SetTitle("Interpolation Functions");
    approxPol->Draw("same");
    approxLag->Draw("same");
    approxNewton->Draw("same");

    legend->Draw();
    
    //Draw two Newton polynomial interpolations with errors
    TCanvas* canvas2 = new TCanvas("canvas2", "Newton approximations & errors", 1200, 900);
    canvas2->Divide(1, 2);

    //Newton interpolations
    canvas2->cd(1);
    gPad->SetGrid();
    TLegend *legend3 = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend3->SetHeader("Newtons Interpolations","C");
    legend3->AddEntry(origin, "original", "lp");
    legend3->AddEntry(approxNewton, "Newton Polynomial n = 7", "lp");
    legend3->AddEntry(approxNewtonMore, "Newton Polynomial n = 8", "lp");

    origin->Draw();
    nodes7->Draw("Psame");
    nodes8->Draw("Psame");
    approxNewton->Draw("same");
    approxNewtonMore->Draw("same");
    legend3->Draw();

    //Errors with evaluation
    canvas2->cd(2);
    TLegend *legend2 = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend2->SetHeader("Method error","C");
    legend2->AddEntry(diff, "Difference n = 7", "lp");
    legend2->AddEntry(diffmore, "Difference n = 8", "lp");
    legend2->AddEntry(error7, "Evaluated error n = 7" , "lp");
    legend2->AddEntry(error, "Evaluated error n = 8" , "lp");
    gPad->SetGrid();

    error7->GetYaxis()->SetRangeUser(0.0, 0.0003);
    error7->SetTitle("Errors");
    error7->Draw("AL");
    diff->Draw("same");
    diffmore->Draw("same");
    error->Draw("same");
    legend2->Draw();

    //Draw Chebyshev distribution
    TCanvas *canvas3 = new TCanvas("canvas3", "Newton approximations & errors", 1200, 900);
    canvas3->Divide(1, 2);
    
    //Newton interpolations
    canvas3->cd(1);
    gPad->SetGrid();
    TLegend *legend4 = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend4->SetHeader("Newton Polynomial with Chebyshev distribution","C");
    legend4->AddEntry(origin, "original", "lp");
    legend4->AddEntry(approxNewtonMore, "Newton Polynomial n = 8", "lp");
    legend4->AddEntry(approxNewtonChebyshev, "Newton Polynomial n = 8 with Chebyshev distribution", "lp");

    origin->Draw();
    nodes8->Draw("Psame");
    nodes9->Draw("Psame");
    approxNewtonMore->Draw("same");
    approxNewtonChebyshev->Draw("same");
    legend4->Draw();
    
    //Errors with evaluation
    canvas3->cd(2);
    gPad->SetGrid();
    TLegend *legend5 = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend5->SetHeader("Method error with Chebyshev distribution","C");
    legend5->AddEntry(diffmore, "Difference n = 8", "lp");
    legend5->AddEntry(diffchebyshev, "Difference n = 8 with Chebyshev distribution", "lp");
    legend5->AddEntry(error, "Evaluated error n = 8", "lp");
    legend5->AddEntry(errorChebyshev, "Evaluated error n = 8 with Chebyshev distribution", "lp");

    error->GetYaxis()->SetRangeUser(0.0, 0.00005);
    error->SetTitle("Errors with Chebyshev distribution");
    error->Draw("AL");
    diffmore->Draw("same");
    diffchebyshev->Draw("same");
    errorChebyshev->Draw("same");
    legend5->Draw();

    return 0;
}