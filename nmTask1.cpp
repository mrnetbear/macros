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
//#define M_PI 3.1415926535


int factorial(int n){
    int res = 1;
    for (int i=1; i<n; i++){
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
            //std::cout << "aj[" << j << "] = " <<aj(j,0) << ";" << std::endl;
        }
        //std::cout << "f[" << i << "] = " << fx2[i] << ";" << std::endl;
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
            //std::cout << " F[" << i << "," << j << "] = " << dividedDifferences[i][j];
        }
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
        //std::cout << "f4[" << i << "] = " << fx4[i] << "; x4[" << i << "] = " << x1[i] << std::endl;
    }
}

void errorEvaluation (TMatrixD &x, TMatrixD &f, double *x1, double *R, int num){
    //Calculate coefficient
    double coefR = 16.0 * sin(M_PI/4.0)*exp(-M_PI/4.0) / factorial(num);
    
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
        //std::cout << "f4[" << i << "] = " << fx4[i] << "; x4[" << i << "] = " << x1[i] << std::endl;
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
    double fx4 [NUM_OF_POINTS];
    newtonPoly(x, f, x1, fx4, NUM_OF_NODS);
    

    //Calculate difference between original and newton values
    double fdiff[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiff[i] = abs(fx4[i] - fx1[i]);
    }

    TMatrixD fmore(NUM_OF_NODS_ACC, 1);
    TMatrixD xmore(NUM_OF_NODS_ACC, 1);


    for (int i=0; i<NUM_OF_NODS_ACC; i++){
        xmore(i,0) = (double)RIGHT_NOD/(double)NUM_OF_NODS_ACC * i;
        fmore(i,0) = sin(xmore(i,0))*exp(-xmore(i,0));
        //std::cout << "xmore[" << i << "] = " << x(i,0) << std::endl;
        //std::cout << "fmore[" << i << "] = " << fmore(i,0) << std::endl;
    }

    double fx5 [NUM_OF_POINTS];
    newtonPoly(xmore, fmore, x1, fx5, NUM_OF_NODS_ACC);

    double fdiffmore[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffmore[i] = abs(fx5[i] - fx1[i]);
    }

    //Calculate R
    double R7[NUM_OF_POINTS];
    double R8[NUM_OF_POINTS];
    errorEvaluation(x, f, x1, R7, NUM_OF_NODS);
    errorEvaluation(xmore, fmore, x1, R8, NUM_OF_NODS_ACC);

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
    approxNewton->SetMarkerColor(kBlack);
    approxNewton->SetLineColor(kBlack);
    approxNewton->SetLineWidth(3);

    //Draw Newton Polynomial
    TGraph* approxNewtonMore = new TGraph (NUM_OF_POINTS, x1, fx5);
    approxNewtonMore->SetMarkerStyle(50);
    approxNewtonMore->SetMarkerColor(kGray);
    approxNewtonMore->SetLineColor(kGray);
    approxNewtonMore->SetLineWidth(1);

    //Draw difference between original and newton values
    TGraph* diff = new TGraph (NUM_OF_POINTS, x1, fdiff);
    diff->SetMarkerStyle(60);
    diff->SetMarkerColor(kBlack);
    diff->SetLineColor(kBlack);
    diff->SetLineWidth(2);

    //Draw difference between original and more newton values
    TGraph* diffmore = new TGraph (NUM_OF_POINTS, x1, fdiffmore);
    diffmore->SetMarkerStyle(60);
    diffmore->SetMarkerColor(kGreen);
    diffmore->SetLineColor(kGreen);
    diffmore->SetLineWidth(2);

    //Draw R7
    TGraph* error7 = new TGraph (NUM_OF_POINTS, x1, R7);
    error7->SetMarkerStyle(60);
    error7->SetMarkerColor(kRed);
    error7->SetLineColor(kRed);
    error7->SetLineWidth(2);

    //Draw R8
    TGraph* error = new TGraph (NUM_OF_POINTS, x1, R8);
    error->SetMarkerStyle(60);
    error->SetMarkerColor(kBlue);
    error->SetLineColor(kBlue);
    error->SetLineWidth(2);

    TCanvas* canvas1 = new TCanvas("canvas1", "Graph and approximations", 900, 600);
    canvas1->SetGrid();

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
    
    TCanvas* canvas2 = new TCanvas("canvas2", "Newton approximations & errors", 1200, 900);
    canvas2->Divide(1, 2);
    canvas2->cd(1);
    gPad->SetGrid();
    TLegend *legend3 = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend3->SetHeader("Newtons Interpolations","C");
    legend3->AddEntry(origin, "original", "lp");
    legend3->AddEntry(approxNewton, "Newton Polynomial n = 7", "lp");
    legend3->AddEntry(approxNewtonMore, "Newton Polynomial n = 8", "lp");

    origin->Draw();
    approxNewton->Draw("same");
    approxNewtonMore->Draw("same");
    legend3->Draw();


    canvas2->cd(2);
    TLegend *legend2 = new TLegend(0.2, 0.2, 0.5, 0.5);
    legend2->SetHeader("Method error","C");
    legend2->AddEntry(diff, "Difference n = 7", "lp");
    legend2->AddEntry(diffmore, "Difference n = 8", "lp");
    legend2->AddEntry(error7, "Evaluated error n = 7" , "lp");
    legend2->AddEntry(error, "Evaluated error n = 8" , "lp");
    gPad->SetGrid();

    error7->Draw("AL");
    diff->Draw("same");
    diffmore->Draw("same");
    error->Draw("same");
    legend2->Draw();
    return 0;
}