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

void fillSplineMatrix(TMatrixD &spline, double *h, int num);
void splineInterpolation(TMatrixD &x, TMatrixD &f, double *x1, double *S, int num);

void fillSplineMatrix(TMatrixD &spline, double *h, int num){
    //Fill the spline matrix
    //Fill the first row
    spline(0,0) = 1.0;
    std::cout << spline(0,0);
    for (int i = 1; i < num; i++){
        spline(0,i) = 0;
        std::cout << "; " << spline(0,i);
    }
    std::cout << std::endl;
    //Fill center
    for (int i=1; i<num-1; i++){
        for (int j=0; j<num; j++){
            if (j==i-1) spline(i,j) = h[i];
            else if (j==i) spline(i,j) = 2 * (h[i] + h[i+1]);
            else if (j==i+1) spline(i,j) = h[i+1];
            else spline(i,j) = 0.0;
            std::cout << spline(i,j);
            if (j<num-1) std::cout << "; ";
        }
        std::cout << std::endl;
    }
    //Fill the last row
    for (int i=0; i<num-1; i++){
        spline(num-1,i) = 0;
        std::cout << spline(num-1,i);
        if (i<num-1) std::cout << "; ";
    }
    spline(num-1,num-1) = 1.0;
    std::cout << spline(num-1,num-1) << std::endl;
}
void splineInterpolation(TMatrixD &x, TMatrixD &f, double *x1, double *S, int num){

    //Calculate the function differences
    double h[num-1];
    for (int i = 0; i < num-1; i++){
        h[i] = x(i+1,0) - x(i,0);
    }

    //Calculate the spline matrix
    TMatrixD spline(num, num);
    fillSplineMatrix(spline, h, num);

    //Calculating devided differences
    TMatrixD divDiff(num, 1);
    divDiff(0,0) = 0;
    divDiff(num-1,0) = 0;
    for (int i = 1; i < num-1; i++){
        divDiff(i,0) = 3.0 * ((f(i+1,0) - f(i,0)) / h[i] - (f(i,0) - f(i-1,0))/h[i-1]);
    }
    //Calculate the C coefficients
    TMatrixD cVector = spline.Invert() * divDiff;
    //Calculate the B & D coefficients
    TMatrixD b(num-1,1);
    TMatrixD d(num-1,1);
    for (int i = 0; i < num-1; i++){
        b(i,0) = (f(i+1,0) - f(i,0)) / h[i] - h[i] * (2.0 * cVector(i,0) + cVector(i+1,0)) / 3.0;
        d(i,0) = (cVector(i+1,0) - cVector(i,0)) / (3.0 * h[i]);
    }

    //Calculate the spline values
    int k = 1;
    for (int i = 0; i < NUM_OF_POINTS; i++){
        if (k < NUM_OF_NODS && x1[i] > x(k,0)){
            k++;
        }
        int j = k-1;
        S[i] = f(j,0) + b(j,0)*(x1[i] - x(j,0)) + cVector(j,0)*(x1[i] - x(j,0))*(x1[i] - x(j,0)) + d(j,0)*(x1[i] - x(j,0))*(x1[i] - x(j,0))*(x1[i] - x(j,0));
    }
}

int interPoly(){

    //Define the nodes of the polynom

    //Number of nodes is 7
    TMatrixD f(NUM_OF_NODS,1);
    TMatrixD x(NUM_OF_NODS,1);


    for (int i=0; i<NUM_OF_NODS; i++){
        //x(i,0) = (double)RIGHT_NOD/(double)NUM_OF_NODS * i;
        x(i, 0) = i * 0.5;
        f(i,0) = sin(x(i,0))*exp(-x(i,0));
    }

    //Generate original function
    double x1 [NUM_OF_POINTS];
    double fx1 [NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        x1[i] = (double)RIGHT_NOD/(double)NUM_OF_POINTS * i;
        fx1[i] = sin(x1[i])*exp(-x1[i]);
    }

    //Calculate Spline interpolation

    double S[NUM_OF_POINTS];
    splineInterpolation(x, f, x1, S, NUM_OF_NODS);

    //Calculate difference between original and spline values
    double fdiffspline[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffspline[i] = abs(S[i] - fx1[i]);
    }

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

    //Draw original function
    TGraph* origin = new TGraph (NUM_OF_POINTS, x1, fx1);
    origin->SetMarkerStyle(20);
    origin->SetMarkerSize(0.1);
    origin->SetMarkerColor(kRed);
    origin->SetLineColor(kRed);
    origin->SetLineStyle(10);
    origin->SetLineWidth(2);
    
    //Draw spline interpolation
    TGraph* spline = new TGraph (NUM_OF_POINTS, x1, S);
    spline->SetMarkerStyle(30);
    spline->SetMarkerSize(0.1);
    spline->SetMarkerColor(kOrange);
    spline->SetLineColor(kOrange);
    spline->SetLineStyle(1);
    spline->SetLineWidth(2);

    //Draw difference between original and spline values
    TGraph* diffspline = new TGraph (NUM_OF_POINTS, x1, fdiffspline);
    diffspline->SetMarkerStyle(60);
    diffspline->SetMarkerColor(kOrange);
    diffspline->SetLineColor(kOrange);
    diffspline->SetLineStyle(1);
    diffspline->SetLineWidth(2);

    TCanvas *canvas4 = new TCanvas("canvas4", "Spline interpolation", 1200, 900);
    canvas4->Divide(1, 2);
    canvas4->cd(1);
    gPad->SetGrid();
    origin->Draw();
    nodes7->Draw("Psame");
    spline->Draw("same");

    canvas4->cd(2);
    gPad->SetGrid();
    diffspline->Draw("AL");

    return 0;
}