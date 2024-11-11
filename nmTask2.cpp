#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
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
void fillSpline3Matrix(TMatrixD &spline3, double *h, int num);
void fillSplineB3Matrix(TMatrixD &B3, double *h, int num);
void splineInterpolation(TMatrixD &x, TMatrixD &f, double *x1, double *S, int num);
void spline3Interpolation(TMatrixD &x, TMatrixD &f, double *x1, double *S3, int num);
void splineB1Interpolation(TMatrixD &x, TMatrixD &f, double *x1, double *B1, int num);

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

void fillSpline3Matrix(TMatrixD &spline3, double *h, int num){
    //Fill the spline3 matrix
    //Fill the first row
    spline3(0,0) = 1/h[0];
    spline3(0,1) = -1/h[0] - 1/h[1];
    spline3(0,2) = 1/h[1];
    std::cout << spline3(0,0) << "; " << spline3(0,1) << "; " << spline3(0,2);
    for (int i = 3; i < num; i++){
        spline3(0,i) = 0;
        std::cout << "; " << spline3(0,i);
    }
    std::cout << std::endl;
    //Fill center
    for(int i = 1; i < num-1; i++){
        for (int j = 0; j < num; j++){
            if (j == i-1) spline3(i,j) = h[j];
            else if (j == i) spline3(i,j) = 2 * (h[j-1] + h[j]);
            else if (j == i+1) spline3(i,j) = h[j-1];
            else spline3(i,j) = 0.0;
            std::cout << spline3(i,j);
            if (j<num-1) std::cout << "; ";
        }
        std::cout << std::endl;
    }
    //Fill the last row
    for (int i = 0; i < num-3; i++){
        spline3(num-1,i) = 0;
        std::cout << spline3(num-1,i) << "; ";
    }
    spline3(num-1,num-3) = 1/h[num-3];
    spline3(num-1,num-2) = -1/h[num-2] - 1/h[num-3];
    spline3(num-1,num-1) = 1/h[num-2];
    std::cout << spline3(num-1,num-3) << "; " << spline3(num-1,num-2) << "; " << spline3(num-1,num-1);
    std::cout << std::endl;
}

void fillSplineB3Matrix(TMatrixD &splineB3, TMatrixD &x, double *h, int num){
    //Fill the splineB3 matrix
    //Fill the first row
    splineB3(0,0) = -2/(h[0]*h[0])+ 3/(h[0]*h[0]) * (x(0,0) + 0.5) / h[0];
    //splineB3(0,0) = 12;
    splineB3(0,1) = -2/(h[0]*h[0]);
    splineB3(0,2) = 1/(h[0]*h[0]) * (2 + (x(0,0) - x(1,0))/h[0]);
    //splineB3(0,2) = 12;
    std::cout << splineB3(0,0) << "; " << splineB3(0,1) << "; " << splineB3(0,2);
    for (int i = 3; i < num+2; i++){
        splineB3(0,i) = 0;
        std::cout << "; " << splineB3(0,i);
    }
    std::cout << std::endl;
    //Fill center
    for(int i = 1; i < num+1; i++){
        for (int j = 0; j < num+2; j++){
            if (j == i-1) {
                if (i == 1) splineB3(i,j) = 1.0/6.0 * pow((2 - (x(0,0) + 0.5)/h[0]),3);
                else splineB3(i,j) = 1.0/6.0 * pow((2 - (x(j,0) - x(j-1,0))/h[0]), 3);
                //else splineB3(i,j) = 5.0/30.0;
            }
            else if (j == i) splineB3(i,j) = 2.0/3.0;
            else if (j == i+1) {
                if (i == num) splineB3(i,j) = 1.0/6.0 * pow((2 - (0.5)/h[0]), 3);
                else splineB3(i,j) = 1.0/6.0 * pow((2 + (x(j-2,0) - x(j-1,0))/h[0]), 3);
                //else splineB3(i,j) = -5.0/30.0;
            }
            else splineB3(i,j) = 0.0;
            std::cout << splineB3(i,j) << "; ";
        }
        std::cout << std::endl;     
    }
    // Fill the last row
    for (int i = 0; i < num-1; i++){
        splineB3(num+1,i) = 0.0;
        std::cout << splineB3(num+1,i) << "; ";
    }
    splineB3(num+1,num-1) = -2/(h[0]*h[0])+ 3/(h[0]*h[0]) * (x(num-1,0) - x(num-2,0)) / h[0];
    //splineB3(num+1,num-1) = 12;
    splineB3(num+1,num) = -2/(h[0]*h[0]);
    splineB3(num+1,num+1) = 1/(h[0]*h[0]) * (2 + (x(num-1,0) - 3.5)/h[0]);
    //splineB3(num+1,num+1) = 12;
    
    std::cout << splineB3(num+1,num-1) << "; " << splineB3(num+1,num) << "; " << splineB3(num+1,num+1) << std::endl;

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

void spline3Interpolation(TMatrixD &x, TMatrixD &f, double *x1, double *S3, int num){
    double h[num-1];
    for (int i = 0; i < num-1; i++){
        h[i] = x(i+1,0) - x(i,0);
    }
    //Calculate the spline matrix
    TMatrixD spline3(num, num);
    fillSpline3Matrix(spline3, h, num);

    //calculate the devided differences
    TMatrixD divDiff(num, 1);
    divDiff(0,0) = 0;
    divDiff(num-1,0) = 0;
    for (int i = 1; i < num-1; i++){
        divDiff(i,0) = 3.0 * ((f(i+1,0) - f(i,0)) / h[i] - (f(i,0) - f(i-1,0))/h[i-1]);
    }

    //Calculate the C coefficients
    TMatrixD cVector = spline3.Invert() * divDiff;
    //Calculate the B & D coefficients
    TMatrixD b(num-1,1);
    TMatrixD d(num-1,1);
    for (int i = 0; i < num-1; i++){
        b(i,0) = (f(i+1,0) - f(i,0)) / h[i] - h[i] * (2.0 * cVector(i,0) + cVector(i+1,0)) / 3.0;
        if (i == 1 || i == num-1){
            d[i] = d[i-1];
            continue;
        }
        d(i,0) = (cVector(i+1,0) - cVector(i,0)) / (3.0 * h[i]);
    }
    std::cout << "Calculation of B&&D coefficients completed" << std::endl;
     //Calculate the spline values
    int k = 1;
    for (int i = 0; i < NUM_OF_POINTS; i++){
        if (k < num && x1[i] > x(k,0)){
            k++;
        }
        int j = k-1;
        S3[i] = f(j,0) + b(j,0)*(x1[i] - x(j,0)) + cVector(j,0)*(x1[i] - x(j,0))*(x1[i] - x(j,0)) + d(j,0)*(x1[i] - x(j,0))*(x1[i] - x(j,0))*(x1[i] - x(j,0));
    }
}

void splineB1Interpolation(TMatrixD &x, TMatrixD &f, double *x1, double *B1, int num){
    TMatrixD splineB1(num,num);

    std::cout << "Calculation of B1 coefficients completed" << std::endl;

    
    int k = 0;
    for (int i = 0; i < NUM_OF_POINTS; i++){
        if (x1[i] > x(k+1,0)) k++;
        B1[i] = 2 * (f(k,0) * (x(k+1,0) - x1[i]) + f(k+1,0) * (x1[i] - x(k,0)));
    }
}

void splineB3Interpolation(TMatrixD &x, TMatrixD &f, double *x1, double *B3, int num){
    //Calculate the function differences
    double h[num-1];
    for (int i = 0; i < num-1; i++){
        h[i] = x(i+1,0) - x(i,0);
    }

    TMatrixD splineB3(num+2, num+2);
    TMatrixD f3(num+2, 1);
    f3(0, 0) = 0.0;
    f3(num+1, 0) = 0.0;
    for (int i = 1; i < num+1; i++){
        f3(i,0) = f(i-1, 0);
    }
    fillSplineB3Matrix(splineB3, x, h, num);

    TMatrixD aj = splineB3.Invert() * f3;

    //Calculate the B3 coefficients
    int k = 0;
    int j = 0;
    for (int i = 0; i < NUM_OF_POINTS; i++){
        if (x1[i] > x(k+1,0)) {k++; j++;}
        if (k == 0) B3[i] = aj(j,0) * 1.0/6.0 * pow((2 - (x1[i] + 0.5)/0.5),3) + aj(j+1,0) * (2.0/3.0 - pow(((x1[i] - x(k, 0))/0.5),2) + 0.5*pow(((x1[i] - x(k, 0))/0.5),3)) + aj(j+2,0) * (2.0/3.0 - pow(((x1[i] - x(k+1, 0))/0.5),2) - 0.5*pow(((x1[i] - x(k+1, 0))/0.5),3)) + aj(j+3,0) * 1/6.0 * pow((2 + (x1[i] - x(k+2,0))/0.5),3);
        else if (k == num - 2) B3[i] = aj(j,0) * 1.0/6.0 * pow((2 - (x1[i] - x(k-1,0))/0.5),3) + aj(j+1,0) * (2.0/3.0 - pow(((x1[i] - x(k, 0))/0.5),2) + 0.5*pow(((x1[i] - x(k, 0))/0.5),3)) + aj(j+2,0) * (2.0/3.0 - pow(((x1[i] - x(k+1, 0))/0.5),2) - 0.5*pow(((x1[i] - x(k+1, 0))/0.5),3)) + aj(j+3,0) * 1/6.0 * pow((2 + (x1[i] - 3.5)/0.5),3);
        else B3[i] = aj(j,0) * 1.0/6.0 * pow((2 - (x1[i] - x(k-1,0))/0.5),3) + aj(j+1,0) * (2.0/3.0 - pow(((x1[i] - x(k, 0))/0.5),2) + 0.5*pow(((x1[i] - x(k, 0))/0.5),3)) + aj(j+2,0) * (2.0/3.0 - pow(((x1[i] - x(k+1, 0))/0.5),2) - 0.5*pow(((x1[i] - x(k+1, 0))/0.5),3)) + aj(j+3,0) * 1/6.0 * pow((2 + (x1[i] - x(k+2, 0))/0.5),3);
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

    std::cout << std::endl;
    //Calculate spline3 interpolation
    double S3[NUM_OF_POINTS];
    spline3Interpolation(x, f, x1, S3, NUM_OF_NODS);

    // Calculate B1 spline interpolation
    double B1[NUM_OF_POINTS];
    splineB1Interpolation(x, f, x1, B1, NUM_OF_NODS);

    // Calculate B3 spline interpolation
    double B3[NUM_OF_POINTS];
    splineB3Interpolation(x, f, x1, B3, NUM_OF_NODS);

    //Calculate difference between original and spline B3 values
    double fdiffB3[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffB3[i] = abs(B3[i] - fx1[i]);
    }
    //Calculate difference between original and spline B1 values
    double fdiffB1[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffB1[i] = abs(B1[i] - fx1[i]);
    }

    //Calculate difference between original and spline values
    double fdiffspline[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffspline[i] = abs(S[i] - fx1[i]);
    }

    //Calculate difference between original and spline3 values
    double fdiffspline3[NUM_OF_POINTS];
    for (int i = 0; i < NUM_OF_POINTS; i++){
        fdiffspline3[i] = abs(S3[i] - fx1[i]);
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

    //Draw spline3 interpolation
    TGraph* spline3 = new TGraph (NUM_OF_POINTS, x1, S3);
    spline3->SetMarkerStyle(22);
    spline3->SetMarkerSize(0.1);
    spline3->SetMarkerColor(kBlue);
    spline3->SetLineColor(kBlue);
    spline3->SetLineStyle(2);
    spline3->SetLineWidth(2);

    //Draw B1 spline interpolation
    TGraph* splineB1 = new TGraph (NUM_OF_POINTS, x1, B1);
    splineB1->SetMarkerStyle(23);
    splineB1->SetMarkerSize(0.1);
    splineB1->SetMarkerColor(kGreen);
    splineB1->SetLineColor(kGreen);
    splineB1->SetLineStyle(3);
    splineB1->SetLineWidth(2);

    //Draw B3 spline interpolation
    TGraph* splineB3 = new TGraph (NUM_OF_POINTS, x1, B3);
    splineB3->SetMarkerStyle(24);
    splineB3->SetMarkerSize(0.1);
    splineB3->SetMarkerColor(kCyan);
    splineB3->SetLineColor(kCyan);
    splineB3->SetLineStyle(4);
    splineB3->SetLineWidth(2);

    //Draw difference between original and spline values
    TGraph* diffspline = new TGraph (NUM_OF_POINTS, x1, fdiffspline);
    diffspline->SetMarkerStyle(60);
    diffspline->SetMarkerColor(kOrange);
    diffspline->SetLineColor(kOrange);
    diffspline->SetLineStyle(1);
    diffspline->SetLineWidth(2);

    //Draw difference between original and spline3 values
    TGraph* diffspline3 = new TGraph (NUM_OF_POINTS, x1, fdiffspline3);
    diffspline3->SetMarkerStyle(60);
    diffspline3->SetMarkerColor(kBlue);
    diffspline3->SetLineColor(kBlue);
    diffspline3->SetLineStyle(2);
    diffspline3->SetLineWidth(2);

    //Draw difference between original and spline B1 values
    TGraph* diffB1 = new TGraph (NUM_OF_POINTS, x1, fdiffB1);
    diffB1->SetMarkerStyle(60);
    diffB1->SetMarkerColor(kGreen);
    diffB1->SetLineColor(kGreen);
    diffB1->SetLineStyle(3);
    diffB1->SetLineWidth(2);

    //Draw difference between original and spline B3 values
    TGraph* diffB3 = new TGraph (NUM_OF_POINTS, x1, fdiffB3);
    diffB3->SetMarkerStyle(60);
    diffB3->SetMarkerColor(kCyan);
    diffB3->SetLineColor(kCyan);
    diffB3->SetLineStyle(4);
    diffB3->SetLineWidth(2);

    

    //Draw all

    TCanvas *canvas1 = new TCanvas("canvas4", "Spline interpolation", 1200, 900);
    canvas1->Divide(1, 2);
    canvas1->cd(1);
    gPad->SetGrid();

    TLegend* leg1 = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg1->SetFillColor(0);
    leg1->AddEntry(origin, "Original function", "L");
    leg1->AddEntry(spline, "Spline interpolation", "L");
    leg1->AddEntry(spline3, "Spline3 interpolation", "L");
    leg1->AddEntry(splineB1, "B1 spline interpolation", "L");
    leg1->AddEntry(splineB3, "B3 spline interpolation", "L");


    origin->Draw();
    nodes7->Draw("Psame");
    spline->Draw("same");
    spline3->Draw("same");
    splineB1->Draw("same");
    splineB3->Draw("same");

    leg1->Draw("same");

    canvas1->cd(2);
    gPad->SetGrid();

    TLegend* leg2 = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg2->SetFillColor(0);
    leg2->AddEntry(diffspline, "Spline Error", "L");
    leg2->AddEntry(diffspline3, "Spline3 Error", "L");
    leg2->AddEntry(diffB1, "B1 spline Error", "L");
    leg2->AddEntry(diffB3, "B3 spline Error", "L");

    diffspline->Draw("AL");
    diffspline3->Draw("same");
    diffB1->Draw("same");
    diffB3->Draw("same");

    leg2->Draw("same");

    return 0;
}