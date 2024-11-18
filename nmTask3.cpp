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

#define NUM_OF_POINTS 1000
#define LEFT_POINT -3.0
#define RIGHT_POINT 3.0

void leftIntegrator(double *result, double h){
        double  x = LEFT_POINT;
        while (x < RIGHT_POINT){
                *result += exp(-x) * h;
                x += h;
                //std::cout << "Left Integrator: " << result << std::endl;
        }
        std::cout << std::this_thread::get_id() << " is done!" << std::endl;
}

void centralIntegrator(double *result, double h){
        double  x = LEFT_POINT;
        while (x < RIGHT_POINT){
                *result += exp(-(x+0.5*h)) * h;
                x += h;
                //std::cout << "Central Integrator: " << result << std::endl;
        }
        std::cout << std::this_thread::get_id() << " is done!" << std::endl;
}

void trapezoidIntegrator(double *result, double h){
        double  x = LEFT_POINT;
        while (x < RIGHT_POINT){
                *result += 0.5 * (exp(-(x)) + exp(-(x+h))) * h;
                x += h;
                //std::cout << "Trapezoid Integrator: " << result << std::endl;
        }
        std::cout << std::this_thread::get_id() << " is done!" << std::endl;
}

void simpsonIntegrator(double *result, double h){
        double  x = LEFT_POINT;
        while (x < RIGHT_POINT){
                *result += (exp(-(x)) + 4*exp(-(x+0.5*h)) + exp(-(x+h))) / 6.0 * h;
                x += h;
                //std::cout << "Simpson Integrator: " << result << std::endl;
        }
        std::cout << std::this_thread::get_id() << " is done!" << std::endl;
}

void builtInIntegrator(double *result){
        TF1 *f  = new TF1("f", "exp(-x)", LEFT_POINT, RIGHT_POINT);
        *result = f->Integral(LEFT_POINT, RIGHT_POINT);
        std::cout << std::this_thread::get_id() << " is done!" << std::endl;
}

void numDiff(){
    //define original functions
    double  x = 1.0,
            f = exp(-x),
            df  = -exp(-x),
            d2f = exp(-x),
            d3f = -exp(-x),
            d4f = exp(-x);

    double  R1a[NUM_OF_POINTS], 
            R1b[NUM_OF_POINTS], 
            R1c[NUM_OF_POINTS], 
            R2a[NUM_OF_POINTS], 
            R2b[NUM_OF_POINTS],
            R1a1[NUM_OF_POINTS],
            R1b1[NUM_OF_POINTS], 
            R1c1[NUM_OF_POINTS],
            R2a1[NUM_OF_POINTS], 
            R2b1[NUM_OF_POINTS],
            h[NUM_OF_POINTS];
    // calculate numerical derivatives
    for (size_t i = 0; i < NUM_OF_POINTS; i++){
        h[i] = 0.0001 * i + 1e-5;
        double  fleft = exp(-(x-h[i])),
                dfleft = -exp(-(x-h[i])),
                d2fleft = exp(-(x-h[i])),
                d3fleft = -exp(-(x-h[i])),
                fright = exp(-(x+h[i])),
                fright2 = exp(-(x+2*h[i])),
                df1a = (fright - f)/h[i],
                df1b = (fright - fleft)/(2*h[i]),
                df1c = (-3.0*f + 4*fright - fright2)/(2*h[i]),
                df2a = (f - 2*fright + fright2)/(h[i]*h[i]),
                df2b = (fleft - 2*f + fright)/(h[i]*h[i]);
        R1a[i] = 0.5 * d2f * h[i];
        R1a1[i] = abs(df - df1a);
        R1b[i] = -1/6.0 * d3f * h[i]*h[i];
        R1b1[i] = df - df1b;
        R1c[i] = -1.0/3.0 * d3f * h[i]*h[i];
        R1c1[i] = abs(df - df1c);
        R2a[i] = -d3fleft * h[i];
        R2a1[i] = d2f - df2a;
        R2b[i] = 1/12.0 * d4f * h[i]*h[i];
        R2b1[i] = abs(d2f - df2b);
    }

    double hlp[2] = {-0.02, 0.05}; 
    double xhelp[2] = {0, 0.1};
    //Draw the errors
    TCanvas *c1 = new TCanvas("c1", "Numerical Derivatives", 800, 600);
    TGraph *g = new TGraph(2, xhelp, hlp);
    g -> SetLineWidth(0);

    TGraph *gR1a = new TGraph(NUM_OF_POINTS, h, R1a);
    gR1a -> SetLineStyle(5);
    gR1a -> SetLineWidth(6);
    gR1a -> SetLineColor(kOrange);

    TGraph *gR1a1 = new TGraph(NUM_OF_POINTS, h, R1a1);
    gR1a1 -> SetLineStyle(1);
    gR1a1 -> SetLineWidth(6);
    gR1a1 -> SetLineColor(kOrange);

    TGraph *gR1b = new TGraph(NUM_OF_POINTS, h, R1b);
    gR1b -> SetLineStyle(6);
    gR1b -> SetLineWidth(6);
    gR1b -> SetLineColor(kGreen);

    TGraph *gR1b1 = new TGraph(NUM_OF_POINTS, h, R1b1);
    gR1b1 -> SetLineStyle(1);
    gR1b1 -> SetLineWidth(6);
    gR1b1 -> SetLineColor(kGreen);

    TGraph *gR1c = new TGraph(NUM_OF_POINTS, h, R1c);
    gR1c -> SetLineStyle(7);
    gR1c -> SetLineWidth(6);
    gR1c -> SetLineColor(kBlue);

    TGraph *gR1c1 = new TGraph(NUM_OF_POINTS, h, R1c1);
    gR1c1 -> SetLineStyle(1);
    gR1c1 -> SetLineWidth(6);
    gR1c1 -> SetLineColor(kBlue);

    TGraph *gR2a = new TGraph(NUM_OF_POINTS, h, R2a);
    gR2a -> SetLineStyle(8);
    gR2a -> SetLineWidth(6);
    gR2a -> SetLineColor(kRed);

    TGraph *gR2a1 = new TGraph(NUM_OF_POINTS, h, R2a1);
    gR2a1 -> SetLineStyle(1);
    gR2a1 -> SetLineWidth(6);
    gR2a1 -> SetLineColor(kRed);

    TGraph *gR2b = new TGraph(NUM_OF_POINTS, h, R2b);
    gR2b -> SetLineStyle(9);
    gR2b -> SetLineWidth(6);
    gR2b -> SetLineColor(kMagenta);

    TGraph *gR2b1 = new TGraph(NUM_OF_POINTS, h, R2b1);
    gR2b1 -> SetLineStyle(1);
    gR2b1 -> SetLineWidth(6);
    gR2b1 -> SetLineColor(kMagenta);

    c1->SetGrid();
    c1->SetLogy();
    c1->SetLogx();

    TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
    leg->SetHeader("Numerical Derivatives' errors","C");
    leg->AddEntry(gR1a, "R_{1a}", "l");
    leg->AddEntry(gR1b, "R_{1b}", "l");
    leg->AddEntry(gR1c, "R_{1c}", "l");
    leg->AddEntry(gR2a, "R_{2a}", "l");
    leg->AddEntry(gR2b, "R_{2b}", "l");
    leg->SetHeader("Real Numerical Derivatives' errors","C");
    leg->AddEntry(gR1a1, "R_{1a1}", "l");
    leg->AddEntry(gR1b1, "R_{1b1}", "l");
    leg->AddEntry(gR1c1, "R_{1c1}", "l");
    leg->AddEntry(gR2a1, "R_{2a1}", "l");
    leg->AddEntry(gR2b1, "R_{2b1}", "l");


    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gR1a);
    mg->Add(gR1b);
    mg->Add(gR1c);
    mg->Add(gR2a);
    mg->Add(gR2b);

    mg->Add(gR1a1);
    mg->Add(gR1b1);
    mg->Add(gR1c1);
    mg->Add(gR2a1);
    mg->Add(gR2b1);

    mg->SetTitle("Numerical Derivatives' Errors");

    mg->Draw("AC");
    //g->Draw("AC");
    gR1a->Draw("same");
    gR1b->Draw("same");
    gR1c->Draw("same");
    gR2a->Draw("same");
    gR2b->Draw("same");
    leg->Draw();

}

void numInt(){
        // Code for numerical integration goes here
        double  leftInt = 0.0, leftInt2 = 0.0, leftInt4 = 0.0,
                ctrInt = 0.0, ctrInt2 = 0.0, ctrInt4 = 0.0,
                trpInt = 0.0, trpInt2 = 0.0, trpInt4 = 0.0,
                simInt = 0.0, simInt2 = 0.0, simInt4 = 0.0,
                bInIntegral = 0.0,
                realInt = 20.036;
        
        double h = 1e-6;

        std::thread th1(leftIntegrator, &leftInt, h);
        std::thread th2(centralIntegrator, &ctrInt, h);
        std::thread th3(trapezoidIntegrator, &trpInt, h);
        std::thread th4(simpsonIntegrator, &simInt, h);
        std::thread th5(builtInIntegrator, &bInIntegral);
        std::thread th6(leftIntegrator, &leftInt2, h/2.0);
        std::thread th7(centralIntegrator, &ctrInt2, h/2.0);
        std::thread th8(trapezoidIntegrator, &trpInt2, h/2.0);
        std::thread th9(simpsonIntegrator, &simInt2, h/2.0);

        th1.join();
        th2.join();
        th3.join();
        th4.join();
        th5.join();
        th6.join();
        th7.join();
        th8.join();
        th9.join();
        
        std::thread th10(leftIntegrator, &leftInt4, h/4.0);
        std::thread th11(centralIntegrator, &ctrInt4, h/4.0);
        std::thread th12(trapezoidIntegrator, &trpInt4, h/4.0);
        std::thread th13(simpsonIntegrator, &simInt4, h/4.0);

        th10.join();
        th11.join();
        th12.join();
        th13.join();

        double qLeft, qCtr, qTrp, qSim;

        qLeft = log2((leftInt2 - leftInt)/(leftInt4 - leftInt2));
        qCtr = log2((ctrInt2 - ctrInt)/(ctrInt4 - ctrInt2));
        qTrp = log2((trpInt2 - trpInt)/(trpInt4 - trpInt2));
        qSim = log2((simInt2 - simInt)/(simInt4 - simInt2));

        double  leftErrorReal[NUM_OF_POINTS],
                ctrErrorReal[NUM_OF_POINTS],
                trpErrorReal[NUM_OF_POINTS],
                simErrorReal[NUM_OF_POINTS];

        std::cout << "Left Riemann Sum:\t" << leftInt << "\t" << leftInt2 << "\t" << leftInt4 << std::endl;
        std::cout << "Central Riemann Sum:\t" << ctrInt << "\t" << ctrInt2 << "\t" << ctrInt4 << std::endl;
        std::cout << "Trapezoidal Rule:\t" << trpInt << "\t" << trpInt2 << "\t" << trpInt4 << std::endl;
        std::cout << "Simpson's Rule:\t\t" << simInt << "\t" << simInt2 << "\t" << simInt4 << std::endl;
        std::cout << "Built-in Integral:\t" << bInIntegral << std::endl;

        std::cout << "Order of Accuracy:\nLeft: " << qLeft << ", Central: " << qCtr << ", Trapezoidal: " << qTrp << ", Simpson: " << qSim << std::endl;

        double h1[NUM_OF_POINTS];
        
        for(int i = 0; i < NUM_OF_POINTS; i++){
                switch(i){
                        case 0:
                                h1[i] = h;
                                leftErrorReal[i] = abs(leftInt - realInt);
                                ctrErrorReal[i] = abs(ctrInt - realInt);
                                trpErrorReal[i] = abs(trpInt - realInt);
                                simErrorReal[i] = abs(simInt - realInt);
                        default:
                                h1[i] = h+double(1.0/(NUM_OF_POINTS*10)*i);
                                double  leftIntP = 0.0,
                                        ctrIntP = 0.0,
                                        trpIntP = 0.0,
                                        simIntP = 0.0;
                                std::thread th1(leftIntegrator, &leftIntP, h1[i]);
                                std::thread th2(centralIntegrator, &ctrIntP, h1[i]);
                                std::thread th3(trapezoidIntegrator, &trpIntP, h1[i]);
                                std::thread th4(simpsonIntegrator, &simIntP, h1[i]);
                                th1.join();
                                th2.join();
                                th3.join();
                                th4.join();
                                leftErrorReal[i] = abs(leftIntP - realInt);
                                ctrErrorReal[i] = abs(ctrIntP - realInt);
                                trpErrorReal[i] = abs(trpIntP - realInt);
                                simErrorReal[i] = abs(simIntP - realInt);
                }
        }

        TCanvas *c1 = new TCanvas("c1","c1",800,600);

        c1->SetGrid();
        c1->SetLogy();
        c1->SetLogx();

        TGraph *gRleft = new TGraph(NUM_OF_POINTS, h1, leftErrorReal);
        gRleft->SetLineStyle(5);
        gRleft->SetLineWidth(6);
        gRleft->SetLineColor(kRed);
        TGraph *gRctr = new TGraph(NUM_OF_POINTS, h1, ctrErrorReal);
        gRctr->SetLineStyle(7);
        gRctr->SetLineWidth(6);
        gRctr->SetLineColor(kBlue);
        TGraph *gRtrp = new TGraph(NUM_OF_POINTS, h1, trpErrorReal);
        gRtrp->SetLineStyle(3);
        gRtrp->SetLineWidth(6);
        gRtrp->SetLineColor(kGreen);
        TGraph *gRsim = new TGraph(NUM_OF_POINTS, h1, simErrorReal);
        gRsim->SetLineStyle(4);
        gRsim->SetLineWidth(6);
        gRsim->SetLineColor(kOrange);

        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gRleft);
        mg->Add(gRctr);
        mg->Add(gRtrp);
        mg->Add(gRsim);
        mg->SetTitle("Error vs. Step Size");

        TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
        leg->SetHeader("Numerical Integration Errors");
        leg->AddEntry(gRleft,"L");
        leg->AddEntry(gRctr,"L");
        leg->AddEntry(gRtrp,"L");
        leg->AddEntry(gRsim,"L");

        mg->Draw("AL");
        leg->Draw();


}