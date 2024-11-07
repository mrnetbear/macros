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
#include "TMultiGraph.h"

#define NUM_OF_POINTS 1000

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
        h[i] = 0.0001 * i + 1e-4;
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
    //leg->AddEntry(gR2a, "R_{2a}", "l");
    //leg->AddEntry(gR2b, "R_{2b}", "l");
    leg->SetHeader("Real Numerical Derivatives' errors","C");
    leg->AddEntry(gR1a1, "R_{1a1}", "l");
    leg->AddEntry(gR1b1, "R_{1b1}", "l");
    leg->AddEntry(gR1c1, "R_{1c1}", "l");
    //leg->AddEntry(gR2a1, "R_{2a1}", "l");
    //leg->AddEntry(gR2b1, "R_{2b1}", "l");


    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gR1a);
    mg->Add(gR1b);
    mg->Add(gR1c);
    //mg->Add(gR2a);
    //mg->Add(gR2b);

    mg->Add(gR1a1);
    mg->Add(gR1b1);
    mg->Add(gR1c1);
    //mg->Add(gR2a1);
    //mg->Add(gR2b1);

    mg->SetTitle("Numerical Derivatives' Errors");

    mg->Draw("AC");
    /*g->Draw("AC");
    gR1a->Draw("same");
    gR1b->Draw("same");
    gR1c->Draw("same");
    gR2a->Draw("same");
    gR2b->Draw("same");*/
    leg->Draw();

    TCanvas *c2 = new TCanvas("c2", "Numerical Derivatives in segond grade Errors", 800, 600);
    c2->SetGrid();
    c2->SetLogy();
    c2->SetLogx();
    c2->cd();

    TLegend *leg2 = new TLegend(0.1, 0.7, 0.48, 0.9);
    leg2->SetHeader("Numerical Derivatives' errors in second grade","C");
    leg2->AddEntry(gR2a, "R_{2a}", "l");
    leg2->AddEntry(gR2b, "R_{2b}", "l");
    leg2->AddEntry(gR2a1, "R_{2a1}", "l");
    leg2->AddEntry(gR2b1, "R_{2b1}", "l");

    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(gR2a);
    mg2->Add(gR2b);

    mg2->Add(gR2a1);
    mg2->Add(gR2b1);

    mg2->SetTitle("Numerical Derivatives' Errors in second grade");
    mg2->Draw("AC");
    leg2->Draw();
}