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

#define NUM_OF_POINTS 1000

void numDiff(){
    //define original functions
    double  x = 1.0,
            f = exp(-x),
            df  = -exp(-x),
            d2f = exp(-x),
            d3f = -exp(-x),
            d4f = exp(-x);

    double R1a[NUM_OF_POINTS], R1b[NUM_OF_POINTS], R1c[NUM_OF_POINTS], R2a[NUM_OF_POINTS], R2b[NUM_OF_POINTS], h[NUM_OF_POINTS];
    // calculate numerical derivatives
    for (size_t i = 0; i < NUM_OF_POINTS; i++){
        h[i] = 0.0001 * i;
        double  fleft = exp(-(x-h[i])),
                dfleft = -exp(-(x-h[i])),
                d2fleft = exp(-(x-h[i])),
                d3fleft = -exp(-(x-h[i])),
                fright = exp(-(x+h[i])),
                df1a = (fright - f)/h[i],
                df1b = (fright - fleft)/(2*h[i]),
                df1c = (-3.0*fleft + 4*f - fright)/(2*h[i]),
                df2a = (fleft - 2*f + fright)/(h[i]*h[i]),
                df2b = (fleft - 2*f + fright)/(h[i]*h[i]);
        R1a[i] = -0.5 * d2f * h[i];
        R1b[i] = -1/6.0 * d3f * h[i]*h[i];
        R1c[i] = -2.0/3.0 * d3fleft * h[i]*h[i];
        R2a[i] = -d3fleft * h[i];
        R2b[i] = -1/12.0 * d4f * h[i]*h[i];
    }

    double hlp[2] = {-0.02, 0.05}; 
    double xhelp[2] = {0, 0.1};
    //Draw the errors
    TCanvas *c1 = new TCanvas("c1", "Numerical Derivatives", 800, 600);
    TGraph *g = new TGraph(2, xhelp, hlp);
    g -> SetLineWidth(0);
    TGraph *gR1a = new TGraph(NUM_OF_POINTS, h, R1a);
    gR1a -> SetLineStyle(5);
    gR1a -> SetLineColor(kOrange);
    TGraph *gR1b = new TGraph(NUM_OF_POINTS, h, R1b);
    gR1b -> SetLineStyle(6);
    gR1b -> SetLineColor(kGreen);
    TGraph *gR1c = new TGraph(NUM_OF_POINTS, h, R1c);
    gR1c -> SetLineStyle(7);
    gR1c -> SetLineColor(kBlue);
    TGraph *gR2a = new TGraph(NUM_OF_POINTS, h, R2a);
    gR2a -> SetLineStyle(8);
    gR2a -> SetLineColor(kRed);
    TGraph *gR2b = new TGraph(NUM_OF_POINTS, h, R2b);
    gR2b -> SetLineStyle(9);
    gR2b -> SetLineColor(kMagenta);

    c1->SetGrid();

    TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
    leg->SetHeader("Numerical Derivatives' errors","C");
    leg->AddEntry(gR1a, "R_{1a}", "l");
    leg->AddEntry(gR1b, "R_{1b}", "l");
    leg->AddEntry(gR1c, "R_{1c}", "l");
    leg->AddEntry(gR2a, "R_{2a}", "l");
    leg->AddEntry(gR2b, "R_{2b}", "l");

    g->Draw("AC");
    gR1a->Draw("same");
    gR1b->Draw("same");
    gR1c->Draw("same");
    gR2a->Draw("same");
    gR2b->Draw("same");

    leg->Draw();

    
}