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

//define step size
const double stepSize = 1e-5;


double f(double x){
        return TMath::Exp(-x);
}

double qGrid(double h){
        return h/pow(1-h,3);
}

void centralIntegrator(double *result, double stepSize){
    double integral = 0;
    for(double xi = stepSize; xi < 1; xi += stepSize){
        double h = qGrid(xi) - qGrid(xi-stepSize);
        integral += f(qGrid(xi)-h/2);
    }
    *result = integral;
        
}

void numIntMore(){
    double ctrInt;
    std::thread th1(centralIntegrator, &ctrInt, stepSize);
    th1.join();

    std::cout << "Numerical Integration Result: " << ctrInt << std::endl;
}