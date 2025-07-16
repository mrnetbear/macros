#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>

#define t095 1.9602


//t-criteria calculation
double calculateTTest(double correlation, int nPoints) {
    if (correlation == 1.0)
        return 0; 
    double t = correlation * TMath::Sqrt(nPoints - 2) / TMath::Sqrt(1 - correlation * correlation);
    return t;
}

void manualLinearFit(double *x, double *y, int nPoints, double &a, double &b) {
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    for (int i = 0; i < nPoints; ++i) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }
    a = (nPoints * sumXY - sumX * sumY) / (nPoints * sumX2 - sumX * sumX);
    b = (sumY - a * sumX) / nPoints;
}

void randomDistribution() {
    TRandom3 rand;
    const int nPoints = 100000;

    //arrays of points
    double x[nPoints], y[nPoints];

    //generating numbers
    for (int i = 0; i < nPoints; ++i) {
        x[i] = rand.Uniform(-5, 5);  
        y[i] = rand.Uniform(-5, 5);  
    }

    //hist construction
    TH2D *h1 = new TH2D("h1", "Random generation", 50, -5, 5, 50, -5, 5);
    for (int i = 0; i < nPoints; ++i) {
        h1->Fill(x[i], y[i]);
    }

    //linear fit calculation
    double a1, b1;
    manualLinearFit(x, y, nPoints, a1, b1);

    //approximmation construction
    TF1 *f1 = new TF1("f1", "[0]+[1]*x", -5, 5);
    f1->SetParameters(a1, b1);

    //covariation and correlation calculation
    double cov1 = 0, corr1 = 0;
    double meanX1 = 0, meanY1 = 0;

    for (int i = 0; i < nPoints; ++i) {
        meanX1 += x[i];
        meanY1 += y[i];
    }
    meanX1 /= nPoints;
    meanY1 /= nPoints;

    for (int i = 0; i < nPoints; ++i) {
        cov1 += (x[i] - meanX1) * (y[i] - meanY1);
    }
    cov1 /= nPoints;

    double stdX1 = 0, stdY1 = 0;

    for (int i = 0; i < nPoints; ++i) {
        stdX1 += (x[i] - meanX1) * (x[i] - meanX1);
        stdY1 += (y[i] - meanY1) * (y[i] - meanY1);
    }

    stdX1 = TMath::Sqrt(stdX1 / nPoints);
    stdY1 = TMath::Sqrt(stdY1 / nPoints);
    corr1 = cov1 / (stdX1 * stdY1);

    double t1 = calculateTTest(corr1, nPoints);

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Case 1: Covariance = " << cov1 << ", Correlation = " << corr1 << ", t-statistic = " << t1 << std::endl;
    std::cout << "t095 = " << (t095 < abs(t1) ? "OK!" : "NOPE!") << std::endl;
    std::cout << "Ap function: y = " << a1 << "x + " << b1 << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    TCanvas *c1 = new TCanvas("c1", "Case 1", 800, 600);
    h1->Draw("colz");
    f1->Draw("same L");
}

void polynomialDistribution(){
    
    TRandom3 rand;
    const int nPoints = 100000;

    //arrays of points
    double x[nPoints], y[nPoints];

    //generating numbers 
    for (int i = 0; i < nPoints; ++i) {
        x[i] = rand.Uniform(-5, 5); 
        y[i] = -5 * x[i] * x[i] - 3 * x[i] + 0.6;
    }

    //hist construction
    TH2D *h2 = new TH2D("h2", "Polynomial generation", 50, -5, 5, 50, -5, 5);
    for (int i = 0; i < nPoints; ++i) {
        h2->Fill(x[i], y[i]);
    }

    //linear fit calculation
    double a2, b2;
    manualLinearFit(x, y, nPoints, a2, b2);

    //approximmation construction
    TF1 *f2 = new TF1("f2", "[0]+[1]*x", -5, 5);
    f2->SetParameters(a2, b2);

    //covariation and correlation calculation
    double cov1 = 0, corr1 = 0;
    double meanX1 = 0, meanY1 = 0;

    for (int i = 0; i < nPoints; ++i) {
        meanX1 += x[i];
        meanY1 += y[i];
    }
    meanX1 /= nPoints;
    meanY1 /= nPoints;

    for (int i = 0; i < nPoints; ++i) {
        cov1 += (x[i] - meanX1) * (y[i] - meanY1);
    }
    cov1 /= nPoints;

    double stdX1 = 0, stdY1 = 0;
    for (int i = 0; i < nPoints; ++i) {
        stdX1 += (x[i] - meanX1) * (x[i] - meanX1);
        stdY1 += (y[i] - meanY1) * (y[i] - meanY1);
    }
    stdX1 = TMath::Sqrt(stdX1 / nPoints);
    stdY1 = TMath::Sqrt(stdY1 / nPoints);
    corr1 = cov1 / (stdX1 * stdY1);

    double t1 = calculateTTest(corr1, nPoints);

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Case 2: Covariance = " << cov1 << ", Correlation = " << corr1 << ", t-statistic = " << t1 << std::endl;
    std::cout << "t095 = " << (t095 < abs(t1) ? "OK!" : "NOPE!") << std::endl;
    std::cout << "Ap function: y = " << a2 << "x + " << b2 << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    TCanvas *c2 = new TCanvas("c2", "Case 2", 800, 600);
    h2->Draw("colz");
    f2->Draw("same L");
}

void triangleDistribution(){
    
    TRandom3 rand;
    const int nPoints = 100000;

    
    double x[nPoints], y[nPoints];

    
    for (int i = 0; i < nPoints; ++i) {
        x[i] = rand.Uniform(0, 10);  
        double z = rand.Uniform(0, 10);
        y[i] = z * x[i];
    }


    TH2D *h3 = new TH2D("h3", "Triangle generation", 50, 0, 10, 50, 0, 10);
    for (int i = 0; i < nPoints; ++i) {
        h3->Fill(x[i], y[i]);
    }

    
    //linear fit calculation
    double a3, b3;
    manualLinearFit(x, y, nPoints, a3, b3);

    //approximmation construction
    TF1 *f3 = new TF1("f3", "[0]+[1]*x", 0, 10);
    f3->SetParameters(a3, b3);

    
    double cov1 = 0, corr1 = 0;
    double meanX1 = 0, meanY1 = 0;

    for (int i = 0; i < nPoints; ++i) {
        meanX1 += x[i];
        meanY1 += y[i];
    }
    meanX1 /= nPoints;
    meanY1 /= nPoints;

    for (int i = 0; i < nPoints; ++i) {
        cov1 += (x[i] - meanX1) * (y[i] - meanY1);
    }
    cov1 /= nPoints;

    double stdX1 = 0, stdY1 = 0;
    for (int i = 0; i < nPoints; ++i) {
        stdX1 += (x[i] - meanX1) * (x[i] - meanX1);
        stdY1 += (y[i] - meanY1) * (y[i] - meanY1);
    }
    stdX1 = TMath::Sqrt(stdX1 / nPoints);
    stdY1 = TMath::Sqrt(stdY1 / nPoints);
    corr1 = cov1 / (stdX1 * stdY1);

    double t1 = calculateTTest(corr1, nPoints);

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Case 3: Covariance = " << cov1 << ", Correlation = " << corr1 << ", t-statistic = " << t1 << std::endl;
    std::cout << "t095 = " << (t095 < abs(t1) ? "OK!" : "NOPE!") << std::endl;
    std::cout << "Ap function: y = " << a3 << "x + " << b3 << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    
    TCanvas *c3 = new TCanvas("c3", "Case 3", 800, 600);
    h3->Draw("colz");
    f3->Draw("same L");
}

void circleDistribution(){
  
    TRandom3 rand;
    const int nPoints = 100000;

   
    double x[nPoints], y[nPoints];

  
    for (int i = 0; i < nPoints; ++i) {
        double theta = rand.Uniform(0, 2 * TMath::Pi());
        x[i] = 2 * TMath::Cos(theta);
        y[i] = 2 * TMath::Sin(theta);
    }

   
    TH2D *h4 = new TH2D("h4", "Circle generation", 50, -5, 5, 50, -5, 5);
    for (int i = 0; i < nPoints; ++i) {
        h4->Fill(x[i], y[i]);
    }

  
    //linear fit calculation
    double a4, b4;
    manualLinearFit(x, y, nPoints, a4, b4);

    //approximmation construction
    TF1 *f4 = new TF1("f3", "[0]+[1]*x", -5, 5);
    f4->SetParameters(a4, b4);

   
    double cov1 = 0, corr1 = 0;
    double meanX1 = 0, meanY1 = 0;
    for (int i = 0; i < nPoints; ++i) {
        meanX1 += x[i];
        meanY1 += y[i];
    }
    meanX1 /= nPoints;
    meanY1 /= nPoints;
    for (int i = 0; i < nPoints; ++i) {
        cov1 += (x[i] - meanX1) * (y[i] - meanY1);
    }
    cov1 /= nPoints;
    double stdX1 = 0, stdY1 = 0;
    for (int i = 0; i < nPoints; ++i) {
        stdX1 += (x[i] - meanX1) * (x[i] - meanX1);
        stdY1 += (y[i] - meanY1) * (y[i] - meanY1);
    }
    stdX1 = TMath::Sqrt(stdX1 / nPoints);
    stdY1 = TMath::Sqrt(stdY1 / nPoints);
    corr1 = cov1 / (stdX1 * stdY1);

    double t1 = calculateTTest(corr1, nPoints);

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Case 4: Covariance = " << cov1 << ", Correlation = " << corr1 << ", t-statistic = " << t1 << std::endl;
    std::cout << "t095 = " << (t095 < abs(t1) ? "OK!" : "NOPE!") << std::endl;
    std::cout << "Ap function: y = " << a4 << "x + " << b4 << std::endl;
    std::cout << "---------------------------------------------" << std::endl;


   
    TCanvas *c4 = new TCanvas("c4", "Case 4", 800, 600);
    h4->Draw("colz");
    f4->Draw("same L");
}

int main() {
    randomDistribution();
    polynomialDistribution(); 
    triangleDistribution();
    circleDistribution();
    return 0;
}