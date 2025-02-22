#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
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

template <class T>
void printVector(std::vector<std::vector<T>> const& vec){
    for (const auto &row : vec) {
        for (const auto &element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

void generateFullRandom(){
    srand(time(NULL));
    std::set<std::pair<int, int>> coordinates;
    TH2I *randHist = new TH2I("Random distribution", "Random distribution", 100, 0, 100, 100, 0, 100);
    TCanvas *c1 = new TCanvas("c1", "Random distribution", 800, 600);
    for (size_t i = 0; i < 1e2; ++i){
        std::pair<int, int> buf (rand() % 100, rand() % 100);
        coordinates.insert(buf);
        //std::cout << buf.first << " " << buf.second << std::endl;
    };

    int iSumX = 0,
        iSumY = 0;
    size_t countier = 0;
    for (auto i : coordinates){
        iSumX += i.first;
        iSumY += i.second;
        ++countier;
        randHist->Fill(i.first, i.second);
    }

    double  duX = iSumX / countier,
            duY = iSumY / countier;
    
    double  dSumXY = 0,
            dSumXX = 0,
            dSumYY = 0;
    for (auto i : coordinates){
        dSumXY += (i.first - duX) * (i.second - duY);
        dSumXX += (i.first - duX) * (i.first - duX);
        dSumYY += (i.second - duY) * (i.second - duY);
    }
    double b1 = dSumXY / dSumXX;
    double b0 = duY - b1 * duX;

    //printVector(point);
    std::cout << "Absolute random generation: y = " << b1 << "x + " << b0 << std::endl;
    randHist->Draw("colz");

}

void generate2ndPol(){
    srand(time(NULL));
    std::set<std::pair<double, double>> coordinates;
    TH2I *pol2nd = new TH2I("2nd polynom", "2nd polynom", 100, 0, 100, 10000, -50000, 1);
    TCanvas *c2 = new TCanvas("c2", "2nd polynom", 800, 600);
    for (size_t i = 0; i < 1e6; ++i){
        double  x = (rand() % 100),
                y = -5*(x*x) - 3*x + 0.6;
        std::pair<double, double> buf (x, y);
        coordinates.insert(buf);
        //std::cout << buf.first << " " << buf.second << std::endl;
    };

    double iSumX = 0,
        iSumY = 0;
    size_t countier = 0;
    for (auto i : coordinates){
        iSumX += i.first;
        iSumY += i.second;
        ++countier;
        pol2nd->Fill(i.first, i.second);
    }

    double  duX = iSumX / countier,
            duY = iSumY / countier;
    
    double  dSumXY = 0,
            dSumXX = 0,
            dSumYY = 0;
    for (auto i : coordinates){
        dSumXY += (i.first - duX) * (i.second - duY);
        dSumXX += (i.first - duX) * (i.first - duX);
        dSumYY += (i.second - duY) * (i.second - duY);
    }
    double b1 = dSumXY / dSumXX;
    double b0 = duY - b1 * duX;

    //printVector(point);
    std::cout << "Polynomial generation: y = " << b1 << "x + " << b0 << std::endl;
    pol2nd->Draw("colz");
}

void generatePseudoLine(){
    srand(time(NULL));
    std::set<std::pair<double, double>> coordinates;
    TH2I *pseudoLine = new TH2I("Triangle", "Triangle", 100, 0, 100, 10000, 0, 10000);
    TCanvas *c3 = new TCanvas("c3", "Triangle", 800, 600);
    for (size_t i = 0; i < 1e6; ++i){
        double  x = (rand() % 100),
                z = (rand() % 100),
                y = z * x;
        std::pair<double, double> buf (x, y);
        coordinates.insert(buf);
        //std::cout << buf.first << " " << buf.second << std::endl;
    };

    double iSumX = 0,
        iSumY = 0;
    size_t countier = 0;
    for (auto i : coordinates){
        iSumX += i.first;
        iSumY += i.second;
        ++countier;
        pseudoLine->Fill(i.first, i.second);
    }

    double  duX = iSumX / countier,
            duY = iSumY / countier;
    
    double  dSumXY = 0,
            dSumXX = 0,
            dSumYY = 0;
    for (auto i : coordinates){
        dSumXY += (i.first - duX) * (i.second - duY);
        dSumXX += (i.first - duX) * (i.first - duX);
        dSumYY += (i.second - duY) * (i.second - duY);
    }
    double b1 = dSumXY / dSumXX;
    double b0 = duY - b1 * duX;

    //printVector(point);
    std::cout << "Triangle generation: y = " << b1 << "x + " << b0 << std::endl;
    pseudoLine->Draw("colz");
} 

void generateCircle(){
    srand(time(NULL));
    std::set<std::pair<double, double>> coordinates;
    TH2I *circle = new TH2I("Circle", "Circle", 200, 0, 4, 200, 0, 4);
    TCanvas *c4 = new TCanvas("c4", "Circle", 800, 800);
    for (size_t i = 0; i < 1e6; ++i){
        double  x = (double)(rand() % 101) / 100.0 * 2.0,
                y = sqrt(4 - x * x);
        std::pair<double, double> buf (x, y);
        coordinates.insert(buf);
        //std::cout << buf.first << " " << buf.second << std::endl;
    };

    double iSumX = 0,
        iSumY = 0;
    size_t countier = 0;
    for (auto i : coordinates){
        iSumX += i.first;
        iSumY += i.second;
        ++countier;
        circle->Fill(i.first, i.second);
    }

    double  duX = iSumX / countier,
            duY = iSumY / countier;
    
    double  dSumXY = 0,
            dSumXX = 0,
            dSumYY = 0;
    for (auto i : coordinates){
        dSumXY += (i.first - duX) * (i.second - duY);
        dSumXX += (i.first - duX) * (i.first - duX);
        dSumYY += (i.second - duY) * (i.second - duY);
    }
    double b1 = dSumXY / dSumXX;
    double b0 = duY - b1 * duX;

    //printVector(point);

    std::cout << "Circle generation: y = " << b1 << "x + " << b0 << std::endl;
    circle->Draw("colz");
}

void generateAll(){

    generateFullRandom();
    generate2ndPol();
    generatePseudoLine();
    generateCircle();

    /*std::thread th1(generateFullRandom);
    std::thread th2(generate2ndPol);
    std::thread th3(generatePseudoLine);
    std::thread th4(generateCircle);

    th1.join();
    th2.join();
    th3.join();
    th4.join();*/
    //system("pause");
}