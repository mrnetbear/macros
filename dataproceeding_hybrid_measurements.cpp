#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
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
#include "TMultiGraph.h"

void timeFileRead(std::string filename, std::vector<double>& data, double& mean) {
    std::cout << "Reading file: " << filename << std::endl;
    std::fstream readfile;
    readfile.open(filename, std::ios::in);
    std::string clipboard;
    while (std::getline(readfile, clipboard)) {
        data.push_back(stod(clipboard));
        mean += stod(clipboard);
        if (std::isnan(mean)) {
            throw std::invalid_argument("Error 1: Invalid data in file: " + filename);
        }
    }
    readfile.close();

    sort(data.begin(), data.end());
    short n = data.size();
    mean = mean / n;
}

std::string nameParse(std::string root, size_t chNum, size_t c, size_t r, size_t i){
    std::string name = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/Data/" + root + "/c" + std::to_string(c) 
            + "/r" + std::to_string(r) + "/C"+ std::to_string(chNum) + "--c"+ std::to_string(c) + "r" 
            + std::to_string(r) + "--0";
    if (i < 10)
        name = name + "000" + std::to_string(i);
    else if (i < 100)
        name = name + "00" + std::to_string(i);
    else if (i < 1000)
        name = name + "0" + std::to_string(i);
    else
        name += std::to_string(i);
    name += ".trc";
    return name;
}

//Data pocessing function
int dataProcessing(){
    //Data reading
    std::string filename;

    std::vector <std::vector<TH1D*> > timeHist;
    size_t nBins = 100;

    std::string base_root;

    for (size_t amm = 0; amm < 4; ++amm) {
        switch (amm){
            case 0: base_root = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybrid_measurements/115_460/c"; break;
            case 1: base_root = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybrid_measurements/120_460/c"; break;
            case 2: base_root = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybrid_measurements/123_450/c"; break;
            case 3: base_root = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybrid_measurements/123_460/c"; break;
        }
        std::vector<TH1D*> timeHistTemp;
        for (size_t c = 5; c < 8; ++c) {
            for (size_t r = 6; r < 9; ++r){
                for (size_t pair = 0; pair < 1; ++pair){
                    filename = base_root + std::to_string(c) + "/r" + std::to_string(r);
                    switch (pair){
                        case 0: filename = filename + "/timeHistDeltaTOA23.txt"; break;
                        case 1: filename = filename + "/timeHistDeltaTOA24.txt"; break;
                        case 2: filename = filename + "/timeHistDeltaTOA34.txt"; break;
                    }
                    double  mean = 0.0;
                    std::vector<double> data;
                    try{timeFileRead(filename, data, mean);}
                    catch (const std::invalid_argument& e) {
                        std::cerr << "Error reading file: " << e.what() << std::endl;
                        return 1;
                    }
                    if (data.empty()) {
                        std::cerr << "Error: No data read from file." << std::endl;
                        return 2;
                    }
                    TH1D* hist = new TH1D();
                    std::string histName = std::to_string(amm) + "c" + std::to_string(c) + "r" + std::to_string(r) + "timeHistDeltaTOA";
                    switch (pair){
                        case 0:histName+="23"; break;
                        case 1:histName+="24"; break;
                        case 2:histName+="34"; break;
                    }
                    hist->SetNameTitle(histName.c_str(), histName.c_str());
                    hist->SetBins(nBins, 8000, 11000);
                    std::string histTitle = "Number" + std::to_string(amm) + "#Delta TOA for c" + std::to_string(c) + "r" + std::to_string(r);
                    hist->SetTitle(histTitle.c_str());
                    hist->SetXTitle("Time, ps");
                    hist->SetYTitle("Entries");

                    for (auto a : data){
                        hist->Fill(a);
                    }
                    timeHistTemp.push_back(hist);
                }
            }
        }
        timeHist.push_back(timeHistTemp);
    }

    TFile* file = new TFile("timeHist.root", "RECREATE");
    for (size_t i = 0; i < timeHist.size(); ++i)
        for (size_t j = 0; j < timeHist[i].size(); ++j)
            timeHist[i][j]->Write();
    file->Close();

    TCanvas* c1 = new TCanvas("c1", "timeHist", 800, 600);
    timeHist[0][0]->SetXTitle("Time, ps");
    timeHist[0][0]->SetYTitle("Entries");
    timeHist[0][0]->Draw();

    return 0;
}