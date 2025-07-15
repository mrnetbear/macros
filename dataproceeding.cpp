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

//Fast-Furier transformation
void FFT(const double* input_data, double* output_dataRe, double* output_dataIm, size_t len){
    for (size_t i = 0; i < len; ++i){
        output_dataRe[i] = 0;
        output_dataIm[i] = 0;
        for (size_t j = 0; j < len; ++j){
            output_dataRe[i] += input_data[j] * cos(2 * M_PI * i * j / len);
            output_dataIm[i] += input_data[j] * sin(2 * M_PI * i * j / len);
        }
    }
}

//Data pocessing function
int dataProcessing(){
    //Data reading
    std::fstream readfile;
    std::string filename = "timeCh2Hist.txt";
    readfile.open(filename, std::ios::in);
    std::vector <double> data;
    std::string clipboard;
    while (std::getline(readfile, clipboard)) {
        data.push_back(stod(clipboard));
    }
    readfile.close();

    //Visualizating data
    TH1D* timeHist = new TH1D("timeHist", "#Delta TOA", 50, -2000, -1300);
    for (auto a : data){
        timeHist->Fill(a);
        std::cout << a << std::endl;
    }
    TCanvas* c1 = new TCanvas("c1", "timeHist", 800, 600);
    timeHist->SetXTitle("Time, ps");
    timeHist->SetYTitle("Entries");
    timeHist->Draw();

    //////////Furier processing//////////
    const size_t num_signals = 1087;
    const size_t signal_length = 1002;
    std::vector<std::vector<double>> signals(num_signals, std::vector<double>(signal_length));
    std::vector<std::vector<double>> furier_transformedRe(num_signals, std::vector<double>(signal_length));
    std::vector<std::vector<double>> furier_transformedIm(num_signals, std::vector<double>(signal_length));


    //Reading data
    std::fstream furier_data;
    for (size_t i = 0; i < num_signals; ++i){
        std::string furier_filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed/C2---" + std::to_string(i+1) + ".txt";
        furier_data.open(furier_filename, std::ios::in);
        std::string strBuf;
        size_t j = 0;
        for (size_t j = 0; j < signal_length && std::getline(furier_data, strBuf); ++j) {
            signals[i][j] = (stod(strBuf));

        }
        furier_data.close();
    }

    for (size_t i = 0; i < signal_length; ++i){
        std::cout << signals[0][i] << std::endl;
    }

    std::cout << "Success! Beginning Grid-Stride loop..." << std::endl; 
    //////////grid-stride loop////////////
    
    // Multi-threaded processing
    const size_t num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (size_t i = 0; i < num_signals; i += num_threads) {
        threads.clear();
        
        for (size_t t = 0; t < num_threads && (i + t) < num_signals; ++t) {
            size_t idx = i + t;
            threads.emplace_back([&, idx]() {
                FFT(signals[idx].data(), 
                    furier_transformedRe[idx].data(), 
                    furier_transformedIm[idx].data(), 
                    signal_length);
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
    }


    //Data visualization
    TH1D *furier_spectreRe = new TH1D("furier_spectreRe", "Furier Re Spectre on Ch2", 100, -3500, 4500);
    TH1D *furier_spectreIm = new TH1D("furier_spectreIm", "Furier Im Spectre on Ch2", 100, -3500, 3500);

    for (size_t i = 0; i < signal_length; ++i){
        furier_spectreRe->Fill(furier_transformedRe[0][i]);
        furier_spectreIm->Fill(furier_transformedIm[0][i]);
    }

    TCanvas *c2 = new TCanvas("c2", "furier Spectre", 1024, 768);
    c2->Divide(2, 1);
    c2->cd(1);
    furier_spectreRe->SetXTitle("Frequency, Hz");
    furier_spectreRe->SetYTitle("Entries");
    furier_spectreRe->Draw();
    c2->cd(2);
    furier_spectreIm->SetXTitle("Frequency, Hz");
    furier_spectreIm->SetYTitle("Entries");
    furier_spectreIm->Draw();


    return 0;
}