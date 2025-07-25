#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
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
#include "FFT/FFT.hpp"

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
        std::cout << "Read value: " << stod(clipboard) << std::endl;
        std::cout << "Write value: " << data[data.size()-1] << std::endl;
    }
    readfile.close();

    sort(data.begin(), data.end());
    short n = data.size();
    mean = mean / n;
    std::cout << "Mean value: " << mean << std::endl;
}

//Data pocessing function
int dataProcessing(){

    std::cout << "Choose state: ";
    short state;
    std::cin >> state;
    //Data reading
    std::string filename;

    TH1D* timeHist;
    size_t nBins = 100;

    size_t num_signals;
    size_t signal_length;
    double tstep;

    std::string furierRootFolder;
    
    try{
        switch (state){
            case 0: //20 degrees
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed/timeCh2HistDeltaTOA.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA@20degrees (board 9)", nBins, 1600, 2600);
                num_signals = 1087;
                signal_length = 1002;
                tstep = 50e-12;
                break;
            case 1: //-30 degrees
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/timeCh3Hist@-30.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@-30/C3---";
                timeHist = new TH1D("timeHist", "#Delta TOA@minus30degrees (board 9)", nBins, 200, 1300);
                num_signals = 2006;
                signal_length = 2002;
                tstep = 50e-12;
                break;
            case 2: //no_hybird_HPKHV180
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@no_hybird_HPKHV180/timeCh2HistDeltaTOAno_hybird_HPKHV180.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@no_hybird_HPKHV180/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 0, 2400);
                num_signals = 3335;
                signal_length = 8002;
                tstep = 50e-12;
                break;
            case 3: //no_hybrid_HPKHV190
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@no_hybird_HPKHV190/timeCh2HistDeltaTOAno_hybird_HPKHV190.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@no_hybird_HPKHV190/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 1200, 2400);
                num_signals = 3335;
                signal_length = 8002;
                tstep = 50e-12;
                break;
            case 45:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV180_HPKHV180/r5/timeHistDeltaTOA23hybridHV180_HPKHV180.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV180_HPKHV180/r5/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 4500, 7000);
                num_signals = 146;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 55:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV175_HPKHV180/r5/timeHistDeltaTOA23hybridHV175_HPKHV180.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV175_HPKHV180/r5/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 4500, 7000);
                num_signals = 150;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 65:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV170_HPKHV180/r5/timeHistDeltaTOA23hybridHV170_HPKHV180.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV170_HPKHV180/r7/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 4500, 7000);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 7:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/20degree_HV190/timeHistDeltaTOA2320degree_HV190.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/20degree_HV190/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 700, 1700);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 8:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/20degree_HV180/timeHistDeltaTOA2320degree_HV180.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/20degree_HV180/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 700, 1700);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 9:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/20degree_HV200/timeHistDeltaTOA2320degree_HV200.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/20degree_HV200/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 700, 1700);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 10:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/minus30degree_HV160/timeHistDeltaTOA23minus30degree_HV160.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybridHV170_HPKHV180/r7/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 700, 1700);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 11:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/minus30degree_HV170/timeHistDeltaTOA23minus30degree_HV170.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/minus30degree_HV170/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", 500, 700, 1700);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            case 12:
                filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@hybrid_measurements/115_460/c5/r6/timeHistDeltaTOA23.txt";
                furierRootFolder = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@calibration_7_18/minus30degree_HV170/C2---";
                timeHist = new TH1D("timeHist", "#Delta TOA", nBins, 8000, 10000);
                num_signals = 502;
                signal_length = 8002;
                tstep = 12e-12;
                break;
            default: 
                throw std::invalid_argument("Error 9: Invalid binary digit: " + std::to_string(state));
        }
    }catch(const std::invalid_argument& e) { //catching incorrect symbol exception 
                std::cerr << e.what() << '\n';
                return 9;
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
    double var = 0.0;

    //Visualizating data
    for (auto a : data){
        timeHist->Fill(a);
    }

    size_t n = data.size();

    for (size_t i = 0; i < 100; ++i){
        double value = timeHist->GetBinCenter(i);
        double weight  = timeHist->GetBinContent(i);
        var += (value - mean) * (value - mean) * weight / n;
    }
    
    std::cout << "Mean = " << mean << "; Var = " << sqrt(var) << std::endl;

    TCanvas* c1 = new TCanvas("c1", "timeHist", 800, 600);
    timeHist->SetXTitle("Time, ps");
    timeHist->SetYTitle("Entries");
    timeHist->Draw();

    //////////Furier processing//////////
    std::vector<std::vector<double>> signals(num_signals, std::vector<double>(signal_length));


    //Reading data
    std::fstream furier_data;
    for (size_t i = 0; i < num_signals; ++i){
        std::string furier_filename = furierRootFolder + std::to_string(i+1) + ".txt";
        furier_data.open(furier_filename, std::ios::in);
        std::string strBuf;
        size_t j = 0;
        for (size_t j = 0; j < signal_length && std::getline(furier_data, strBuf); ++j) {
            signals[i][j] = (stod(strBuf));
        }
        furier_data.close();
    }

    std::cout << "Success! Beginning Grid-Stride loop..." << std::endl; 

    //////////grid-stride loop////////////
    
    // Multi-threaded processing
    //std::vector<std::vector<Complex>> fft_results(num_signals);
    std::vector<std::vector<Complex>> fft_results(1);
    std::vector<std::vector<double>> fft_signal_test;
    fft_signal_test.push_back(signals[0]);
    const size_t num_threads = std::thread::hardware_concurrency();
    //ParallelFFT(signals, fft_results, num_threads);
    ParallelFFT(fft_signal_test, fft_results, num_threads);

    //Data visualization
    double sampling_rate = 1 / tstep;

    TH1D *amplitude_spectrum = new TH1D(
        "amplitude", 
        "Amplitude Spectrum; Frequency (Hz); Amplitude", 
        100, 0, sampling_rate/64
    );

    for (size_t i = 0; i < signal_length/2; ++i) {
        double freq = i * sampling_rate / signal_length;
        double amp = std::abs(fft_results[0][i]) / signal_length;
        (amp < 1) ? amp = 0 : amp;
        amplitude_spectrum->Fill(freq, amp);
    }
    TCanvas *c2 = new TCanvas("c2", "Amplitude Spectrum", 800, 600);
    amplitude_spectrum->Draw("hist");


    return 0;
}