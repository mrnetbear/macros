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

using Complex = std::complex<double>;


//Fast-Furier transformation
void FFT(std::vector<Complex>& x){
    const size_t N = x.size();
    if (N <= 1) return;
    
    std::vector<Complex> even(N/2), odd(N/2);
    for (size_t i = 0; i < N / 2; ++i){
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    FFT(even);
    FFT(odd);

    for (size_t k = 0; k < N / 2; ++k) {
        Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

void ProcessSignal( const std::vector<double>& input,
    std::vector<Complex>& output,
    size_t signal_length_padded
){
    //Adding zero values to functuon to achieve 2^n values
    std::vector<Complex> x(signal_length_padded, 0.0);
    for (size_t i = 0; i < input.size(); ++i) {
        x[i] = Complex(input[i], 0.0);
    }

    FFT(x);

    output = x;
}

void ParallelFFT(
    const std::vector<std::vector<double>>& signals,
    std::vector<std::vector<Complex>>& fft_results,
    size_t num_threads
){
    const size_t num_signals = signals.size();
    fft_results.resize(num_signals);

    size_t N = signals[0].size();
    size_t N_padded = 1;
    while (N_padded < N) N_padded <<= 2;

    //////////////grid-stride loop//////////////

    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_signals; i += num_threads) {
        for (size_t t = 0; t < num_threads && (i + t) < num_signals; ++t) {
            size_t idx = i + t;
            threads.emplace_back(
                ProcessSignal,
                std::cref(signals[idx]),
                std::ref(fft_results[idx]),
                N_padded
            );
        }

        // thread join
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
        threads.clear();
    }
}

//Data pocessing function
int dataProcessing(){
    //Data reading
    std::fstream readfile;
    std::string filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/timeCh2HistDeltaTOA.txt"; //20 degrees
    //std::string filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/timeCh3Hist@-30.txt"; //-30 degrees

    double  mean = 0.0,
            meanSquare = 0.0;
    readfile.open(filename, std::ios::in);
    std::vector <double> data;
    std::string clipboard;
    while (std::getline(readfile, clipboard)) {
        data.push_back(stod(clipboard));
        mean += stod(clipboard);
        meanSquare += stod(clipboard) * stod(clipboard);
    }
    readfile.close();

    sort(data.begin(), data.end());
    short n = data.size();
    mean = mean / n;

    double var = 0.0;

    //Visualizating data
    size_t nBins = 100;
    TH1D* timeHist = new TH1D("timeHist", "#Delta TOA", nBins, -2000, -1300); //20 degrees
    //TH1D* timeHist = new TH1D("timeHist", "#Delta TOA", 50, 200, 1300); //-30 degrees
    for (auto a : data){
        timeHist->Fill(a);
    }

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
    const size_t num_signals = 1087; //20 degrees
    const size_t signal_length = 1002; //20 degrees
    //const size_t num_signals = 2007; //-30 degrees
    //const size_t signal_length = 748; // -30 degrees
    std::vector<std::vector<double>> signals(num_signals, std::vector<double>(signal_length));


    //Reading data
    std::fstream furier_data;
    for (size_t i = 0; i < num_signals; ++i){
        std::string furier_filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed/C2---" + std::to_string(i+1) + ".txt"; //20 degrees
        //std::string furier_filename = "/Users/mcsquare/Documents/Работа/2024-2025/FuSEP2025/code/Data_to_Proceed@-30/C3---" + std::to_string(i+1) + ".txt"; //-30 degrees
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
    std::vector<std::vector<Complex>> fft_results(num_signals);
    const size_t num_threads = std::thread::hardware_concurrency();
    ParallelFFT(signals, fft_results, num_threads);

    //Data visualization
    double sampling_rate = 2e10;
    TH1D *amplitude_spectrum = new TH1D(
        "amplitude", 
        "Amplitude Spectrum; Frequency (Hz); Amplitude", 
        signal_length / 2, 0, sampling_rate / 2
    );

    for (size_t i = 0; i < signal_length / 2; ++i) {
        double freq = i * sampling_rate / signal_length;
        double amp = std::abs(fft_results[0][i]);
        amplitude_spectrum->Fill(freq, amp);
    }
    TCanvas *c2 = new TCanvas("c2", "Amplitude Spectrum", 800, 600);
    amplitude_spectrum->Draw("hist");


    return 0;
}