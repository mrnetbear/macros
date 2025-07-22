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
#include "FFT/Fft.h"

void testFunction(){
    std::vector<double> signal;
    double sampling_frequency = 256*128;
    double signal_frequency = 256;
    size_t N = sampling_frequency / signal_frequency;
    double t_step = 1 / sampling_frequency;
    double f_step = sampling_frequency / N;
    double time[N], sinus[N], sinus_dad[N], freq[N];
    for (size_t i = 0; i < N; ++i){
        time[i] = i * t_step;
        freq[i] = i * f_step;
        signal.push_back(sin(2 * M_PI * signal_frequency * time[i]) + sin(M_PI * signal_frequency * time[i]));
        std::cout << "y[" << i << "] = " << signal[i] << std::endl;
        sinus[i] = signal[i];
    }   

    std::cout << "=========================" << std::endl;
    std::vector<Complex> fft_result;
    ProcessSignal(signal, fft_result, signal.size());
    
    for(int i = 0; i < 2*N; i+=2){
        sinus_dad[i] = signal[i];
        sinus_dad[i+1] = 0;
    }

    fft(sinus_dad, N);

    TCanvas *c1 = new TCanvas("c1", "Spectres",600 , 1024);
    c1->Divide(2, 2);
    TGraph *time_spectrum = new TGraph(N, time, sinus);
    c1->cd(1);
    time_spectrum->SetTitle("Amplitude vs time");
    time_spectrum->Draw("ALP");

    c1->cd(2);

    double fft_my[N], fft_dad[N];
    for(int i = 0; i < N; ++i){
        fft_my[i] = (i) ? 2 * abs(fft_result[i]) / N : abs(fft_result[i])/N;
    }

    TGraph *freq_spectrum_my = new TGraph(N/2, freq, fft_my);
    freq_spectrum_my->SetTitle("Amplitude vs frequency (mine)");
    freq_spectrum_my->Draw("ALP");

    c1->cd(3);
    //TGraph *freq_spectrum = new TGraph(100, freq, );

    TGraph *freq_spectrum_dad = new TGraph(N/2, freq, sinus_dad);
    freq_spectrum_my->SetTitle("Amplitude vs frequency (dad's)");
    freq_spectrum_my->Draw("ALP");
    TH1D *amplitude_spectrum = new TH1D(
        "amplitude", 
        "Amplitude Spectrum; Frequency (Hz); Amplitude", 
        100, 0, 1000
    );

    for (size_t i = 0; i < N/2; ++i) {
        //std::cout << fft_result[i] << std::endl;
        amplitude_spectrum->Fill(freq[i],fft_my[i]);
    }

    c1->cd(4);
    amplitude_spectrum->Draw("hist");

    std::cout << "========================" << std::endl;

    for (size_t i = 0; i < N; i+=2){
        //std::cout << sinus_dad[i] << "; " << sinus_dad[i+1] << std::endl;
    }

    //amplitude_spectrum -> Draw();
}