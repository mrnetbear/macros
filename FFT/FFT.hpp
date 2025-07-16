#ifndef FFT_HPP
#define FFT_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <thread>
#include <chrono>
#include <math.h>

using Complex = std::complex<double>;

void FFT(std::vector<Complex>& x);
void ProcessSignal( const std::vector<double>& input,
    std::vector<Complex>& output,
    size_t signal_length_padded
);
void ParallelFFT(
    const std::vector<std::vector<double>>& signals,
    std::vector<std::vector<Complex>>& fft_results,
    size_t num_threads
);


#endif //FTT_HPP