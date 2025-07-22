#include "FFT.hpp"

//Fast-Furier transformation
void FFT(std::vector<Complex>& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    //Separating data
    std::vector<Complex> even(N/2), odd(N/2);
    for (size_t i = 0; i < N/2; ++i) {
        even[i] = x[2*i];
        odd[i] = x[2*i + 1];
    }

    //Recource FFT call
    FFT(even);
    FFT(odd);

    //FFT procedure
    for (size_t k = 0; k < N/2; ++k) {
        Complex t = std::exp(Complex(0, -2 * M_PI * k / N)) * odd[k];
        x[k] = even[k] + t;
        if (abs(x[k].real()) < 1e-10) x[k].real(0);
        if (abs(x[k].imag()) < 1e-10) x[k].imag(0);
        x[k + N/2] = even[k] - t;
        if (abs(x[k + N/2].real()) < 1e-10) x[k + N/2].real(0);
        if (abs(x[k + N/2].imag()) < 1e-10) x[k + N/2].imag(0);
    }
}

//one-thread FFT call
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
) {
    const size_t num_signals = signals.size();
    fft_results.resize(num_signals);

    size_t N = signals[0].size();
    size_t N_padded = 1;
    while (N_padded < N) N_padded <<= 1;

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (size_t t = 0; t < num_threads; ++t) {
        threads.emplace_back([&, t, N_padded]() {
            for (size_t i = t; i < num_signals; i += num_threads) {
                ProcessSignal(signals[i], fft_results[i], N_padded);
            }
        });
    }

    for (auto& thread : threads) {
        if (thread.joinable()) thread.join();
    }
}

