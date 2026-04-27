#pragma once
#include <vector>

namespace dsp
{
    // This function converts a power value to decibels (dB).
    double db(double x);
}

namespace windowing
{
    // This function applies a Hamming window to the input data.
    void hamming_window(double* data, size_t size); // double* data is a pointer to the input array, size is the number of elements in the array
    void hann_window(double* data, size_t size); // double* data is a pointer to the input array, size is the number of elements in the array
}

namespace ts
{
    // Generate a sinusoidal signal with a given n_taps, n_chan, n_windows, frequency, and complex phasor or not.
    std::vector<std::complex<double>> generate_sinusoidal(int n_taps, int n_chan, int n_windows, double omega, bool include_noise, bool complex_phasor, int seed = 42);
}