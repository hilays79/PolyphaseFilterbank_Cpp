#pragma once
#include <vector>
#include <fftw3.h>

namespace misc
{
    // This function converts a power value to decibels (dB).
    double db(double x);

    inline int index_2d_to_1d(int i, int j, int n_cols) {
        return i * n_cols + j;
    }
}

namespace windowing
{
    std::vector<double> generate_win_coeffs(int n_taps, int n_chan);
}

namespace ts
{
    // Generate a sinusoidal signal with a given n_taps, n_chan, n_windows, frequency, and complex phasor or not.
    std::vector<std::complex<double>> generate_sinusoidal(int n_taps, int n_chan, int n_windows, double omega, bool include_noise, bool complex_phasor, int seed = 42);
}