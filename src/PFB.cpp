#include <iostream>
# include <vector>
#include <complex>
#include "dsp.hpp"

std::complex<double> filtering(std::complex<double> signal, int n_taps, int n_chan, int n_windows)
{
    
    return signal; // Return the input signal unchanged for now.
}
int main()
{
    int n_taps = 4; // Number of taps in the filter
    int n_chan = 16; // Number of channels
    int n_windows = 1; // Number of windows
    double omega = M_PI / 2; // Frequency of the sinusoidal signal
    bool include_noise = false; // Whether to include noise in the generated signal
    bool complex_phasor = true; // Whether to generate a complex phasor or a real sinusoidal signal
    std::vector<std::complex<double>> signal = ts::generate_sinusoidal(n_taps, n_chan, n_windows, omega, include_noise, complex_phasor);

    return 0;
}