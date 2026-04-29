#include <vector>
#include <cmath>
#include <complex>

namespace windowing {
    void normalise_win(std::vector<double>& win) // The & means vector passed by reference, so the original vector will be modified, and no copy will be made.
    { 
        // Calculate Processing Gain: Sum of (c_i^2)
        double pg = 0.0;
        for (int i = 0; i < win.size(); ++i) {
            pg += win[i] * win[i];
        }

        double norm_factor = std::sqrt(pg);
        // Normalise the window coefficients
        for (int i = 0; i < win.size(); ++i) {
            win[i] /= norm_factor;
        }
    }

    std::vector<double> generate_hamming(int n_taps, int n_chan) {
        int n_samples = n_taps * n_chan;
        std::vector<double> win(n_samples);
        for (int i = 0; i < n_samples; ++i) {
            win[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (n_samples - 1));
        }
        return win;
    }

    std::vector<double> generate_sinc(int n_taps, int n_chan) {
        int n_samples = n_taps * n_chan;
        std::vector<double> sinc(n_samples);
        double fc = 1.0 / (2.0 * n_chan); // Normalized cutoff frequency (1/n_chan)
        double center = (n_samples - 1) / 2.0;

        for (int i = 0; i < n_samples; ++i) {
            double x = i - center;
            if (std::abs(x) < 1e-9) { // This is probably never true as n_chan is even, but included for completeness
                sinc[i] = 2.0 * fc;
            } else {
                sinc[i] = std::sin(2.0 * M_PI * fc * x) / (M_PI * x);
            }
        }
        return sinc;
    }

    std::vector<double> generate_win_coeffs(int n_taps, int n_chan) {
        // Factory pattern: Generate both and combine
        std::vector<double> win = generate_hamming(n_taps, n_chan);
        std::vector<double> sinc = generate_sinc(n_taps, n_chan);

        // Element-wise multiplication (Equivalent to win *= sinc in Python)
        for (size_t i = 0; i < win.size(); ++i) {
            win[i] *= sinc[i];
        }
        normalise_win(win); // Normalise the window coefficients
        return win;
    }
}