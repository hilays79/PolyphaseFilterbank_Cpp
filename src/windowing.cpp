#include <cmath>

namespace windowing
{
    void FIR_sinc(double* data, int n_taps, int n_chan)
    {
        int n_samples = n_taps*n_chan;
        double fc = 1/n_chan; // Cutoff frequency for bandpass filter
        double center = (n_samples-1) / 2; // Center of the filter
        for (int i = 0; i < n_samples; ++i)
        {
            data[i] = std::sin(2 * M_PI * fc * (i - center)) / (M_PI * (i - center)); // Sinc function for FIR filter coefficients
        }
    }

    // This function applies a Hamming window to the input data.
    void hamming_window(double* data, size_t size) // double* data is a pointer to the input array, size is the number of elements in the array
    {
        for (size_t i = 0; i < size; ++i)
        {
            data[i] *= 0.54 - 0.46 * std::cos(2 * M_PI * i / (size - 1)); // M_PI is the value of pi.
        }
    }

    void hann_window(double* data, size_t size) // double* data is a pointer to the input array, size is the number of elements in the array
    {
        for (size_t i = 0; i < size; ++i)
        {
            data[i] *= 0.5 * (1 - std::cos(2 * M_PI * i / (size - 1))); // M_PI is the value of pi.
        }
    }
}