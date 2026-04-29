#include <iostream>
# include <vector>
#include <complex>
#include <fftw3.h>
#include "dsp.hpp"
#include "dada_io.hpp"
#include <chrono>

std::vector<std::complex<double>> filtering(std::vector<std::complex<double>>& signal, int n_taps, int n_chan, int n_windows)
{
    // setup start
    std::cout << "Starting filtering..." << std::endl;
    std::vector<double> win_coeffs = windowing::generate_win_coeffs(n_taps, n_chan);
    int n_time_blocks = n_taps*n_windows - n_taps + 1; // Number of time blocks we can convolve with the window coefficients

    // Each time block will have n_chan samples, and we will have n_time_blocks of them. So the output will be n_chan x n_time_blocks.
    std::vector<std::complex<double>> filtered_signal(n_time_blocks * n_chan); // Initialize the filtered signal with zeros

    std::cout << "Filtered signal size: " << filtered_signal.size() << std::endl;

    //setup end
    //execution start
    // Now perform the convolution of the input signal with the window coefficients.
    for (int n_t = 0; n_t < n_time_blocks; ++n_t) { // Loop over time blocks
        for (int n_c = 0; n_c < n_chan; ++n_c) { // Loop over channels
            std::complex<double> tap_sum = (0.0, 0.0); // Initialize the sum for the tap this time block and channel
            for (int m = 0; m < n_taps; ++m) { // Loop over taps

                int s_index = misc::index_2d_to_1d(n_t + m, n_c, n_chan); // Calculate the index for the input signal
                int w_index = misc::index_2d_to_1d(m, n_c, n_chan); // Calculate the index for the window coefficients
                tap_sum += signal[s_index] * win_coeffs[w_index]; // Accumulate the weighted sum for this tap
                
            }
            filtered_signal[misc::index_2d_to_1d(n_t, n_c, n_chan)] = tap_sum; // Store the result in the filtered signal
        }
    }
    // execution end
    std::cout << "Filtering completed." << std::endl;
    return filtered_signal;
}

std::vector<std::complex<double>> FFT(std::vector<std::complex<double>>& filtered_signal, int n_taps, int n_chan, int n_windows)
{
    // setup start
    std::cout << "Starting FFT..." << std::endl;
    
    // Calculate the number of time blocks to know how many FFTs we need to do
    int n_time_blocks = n_taps * n_windows - n_taps + 1; 

    // Copy the filtered signal. We will perform the FFT on this new vector.
    std::vector<std::complex<double>> x_pfb = filtered_signal;

    // Cast the std::complex pointer to the fftw_complex pointer required by FFTW (which is a C library)
    auto* data_ptr = reinterpret_cast<fftw_complex*>(x_pfb.data());

    // Create a plan for the FFT. We fft of size n_chan.
    fftw_plan plan = fftw_plan_dft_1d(n_chan, data_ptr, data_ptr, FFTW_FORWARD, FFTW_ESTIMATE);
    // POSSIBLE OPTIMISATION: A single plan can be created here using fftw_plan_many_dft.
    // POSSIBLE OPTIMISATION: FFTW_ESTIMATE can be replaced with FFTW_MEASURE for better performance at the cost of longer planning time.

    // setup end
    // Now perform the FFT for each time block. Each time block has n_chan samples, and we will have n_time_blocks of them.
    for (int n_t = 0; n_t < n_time_blocks; ++n_t) {
        
        // Find the 1D starting index for the current time block (equivalent to row start)
        int offset = misc::index_2d_to_1d(n_t, 0, n_chan);
        
        // Get the pointer for this specific row
        fftw_complex* row_ptr = &data_ptr[offset];
        
        // Execute the FFT in-place for this time block
        fftw_execute_dft(plan, row_ptr, row_ptr);
    }
    // Clean up the FFTW plan
    // execution end
    fftw_destroy_plan(plan);
    std::cout << "FFT completed." << std::endl;
    return x_pfb;
}

std::vector<double> PSD(std::vector<std::complex<double>>& x_pfb, int n_taps, int n_chan, int n_windows, int n_integrations=1)
{
    std::cout << "Starting PSD calculation..." << std::endl;
    int n_time_blocks = n_taps * n_windows - n_taps + 1;
    // Trim for integration (matches valid_length calculation in Python)
    int valid_time_blocks = (n_time_blocks / n_integrations) * n_integrations; // This ensures we only consider complete integration blocks
    int n_integrated_blocks = valid_time_blocks / n_integrations; // Number of blocks after integration
    std::vector<double> psd(n_integrated_blocks * n_chan);
    // execution start
    for (int i = 0; i < valid_time_blocks; ++i) {
        int ind_integration_block = i / n_integrations; // Determine which integration block this time block belongs to
        for (int j = 0; j < n_chan; ++j) {
            int index_integration = misc::index_2d_to_1d(ind_integration_block, j, n_chan);
            int index_time_block = misc::index_2d_to_1d(i, j, n_chan);
            psd[index_integration] += std::norm(x_pfb[index_time_block])/n_integrations; // Power is the squared magnitude of the complex number
            // POSSIBLE OPTIMISATION: number of divisions could be reduced by another loop.
        }
    }
    // execution end
    std::cout << "PSD calculation completed." << std::endl;
    return psd;
}

std::vector<double> PFB_filterbank(std::vector<std::complex<double>>& signal, int n_taps, int n_chan, int n_windows, int n_integrations=1)
{
    std::vector<std::complex<double>> filtered_signal = filtering(signal, n_taps, n_chan, n_windows);
    std::vector<std::complex<double>> x_pfb = FFT(filtered_signal, n_taps, n_chan, n_windows);
    std::vector<double> psd = PSD(x_pfb, n_taps, n_chan, n_windows, n_integrations);
    return psd;
}

int main(int argc, char* argv[]) {
    int M = 4, P = 256;
    
    // 1. Allow Python to pass the Window size via command line
    int W = 100;
    if (argc > 1) {
        W = std::stoi(argv[1]);
    }

    double freq = 1.0;
    int nbit = 64, ndim_out = 1; 
    bool include_noise = false;
    std::string signal_type = "complex_phasors"; // or "complex_phasors" or "dirac_deltas" or "sinusoidals"
    int delta_period = 257, delta_start = 0;

    // 2. Add strict timing around the math function ONLY
    auto my_pfb = [](std::vector<std::complex<double>>& d, int m, int p, int w) {
        auto start = std::chrono::high_resolution_clock::now();
        
        auto result = PFB_filterbank(d, m, p, w); // The actual math
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        
        // Print this exact flag so Python can scrape the time
        std::cout << "CPP_MATH_TIME:" << diff.count() << "\n";
        
        return result;
    };

    try {
        dada::run_pipeline<std::complex<double>, double>(
            my_pfb, signal_type, M, P, W, ndim_out, nbit, include_noise, freq, delta_period, delta_start
        );
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << "\n";
    }

    return 0;
}