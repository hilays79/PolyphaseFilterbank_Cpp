#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include "dsp.hpp"
#include "dada_io.hpp"
#include <chrono>

std::vector<std::complex<double>> filtering(std::vector<std::complex<double>>& signal, int n_taps, int n_chan, int n_windows, double& setup_time, double& exec_time)
{
    // --- SETUP START ---
    auto s_start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> win_coeffs = windowing::generate_win_coeffs(n_taps, n_chan);
    int n_time_blocks = n_taps*n_windows - n_taps + 1; 
    std::vector<std::complex<double>> filtered_signal(n_time_blocks * n_chan); 

    auto s_end = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration<double>(s_end - s_start).count();
    // --- SETUP END ---

    // --- EXECUTION START ---
    auto e_start = std::chrono::high_resolution_clock::now();
    
    for (int n_t = 0; n_t < n_time_blocks; ++n_t) {
        // Calculate the base index (column 0) for the output
        int out_offset = misc::index_2d_to_1d(n_t, 0, n_chan);
        
        for (int m = 0; m < n_taps; ++m) {
            // Calculate the base index (column 0) for the window and signal
            int w_offset = misc::index_2d_to_1d(m, 0, n_chan);
            int s_offset = misc::index_2d_to_1d(n_t + m, 0, n_chan);
            
            for (int n_c = 0; n_c < n_chan; ++n_c) {
                // Simply add the column offset (n_c) to the hoisted row offsets
                filtered_signal[out_offset + n_c] += signal[s_offset + n_c] * win_coeffs[w_offset + n_c];
            }
        }
    }
    
    auto e_end = std::chrono::high_resolution_clock::now();
    exec_time += std::chrono::duration<double>(e_end - e_start).count();
    // --- EXECUTION END ---

    return filtered_signal;
}

std::vector<std::complex<double>> FFT(std::vector<std::complex<double>>& filtered_signal, int n_taps, int n_chan, int n_windows, double& setup_time, double& exec_time)
{
    // --- SETUP START ---
    auto s_start = std::chrono::high_resolution_clock::now();
    
    int n_time_blocks = n_taps * n_windows - n_taps + 1; 
    std::vector<std::complex<double>> x_pfb = filtered_signal;
    auto* data_ptr = reinterpret_cast<fftw_complex*>(x_pfb.data());
    // Create a single plan for multiple FFTs
    int n[] = {n_chan};
    fftw_plan plan = fftw_plan_many_dft(1, n, n_time_blocks,
                                        data_ptr, NULL, 1, n_chan,
                                        data_ptr, NULL, 1, n_chan,
                                        FFTW_FORWARD, FFTW_ESTIMATE);
    
    auto s_end = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration<double>(s_end - s_start).count();
    // --- SETUP END ---

    // --- EXECUTION START ---
    auto e_start = std::chrono::high_resolution_clock::now();
    
    fftw_execute(plan);
    
    auto e_end = std::chrono::high_resolution_clock::now();
    exec_time += std::chrono::duration<double>(e_end - e_start).count();
    // --- EXECUTION END ---

    fftw_destroy_plan(plan);
    return x_pfb;
}

std::vector<double> PSD(std::vector<std::complex<double>>& x_pfb, int n_taps, int n_chan, int n_windows, int n_integrations, double& setup_time, double& exec_time)
{
    // --- SETUP START ---
    auto s_start = std::chrono::high_resolution_clock::now();
    
    int n_time_blocks = n_taps * n_windows - n_taps + 1;
    int valid_time_blocks = (n_time_blocks / n_integrations) * n_integrations; 
    int n_integrated_blocks = valid_time_blocks / n_integrations; 
    std::vector<double> psd(n_integrated_blocks * n_chan);
    
    auto s_end = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration<double>(s_end - s_start).count();
    // --- SETUP END ---

    // --- EXECUTION START ---
    auto e_start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < valid_time_blocks; ++i) {
        int ind_integration_block = i / n_integrations; 
        for (int j = 0; j < n_chan; ++j) {
            int index_integration = misc::index_2d_to_1d(ind_integration_block, j, n_chan);
            int index_time_block = misc::index_2d_to_1d(i, j, n_chan);
            psd[index_integration] += std::norm(x_pfb[index_time_block])/n_integrations; 
        }
    }
    
    auto e_end = std::chrono::high_resolution_clock::now();
    exec_time += std::chrono::duration<double>(e_end - e_start).count();
    // --- EXECUTION END ---

    return psd;
}

std::vector<double> PFB_filterbank(std::vector<std::complex<double>>& signal, int n_taps, int n_chan, int n_windows, int n_integrations=1)
{
    double setup_time = 0.0;
    double exec_time = 0.0;

    // Pass the timing variables by reference so they accumulate across all three functions
    std::vector<std::complex<double>> filtered_signal = filtering(signal, n_taps, n_chan, n_windows, setup_time, exec_time);
    std::vector<std::complex<double>> x_pfb = FFT(filtered_signal, n_taps, n_chan, n_windows, setup_time, exec_time);
    std::vector<double> psd = PSD(x_pfb, n_taps, n_chan, n_windows, n_integrations, setup_time, exec_time);

    // Print the sub-times so Python can scrape them
    std::cout << "CPP_SETUP_TIME:" << setup_time << "\n";
    std::cout << "CPP_EXEC_TIME:" << exec_time << "\n";

    return psd;
}

int main(int argc, char* argv[]) {
    int M = 4, P = 256;
    int W = 100;
    if (argc > 1) {
        W = std::stoi(argv[1]);
    }

    double freq = 1.0;
    int nbit = 64, ndim_out = 1; 
    bool include_noise = false;
    std::string signal_type = "complex_phasors"; 
    int delta_period = 257, delta_start = 0;

    auto my_pfb = [](std::vector<std::complex<double>>& d, int m, int p, int w) {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = PFB_filterbank(d, m, p, w); 
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        
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