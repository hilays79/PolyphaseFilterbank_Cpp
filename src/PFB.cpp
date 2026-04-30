#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include "dsp.hpp"
#include "dada_io.hpp"
#include "FFTW.hpp"
#include <chrono>

template <typename T> // Template function declarations where T is either double or float.
std::vector<std::complex<T>> filtering(std::vector<std::complex<T>>& signal, int n_taps, int n_chan, int n_windows, double& setup_time, double& exec_time)
{
    auto s_start = std::chrono::high_resolution_clock::now();
    
    std::vector<T> win_coeffs = windowing::generate_win_coeffs<T>(n_taps, n_chan); 
    int n_time_blocks = n_taps*n_windows - n_taps + 1; 
    std::vector<std::complex<T>> filtered_signal(n_time_blocks * n_chan); 

    auto s_end = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration<double>(s_end - s_start).count();

    auto e_start = std::chrono::high_resolution_clock::now();
    
    for (int n_t = 0; n_t < n_time_blocks; ++n_t) {
        int out_offset = misc::index_2d_to_1d(n_t, 0, n_chan);
        for (int m = 0; m < n_taps; ++m) {
            int w_offset = misc::index_2d_to_1d(m, 0, n_chan);
            int s_offset = misc::index_2d_to_1d(n_t + m, 0, n_chan);
            for (int n_c = 0; n_c < n_chan; ++n_c) {
                filtered_signal[out_offset + n_c] += signal[s_offset + n_c] * win_coeffs[w_offset + n_c];
            }
        }
    }
    
    auto e_end = std::chrono::high_resolution_clock::now();
    exec_time += std::chrono::duration<double>(e_end - e_start).count();

    return filtered_signal;
}

template <typename T>
void FFT(std::vector<std::complex<T>>& filtered_signal, int n_taps, int n_chan, int n_windows, double& setup_time, double& exec_time)
{
    auto s_start = std::chrono::high_resolution_clock::now();
    int n_time_blocks = n_taps * n_windows - n_taps + 1; 

    // Use the alias for the pointer cast (required)
    auto* data_ptr = reinterpret_cast<typename FFTWWrapper<T>::complex_type*>(filtered_signal.data());
    
    int n[] = {n_chan};

    // Use 'auto' for the plan (it deduces fftw_plan or fftwf_plan automatically)
    auto plan = FFTWWrapper<T>::plan_many_dft(1, n, n_time_blocks,
                                              data_ptr, NULL, 1, n_chan,
                                              data_ptr, NULL, 1, n_chan,
                                              FFTW_FORWARD, FFTW_ESTIMATE);
    
    auto s_end = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration<double>(s_end - s_start).count();

    auto e_start = std::chrono::high_resolution_clock::now();
    
    FFTWWrapper<T>::execute(plan);
    
    auto e_end = std::chrono::high_resolution_clock::now();
    exec_time += std::chrono::duration<double>(e_end - e_start).count();

    FFTWWrapper<T>::destroy_plan(plan);
}

template <typename T>
std::vector<T> PSD(std::vector<std::complex<T>>& x_pfb, int n_taps, int n_chan, int n_windows, int n_integrations, double& setup_time, double& exec_time)
{
    auto s_start = std::chrono::high_resolution_clock::now();
    
    int n_time_blocks = n_taps * n_windows - n_taps + 1;
    int valid_time_blocks = (n_time_blocks / n_integrations) * n_integrations; 
    int n_integrated_blocks = valid_time_blocks / n_integrations; 
    std::vector<T> psd(n_integrated_blocks * n_chan);
    
    auto s_end = std::chrono::high_resolution_clock::now();
    setup_time += std::chrono::duration<double>(s_end - s_start).count();

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

    return psd;
}

template <typename T>
std::vector<T> PFB_filterbank(std::vector<std::complex<T>>& signal, int n_taps, int n_chan, int n_windows, int n_integrations=1)
{
    double setup_time = 0.0;
    double exec_time = 0.0;

    std::vector<std::complex<T>> filtered_signal = filtering<T>(signal, n_taps, n_chan, n_windows, setup_time, exec_time);
    FFT<T>(filtered_signal, n_taps, n_chan, n_windows, setup_time, exec_time); // In-place FFT, so we don't need to capture the return value.
    std::vector<T> psd = PSD<T>(filtered_signal, n_taps, n_chan, n_windows, n_integrations, setup_time, exec_time);

    std::cout << "CPP_SETUP_TIME:" << setup_time << "\n";
    std::cout << "CPP_EXEC_TIME:" << exec_time << "\n";

    return psd;
}

int main(int argc, char* argv[]) {
    int M = 4, P = 256;
    int W = 100;
    
    // --- COMMAND LINE ARGUMENTS ---
    int in_NBIT = 64;
    int out_NBIT = 32;

    if (argc > 1) W = std::stoi(argv[1]);
    if (argc > 2) in_NBIT = std::stoi(argv[2]);
    if (argc > 3) out_NBIT = std::stoi(argv[3]);

    double freq = 1.0;
    int ndim_out = 1; 
    bool include_noise = false;
    std::string signal_type = "complex_phasors"; 
    int delta_period = 257, delta_start = 0;

    try {
        // SCENARIO 1: 64-bit Input -> 64-bit Output
        if (in_NBIT == 64 && out_NBIT == 64) {
            std::cout << "Reading 64-bit | Math 64-bit\n";
            auto my_pfb = [](std::vector<std::complex<double>>& d, int m, int p, int w) {
                auto start = std::chrono::high_resolution_clock::now();
                auto result = PFB_filterbank<double>(d, m, p, w); 
                auto end = std::chrono::high_resolution_clock::now();
                std::cout << "CPP_MATH_TIME:" << std::chrono::duration<double>(end - start).count() << "\n";
                return result;
            };
            dada::run_pipeline<std::complex<double>, double>(my_pfb, signal_type, in_NBIT, out_NBIT, M, P, W, ndim_out, include_noise, freq, delta_period, delta_start);
        } 
        
        // SCENARIO 2: 64-bit Input -> 32-bit Output (Downcast)
        else if (in_NBIT == 64 && out_NBIT == 32) {
            std::cout << "Reading 64-bit | Math 32-bit (Downcasting)\n";
            auto my_pfb = [](std::vector<std::complex<double>>& d, int m, int p, int w) {
                std::vector<std::complex<float>> d_float(d.begin(), d.end()); // Safe Downcast
                auto start = std::chrono::high_resolution_clock::now();
                auto result = PFB_filterbank<float>(d_float, m, p, w); 
                auto end = std::chrono::high_resolution_clock::now();
                std::cout << "CPP_MATH_TIME:" << std::chrono::duration<double>(end - start).count() << "\n";
                return result;
            };
            dada::run_pipeline<std::complex<double>, float>(my_pfb, signal_type, in_NBIT, out_NBIT, M, P, W, ndim_out, include_noise, freq, delta_period, delta_start);
        }

        // SCENARIO 3: 32-bit Input -> 32-bit Output (Native 32-bit)
        else if (in_NBIT == 32 && out_NBIT == 32) {
            std::cout << "Reading 32-bit | Math 32-bit\n";
            auto my_pfb = [](std::vector<std::complex<float>>& d, int m, int p, int w) {
                auto start = std::chrono::high_resolution_clock::now();
                auto result = PFB_filterbank<float>(d, m, p, w); 
                auto end = std::chrono::high_resolution_clock::now();
                std::cout << "CPP_MATH_TIME:" << std::chrono::duration<double>(end - start).count() << "\n";
                return result;
            };
            // Notice the template type here is now std::complex<float> for the INPUT!
            dada::run_pipeline<std::complex<float>, float>(my_pfb, signal_type, in_NBIT, out_NBIT, M, P, W, ndim_out, include_noise, freq, delta_period, delta_start);
        }
        
        // Error trap for upcasting (32->64) or invalid NBITs
        else {
            std::cerr << "Fatal Error: Unsupported NBIT combination. Input: " << in_NBIT << ", Output: " << out_NBIT << "\n";
            return 1;
        }

    } catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << "\n";
    }

    return 0;
}