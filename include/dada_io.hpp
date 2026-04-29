// This code has not been carefully inspected by the author.

#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

namespace dada {

    namespace fs = std::filesystem;

    // --- 1. DATA STRUCTURE ---
    template <typename T>
    struct PFBInput {
        std::vector<T> data;
        int nchan;
        int npol;
        int n_time;
    };

    // --- 2. PARSER (READ) ---
    template <typename T>
    PFBInput<T> read_dada_for_pfb(const std::string& path) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("Could not open file: " + path);
        
        size_t data_bytes = static_cast<size_t>(f.tellg()) - 4096;
        f.seekg(0, std::ios::beg);

        std::vector<char> h_buf(4096);
        f.read(h_buf.data(), 4096);

        int nchan = 1, npol = 1; 
        std::istringstream iss(std::string(h_buf.data()));
        for (std::string k, v, line; std::getline(iss, line) && line[0] != '\0';) {
            if (std::istringstream(line) >> k >> v) {
                if (k == "NCHAN") nchan = std::stoi(v);
                else if (k == "NPOL") npol = std::stoi(v);
            }
        }

        PFBInput<T> out;
        out.nchan = nchan;
        out.npol = npol;
        out.data.resize(data_bytes / sizeof(T));
        f.read(reinterpret_cast<char*>(out.data.data()), data_bytes);
        out.n_time = out.data.size() / (nchan * npol);

        return out;
    }

    // --- 3. HELPER: BUILD FILEPATHS ---
    inline std::string build_filepath(bool is_input, const std::string& type, int M, int P, int W, bool noise, double freq, int d_per, int d_start) {
        std::string base = is_input ? "/Users/hilays79/Fourier_Space/Data/input_files/" : "/Users/hilays79/Fourier_Space/Data/output_files/c++/";
        std::ostringstream oss;
        std::string noise_str = noise ? "True" : "False";

        if (type == "dirac_deltas") {
            oss << type << "_d" << d_per << "_s" << d_start << "_noise" << noise_str;
        } else {
            oss << type << "_freq" << freq << "_M" << M << "_P" << P << "_W" << W << "_noise" << noise_str;
        }
        return base + type + "/" + oss.str() + ".dada";
    }

    // --- 4. SAVER (WRITE) ---
    template <typename T>
    void save_dada(const std::vector<T>& data, int nchan, int ndim, int nbit, const std::string& path) {
        fs::create_directories(fs::path(path).parent_path()); 
        
        std::string hdr = "HDR_VERSION 1.0\nHDR_SIZE 4096\nNCHAN " + std::to_string(nchan) +
                          "\nNPOL 1\nNDIM " + std::to_string(ndim) + "\nNBIT " + std::to_string(nbit) + "\n";
        hdr.resize(4096, '\0'); 
        
        std::ofstream f(path, std::ios::binary);
        f.write(hdr.data(), 4096);
        f.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(T));
    }

    // --- 5. THE END-TO-END PIPELINE ---
    template <typename InT, typename OutT, typename Func>
    std::string run_pipeline(Func pfb_func, const std::string& type, int M, int P_out, int W, int ndim_out, int nbit, bool noise, double freq = 1.0, int d_per = 0, int d_start = 0) {
        
        std::string in_path = build_filepath(true, type, M, P_out, W, noise, freq, d_per, d_start);
        
        if (!fs::exists(in_path)) {
            std::cerr << "\n[ERROR] GENERATE BINARY USING PYTHON. Missing file:\n" << in_path << "\n\n";
            exit(1); 
        }

        std::cout << "Reading input: " << in_path << "\n";
        auto input = read_dada_for_pfb<InT>(in_path);
        
        std::cout << "Running PFB...\n";
        std::vector<OutT> out_data = pfb_func(input.data, M, P_out, W);

        std::string out_path = build_filepath(false, type, M, P_out, W, noise, freq, d_per, d_start);
        save_dada(out_data, P_out, ndim_out, nbit, out_path);
        
        std::cout << "Saved output: " << out_path << "\nPipeline complete!\n";
        return out_path;
    }

} // namespace dada