#include <cmath>

namespace misc
{
    // This function converts a power value to decibels (dB).
    double db(double x)
    {
    return 10.0 * std::log10(x);
    }

    int index_2d_to_1d(int i, int j, int n_cols)
    {
        return i * n_cols + j;
    }
}