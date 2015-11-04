#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

#include "common.h"
#include "Adaptive2DClenshaw.h"

using namespace bemquad;

uint GLOBAL_FUNC_CALLS = 0;

struct Functor
{
    std::complex<double> operator()(const double x, const double y)
    {
        GLOBAL_FUNC_CALLS += 1;
        
        const double cx = -1.0;
        const double cy = -1.0;
        const double k = 101.3;
        const double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
        const auto i = std::complex<double>(0.0, 1.0);
        return std::exp(i * k * r) / (r * r);
    }
};

int main()
{
    double tol = 1.0e-2;
    while(tol > 1.0e-7) {
        GLOBAL_FUNC_CALLS = 0;
        std::cout << "Result = " << adaptive2DClenshaw(0.0, PI, 0.0, PI, Functor(), tol) << "\n";
        std::cout << "Converged to result with " << GLOBAL_FUNC_CALLS << " function calls\n";
        tol *= 0.1;
    }
    std::cout << "done!\n";
    return 0;
}
