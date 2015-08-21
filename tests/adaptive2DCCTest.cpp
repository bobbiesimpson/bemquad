#include <iostream>
#include <iomanip>
#include <fstream>

#include "common.h"
#include "Adaptive2DClenshaw.h"

using namespace bemquad;

uint GLOBAL_FUNC_CALLS = 0;

struct Functor
{
    double operator()(const double x, const double y)
    {
        const double cx = 0.6;
        const double cy = -0.01;
        const double r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
        GLOBAL_FUNC_CALLS += 1;
        return 1.0 / r;
        //return std::exp(std::cos(r * PI));
        return std::exp(std::cos(r * PI)) / r;
    }
};

int main()
{
    std::ofstream ofs("/Users/Robert/Desktop/adaptive2d.dat");
    ofs << std::setprecision(15);
    double tol = 1.0e-2;
    while(tol > 1.0e-7) {
        GLOBAL_FUNC_CALLS = 0;
        const auto r = adaptive2DClenshaw(0.0, 1.0, 0.0, 1.0, Functor(), tol);
        std::cout << std::setprecision(15);
        std::cout << "Converged to result: " << r << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
        ofs << GLOBAL_FUNC_CALLS << "\t" << r << "\n";
        tol *= 0.1;
    }
    std::cout << "done!\n";
    return 0;
}
