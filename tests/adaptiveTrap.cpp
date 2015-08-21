#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "Adaptive1DTrapezoidal.h"
#include "Adaptive1DClenshaw.h"
#include "Adaptive1DGauss.h"

#include "IQuadrature.h"

using namespace bemquad;

uint GLOBAL_FUNC_CALLS = 0;

struct Functor 
{
    double operator()(const double x)
    {
        GLOBAL_FUNC_CALLS += 1;
        //return 1.0;
        return std::exp(std::cos(100.0* PI * x)) * 1.0 / (std::pow((1.004 - x),2));
        //return 1.0 / std::pow((300.0 - x), 2.0);
        //return 1.0 / std::pow((1.004 -x), 2);
    }
};

double error(const double x)
{
  //  const double exact = 249.500998003992;
   // const double exact = 2.5321317555040;
    //const double exact = 2.6020035998653;
//   const double exact = 94.732334765391;
  //const double exact = 611.71659692811;
    const double exact = 416.84975333499;
   //const double exact = 421.67517121718;
        //const double exact = 14.0979937514716;
    return std::abs(exact - x) / std::abs(exact);
}

int main(int argc, char* argv[]) 
{
    std::ofstream trap_ofs("trapezoidal.dat");
    std::ofstream cc_ofs("clenshaw.dat");
    std::ofstream gauss_ofs("gauss.dat");
    if(!trap_ofs || !cc_ofs || !gauss_ofs) {
        throw std::runtime_error("Can't open file for writing");
    }
    

    std::cout << std::setprecision(15);

//    double subresult = 0.0;
//    const double h = 2.0 / 10;
//    double lower = -1.0; double upper = lower + h;
//    while(0.5 * (upper + lower) < 1.0) {
//        subresult += adaptive1DClenshaw(lower, upper, Functor(), 1.0e-2);
//        lower += h;
//        upper += h;
//    }
//    std::cout << GLOBAL_FUNC_CALLS << "\t subresult = " << subresult << "\n";
//    GLOBAL_FUNC_CALLS = 0;
    
//    std::cout << adaptive1DClenshaw(-1.0, 1.0, Functor(), 1.0e-12, 100, 5) << "\n";
//    std::cout << GLOBAL_FUNC_CALLS << "\n"; GLOBAL_FUNC_CALLS = 0;
//    std::cout << adaptive1DGauss(-1.0, 1.0, Functor(), 1.0e-3, 1000) << "\n";
//    std::cout << GLOBAL_FUNC_CALLS << "\n"; GLOBAL_FUNC_CALLS = 0;
    
    //  Simple convergence testing
//    Functor f;
//    const uint nmax = 1000;
//    for(uint n = 2; n < nmax; n += 1) {
//        GLOBAL_FUNC_CALLS = 0;
//        double ccresult = 0.0;
//        for(const auto& p : clenshawCurtisRule(n)) {
//            ccresult += f(p.first) * p.second;
//        }
//        cc_ofs << GLOBAL_FUNC_CALLS << "\t" << error(ccresult) << "\n";
//        std::cout << "CC result: " << ccresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
//        
//        GLOBAL_FUNC_CALLS = 0;
//        double tresult = 0.0;
//        for(const auto& p : trapezoidalRule(n)) {
//            tresult += f(p.first) * p.second;
//        }
//        trap_ofs << GLOBAL_FUNC_CALLS << "\t" << error(tresult) << "\n";
//        std::cout << "T result: " << tresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
//        
//        GLOBAL_FUNC_CALLS = 0;
//        double gresult = 0.0;
//        for(const auto& p : gaussLegendreRule(n)) {
//            gresult += f(p.first) * p.second;
//        }
//        gauss_ofs << GLOBAL_FUNC_CALLS << "\t" << error(gresult) << "\n";
//        std::cout << "G result: " << gresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
//    }
    
    // Adaptive testings
    double tol = 1.0e-2;
    while(tol > 1.0e-11) {
        GLOBAL_FUNC_CALLS = 0;
        const double ccresult = adaptive1DClenshaw(-1.0, 1.0, Functor(), tol, 1000);
        cc_ofs << GLOBAL_FUNC_CALLS << "\t" << error(ccresult) << "\n";
        std::cout << "CC result: " << ccresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
        
        GLOBAL_FUNC_CALLS = 0;
        const double tresult = adaptive1DTrapezoidal(-1.0, 1.0, Functor(), tol, 1000);
        trap_ofs << GLOBAL_FUNC_CALLS << "\t" << error(tresult) << "\n";
        std::cout << "Trapezoidal result: " << tresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
        
        GLOBAL_FUNC_CALLS = 0;
        const double gresult = adaptive1DGauss(-1.0, 1.0, Functor(), tol, 1000);
        gauss_ofs << GLOBAL_FUNC_CALLS << "\t" << error(gresult) << "\n";
        std::cout << "Gauss result: " << gresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
        tol*= 0.1;
    }
    //    }
    
//    while(tol > 1.0e-7) {
//        const double tresult = adaptive1DTrapezoidal(-1.0, 1.0, Functor(), tol, 50);
//        trap_ofs << GLOBAL_FUNC_CALLS << "\t" << error(tresult) << "\n";
//        std::cout << "Trapezoidal result: " << tresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
//        GLOBAL_FUNC_CALLS = 0;
//        
//        const double gresult = adaptiveGaussLegendre(Functor(), tol, 1000);
//        gauss_ofs << GLOBAL_FUNC_CALLS << "\t" << error(gresult) << "\n";
//        std::cout << "Gauss result: " << gresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
//        GLOBAL_FUNC_CALLS = 0;
//        
//        const double ccresult = adaptive1DClenshaw(-1.0, 1.0, Functor(), tol);
//        cc_ofs << GLOBAL_FUNC_CALLS << "\t" << error(ccresult) << "\n";
//        std::cout << "CC result: " << ccresult << " with " << GLOBAL_FUNC_CALLS << " function calls\n";
//        
//        tol *= 0.1;
//    }
    std::cout << "done\n";
    
    
//    std::cout << std::setprecision(15) << adaptive1DTrapezoidal(-1.0, 1.0, Functor(), 1.0e-6, 30,2) << "\n";
//    std::cout << "func calls = " << GLOBAL_FUNC_CALLS << "\n";
//    
//    GLOBAL_FUNC_CALLS = 0;
//    std::cout << "GL converged to " << adaptiveGaussLegendre(Functor(), 1.0e-6, 1000) << " with " << GLOBAL_FUNC_CALLS <<
//    " func calls\n";
//    
//    GLOBAL_FUNC_CALLS = 0;
//    std::cout << std::setprecision(15) << adaptive1DClenshaw(-1.0, 1.0, Functor(), 1.0e-6) << "\n";
//    std::cout << GLOBAL_FUNC_CALLS << "\n";
    
	return 0;
}
