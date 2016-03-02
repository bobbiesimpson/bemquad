#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <complex>
#include <map>

#include "common.h"

uint GLB_FUNC_CALLS = 0;
uint GLB_TFUNC_CALLS = 0;

const double PI = 4.0 * std::atan(1.0);

std::map<double, double> fcall_map;

double f(const double x) 
{
	//return std::exp(std::cos(200.0 * x));
	//return 1.0 / (x + 100.0);
	//return x;

	GLB_FUNC_CALLS += 1;
	//return std::exp(std::cos(0.01 * x)) / (x + 100.0);
	//return 1.0 / (x + 100.0);
    //return std::exp(std::cos(20.0 * x)) * 1.0 / (20.0 + x);
    return 1.0 / std::pow((1.004 - x), 2.0);
}

double ft(const double x)
{
	GLB_TFUNC_CALLS += 1;
	//return x;
	//return std::exp( std::cos(0.01 * x)) / (x + 100.0);
	//return 1.0 / (100.0 + x);
	//return std::exp(std::cos(20.0 * x)) * 1.0 / (20.0 + x);
    return 1.0 / std::pow((1.004 - x), 2.0);
}


int main(int argc, char* argv[]) 
{
	std::cout << "Performing quadrature with double exponential formula....\n";
	std::cout << std::setprecision(15);
	const double hmax = 3.7;

	const double tol = 1.0e-7;
	
	double prevsum = 0.0;
	double sum = 0.0;

	
	const double a = -1.0; // intervals
	const double b = 1.0;
	const uint N = 20; // max number of sample points

	/// First the trapezoidal rule
	std::ofstream ofs("trapezoidal_results.out");
	if(!ofs)
		throw std::runtime_error("Cannot open file for writing");
	ofs << std::setprecision(15);
	
	double tsum = 0.0;
	double prev_tsum = 0.0;
	for(uint n = 1; n < N; ++n) {

		if(1 == n) {
			tsum += ft(0.5 * (b + a)) * (b - a);
			std::cout << tsum << "\n";
		}
		else if(2 == n) {
            tsum = 0.25 * (b - a) * (ft(a) + ft(b)) + 0.5 * prev_tsum;
			std::cout << tsum << "\n";
		}
		else {
			double tempsum = 0.0;
			const double h = (b - a) / std::pow(2.0, n -1);
			//std::cout << "h = " << h << "\n";
			double del = a + h;
			while(del < b) {
				//std::cout << del << "\n";
				tempsum += h * ft(del);
				del += 2.0 * h;
			}
			tsum = tempsum + 0.5 * prev_tsum;
			std::cout << tsum << "\n";
			if(std::abs(tsum - prev_tsum)/std::abs(prev_tsum) < tol) {
				std::cout << "Converged to result = " << tsum << " with n = " << n << " and " << GLB_TFUNC_CALLS << " functino calls\n";
				//break;
			}
		}
		ofs << GLB_TFUNC_CALLS << "\t" << tsum << "\n";
		prev_tsum = tsum;
	}

	std::cout << "DE INTEGRATION\n\n";
	std::ofstream ofs2("dbl_exp_results.out");
	if(!ofs2)
		throw std::runtime_error("Cannot open file for writing");
	ofs2 << std::setprecision(15);
	
	for(uint n = 1; n < N; ++n) {
		if(1 == n) {
			sum += hmax * 2.0 * (b - a) * 0.25 * f(0.5 * (b - a));
			prevsum = sum;
		}
		else {
			uint it = 1;
			for(uint j = 1; j < n - 1; ++j) it <<= 1;
			double th = hmax / it;
			double t = 0.5 * th;
			double tempsum = 0.0;
			for(uint j = 0; j < it; ++j) {
				double q = std::exp(-2.0 * std::sinh(t));
				double del = (b - a) * q / (1.0 + q);
				double fact = q / ( (1.0 + q) * (1.0 + q))  * std::cosh(t);
				tempsum += fact * (f(a + del) + f(b - del));
				t += th;
			}
			prevsum = sum;
			sum = 0.5 * sum + (b - a) * th * tempsum;
			std::cout << "n = " << n << " sum = " << sum << "\n";
			if(std::abs(sum - prevsum)/std::abs(prevsum) < tol) {
				std::cout << "Converged to result = " << sum << " with n = " << n << " and " << GLB_FUNC_CALLS << " functino calls\n";
				//break;
			}
		}
		ofs2 << GLB_FUNC_CALLS << "\t" << sum << "\n";
	}
	return 0;
}
