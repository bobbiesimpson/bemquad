#ifndef BEMQUAD_COMMON_H
#define BEMQUAD_COMMON_H

#include <vector>
#include <cmath>
#include <utility>

namespace bemquad 
{
	typedef unsigned int uint; //< Unsigned integer typedef
	typedef std::vector<std::vector<double>> MatrixD; // Matrix of doubles
	typedef std::vector<std::pair<double, double>> DPairVec;

	/// Pi, of course
	const double PI = std::atan(1.0) * 4.0;

	/// Default number of quadrature points
	const uint DEFAULT_QUADRATURE = 4;

	/// Enumeration of local parent coordinates
	enum class ParentDir
	{
		U = 0, V
	};
	
}

#endif
