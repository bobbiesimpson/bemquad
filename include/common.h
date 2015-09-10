#ifndef BEMQUAD_COMMON_H
#define BEMQUAD_COMMON_H

#include <vector>
#include <cmath>
#include <utility>

namespace bemquad 
{
	typedef unsigned int uint; //< Unsigned integer typedef
    
    typedef std::pair<double, double> DPair;
    
	typedef std::vector<std::vector<double>> MatrixD; // Matrix of doubles
    
	typedef std::vector<std::pair<double, double>> DPairVec;

	/// Pi, of course
	const double PI = std::atan(1.0) * 4.0;

	/// Default number of quadrature points
	const uint DEFAULT_QUADRATURE = 4;

    /// Enumeration of adaptivity for quadrature routines.
    enum class AdaptivityType {
        LOCAL = 0,
        GLOBAL
    };
	/// Enumeration of local parent coordinates
	enum class ParentDir
	{
		U = 0, V
	};
    
    const double DEFAULT_TOLERANCE = 1.0e-10;
    
    template< class T >
	bool essentiallyEqual(T a, T b, T epsilon)
	{
		return fabs( a - b ) <= ( ( fabs( a ) > fabs( b ) ? fabs( b ) : fabs( a ) ) * epsilon );
	}
    
    /// Trim a value to tval if |v| > tval and v == tval
    template< class T >
    void trimToParentInterval(T& a, T& b, const T& epsilon)
    {
        if(essentiallyEqual(std::abs(a), 1.0, epsilon))
            a = a / std::abs(a);
        if(essentiallyEqual(std::abs(b), 1.0, epsilon))
            b = b / std::abs(b);
    }
    
	
}

#endif
