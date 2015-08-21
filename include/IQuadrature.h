#ifndef BEMQUAD_IQUADRATURE_H
#define BEMQUAD_IQUADRATURE_H

#include "common.h"
#include <cstddef>
#include <stdexcept>
#include <cassert>

namespace bemquad 
{

	/// Function for computing Gauss-Legendre points and weights
	DPairVec gaussLegendreRule(const uint npts);

	/// Rule for computing Clenshaw Curtis rule
	DPairVec clenshawCurtisRule(const uint npts);

	/// Trapezoidal rule
	DPairVec trapezoidalRule(const uint npts);
    
    /// Scale the quadrature points and weights to the interval [a,b]
    DPairVec scalePts(const DPairVec& pvec, const double a, const double b);
	
	///
	/// Simple struct for holding a two-dimensional quadrature point
	///
	
	struct QPt
	{
		/// Default constructor
		QPt(const double uin = 0.0,
			const double vin = 0.0) : u(uin), v(vin) {}

		/// U-coordinate
		double u;

		/// V-coordinate
		double v;
	};

	///
	/// Enumeration of available quadrature rules
	///

	enum class QRule 
	{
		GAUSS_LEGENDRE,
		CLENSHAW_CURTIS,
		TRAPEZOIDAL
	};
	
	///
	/// An abstract class that represents an iterator for quadrature over a
	/// two-dimensional manifold.
	///
	
	class IQuadrature 
	{
		public:

		/// Construct with equal points in each parametric direction
		IQuadrature(const uint npts = DEFAULT_QUADRATURE,
					const QRule rule = QRule::GAUSS_LEGENDRE) :
		IQuadrature(npts, npts, rule) {}

		/// Construct with different numbers of points in each parametric direction
		IQuadrature(const uint nu,
					const uint nv,
					const QRule rule) :
		mCurrentI(0),
		mPointN(nu * nv)
		{
			init(nu,nv,rule);
		}

		/// Const quadrature getter
		const DPairVec& quadrature(const ParentDir dir) const
		{
			if(ParentDir::U == dir)
				return mQuadratureU;
			else {
				assert(ParentDir::V == dir);
				return mQuadratureV;
			}
		}

		/// Get the current index
		uint currentIndex() const { return mCurrentI; }
		
		/// Reset the quadrature index
		void reset() { mCurrentI = 0; }

		/// Check if we've finished iterating
		bool isDone() const { return mCurrentI == mPointN; }

		/// Quadrature point getter
		QPt get() const;

		/// Weight getter
		double getWeight() const;

		/// Prefix increment
		IQuadrature& operator++()
		{
			++mCurrentI;
			return *this;
		}

		private:
		
		/// Initialise quadrature
		void init(const uint nu,
				  const uint nv,
				  const QRule rule)
		{
			DPairVec& u_quad = quadrature(ParentDir::U);
			DPairVec& v_quad = quadrature(ParentDir::V);
			switch(rule) {
				case QRule::GAUSS_LEGENDRE:
					u_quad = gaussLegendreRule(nu);
					if(nu == nv)
						v_quad = u_quad;
					else
						v_quad = gaussLegendreRule(nv);
					break;
				case QRule::CLENSHAW_CURTIS:
					u_quad = clenshawCurtisRule(nu);
					if(nu == nv)
						v_quad = u_quad;
					else
						v_quad = clenshawCurtisRule(nv);
					break;
				case QRule::TRAPEZOIDAL:
					u_quad = trapezoidalRule(nu);
					if(nu == nv)
						v_quad = u_quad;
					else
						v_quad = trapezoidalRule(nv);
					break;
				default:
					throw std::runtime_error("Bad quadrature rule specified.");
					break;
			}
		}

		/// Non-const quadrature getter
		DPairVec& quadrature(const ParentDir dir)
		{
			if(ParentDir::U == dir)
				return mQuadratureU;
			else {
				assert(ParentDir::V == dir);
				return mQuadratureV;
			}
		}
		
		/// Clear all stored points and weights
		void clear()
		{
			quadrature(ParentDir::U).clear();
			quadrature(ParentDir::V).clear();
			mPointN = 0;
		}

		
		/// Get number of quadrature points in specified direction
		std::size_t pointN(const ParentDir dir) const
		{
			switch(dir) {
				case ParentDir::U:
					return mQuadratureU.size();
				case ParentDir::V:
					return mQuadratureV.size();
				default:
					throw std::runtime_error("Bad direction specified in pointN()");
					return 0;
			}
		}

		/// Get quadrature point given direction and local index
		double currentPt(const ParentDir dir, uint lindex) const
		{
			return quadrature(dir).at(lindex).first;
		}
		
		/// Get quadrature weight given direction and local index
		double currentWeight(const ParentDir dir, uint lindex) const
		{
			return quadrature(dir).at(lindex).second;
		}

		/// Function to map a current global quadrature index to local indices
		/// in U- and V- direction
		std::pair<uint, uint> currentLocalIndices() const
		{
			const uint index = currentIndex();
			return std::make_pair(index / pointN(ParentDir::U),
								  index % pointN(ParentDir::U));
		}

		/// Quadrature points and weights - U direction
		std::vector<std::pair<double, double>> mQuadratureU;

		/// Quadrature points and weight - V direction
		std::vector<std::pair<double, double>> mQuadratureV;

		/// Current index
		uint mCurrentI;

		/// Total number of quadrature points
		uint mPointN;
		
	};
	
}

#endif
