#ifndef BEMQUAD_POINT3D_H
#define BEMQUAD_POINT3D_H

#include <vector>
#include <cassert>
#include <cstddef>

#include "common.h"

namespace bemquad
{
	///
	/// Representation of a three dimensional point
	///
	
	class Point3d 
	{
		public:

		/// Static member variable that defines the number of components
		static const std::size_t SIZE;
		
		/// Default constructor
		Point3d(const double x = 0.0,
				const double y = 0.0,
				const double z = 0.0) :
		mData{x, y, z} {}

		/// Getter
		double get(uint i) const;

		/// Number of components
		uint size() const { return 3; }
		
		protected:

		private:
		
		/// The 'raw' data
		std::vector<double> mData;

	};

	/// Get distance between two points
	double dist(const Point3d& p1, const Point3d& p2);
	
}

#endif
