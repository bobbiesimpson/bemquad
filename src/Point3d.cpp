#include <cmath>
#include "Point3d.h"

namespace bemquad 
{
	const size_t Point3d::SIZE = 3;
	
	double Point3d::get(uint i) const
	{
		assert(i < size());
		return mData[i];
	}

	double dist(const Point3d& p1, const Point3d& p2)
	{
		double temp = 0.0;
		for(uint i = 0; i < Point3d::SIZE; ++i) 
			temp += ( p2.get(i) - p1.get(i) ) * ( p2.get(i) - p1.get(i) );
		return std::sqrt(temp);
	}
}
