#include "MaterialParams.h"
#include <fstream>

namespace bemquad 
{
	std::ostream& operator<<(std::ostream& ost, const MaterialParam& m)
	{
		for(const auto& p : m.data())
			ost << "{" << p.first << ", " << p.second << "} ";
		return ost;
	}
	
}
