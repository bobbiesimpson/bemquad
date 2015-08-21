#ifndef BEMQUAD_MATERIAL_PARAM_H
#define BEMQUAD_MATERIAL_PARAM_H

#include <map>
#include <utility>
#include <algorithm>
#include <string>

namespace bemquad 
{
	///
	/// A singleton class for holding material parameters used by kernels.
	/// Examples include wavenumber, Young's modulus, Poisson's ratio.
	///
	
	class MaterialParam 
	{
		public:

		/// Get the one and only instance
		static MaterialParam& Instance()
		{
			static MaterialParam theInstance;
			return theInstance;
		}

		/// Add a material parameter with a string identifier
		bool addParam(const std::string& name,
					  const double val)
		{
			auto result = data().insert(std::make_pair(name, val));
			return result.second;
		}

		/// Get the parameter. Throws a runtime error if not found
		double getParam(const std::string& name)
		{
			auto search = data().find(name);
			if(search == data().end())
				throw std::runtime_error("Cannot find parameter with name: " + name);
			return search->second;
		}

		private:

		/// Typedef of data type
		typedef std::map<std::string, double> MapType;
		
		/// Constructor (made private)
		MaterialParam() = default;

		/// Destructor (made private)
		~MaterialParam()  = default;

		/// Copy constructor (made private)
		MaterialParam(const MaterialParam& m)  = default;

		/// Move constructor (made private)
		MaterialParam(MaterialParam&& m) = default;

		/// Assignment copy constructor (made private)
		MaterialParam& operator=(const MaterialParam& m) = default;

		/// Move assignment (made private)
		MaterialParam& operator=(MaterialParam&& m) = default;

		/// Data accessor
		MapType& data() { return mData; }

		/// Const data accessor
		const MapType& data() const { return mData; }
		
		/// Map of material parameters
		std::map<std::string, double> mData;

		/// Overload output operator
		friend std::ostream& operator<<(std::ostream& ost, const MaterialParam& p);
	};
}

#endif
