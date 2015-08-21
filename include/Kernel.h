#ifndef BENQUAD_KERNEL_H
#define BENQUAD_KERNEL_H

#include <complex>

#include "Element.h"
#include "Point3d.h"
#include "common.h"
#include "MaterialParams.h"

namespace bemquad 
{
	///
	/// The class which represent a kernel for boundary element analysis
	///

	template<typename T>
	class Kernel 
	{
		public:

		/// Virtual destructor
		virtual ~Kernel() = default;
		
		/// Typedef for return type
		typedef T returnType;

		/// Evaluate the kernel
		returnType eval(const Element& el,
						const double uf,
						const double vf,
						const Point3d& spt) const
		{
			return evalImpl(el, spt);
		}
		

		/// Evaluate the kernel given local parent coordinates
		returnType eval(const Element& el,
						const double uf,
						const double vf,
						const double us,
						const double vs) const
		{
			return evalImpl(el, uf, vf, el.eval(us,vs));
		}
		
		private:

		/// Kernel evaluation implementation
		virtual returnType evalImpl(const Element& el,
									const double uf,
									const double vf,
									const Point3d& spt) const = 0;
	};

	///
	/// Helmholtz kernel
	///
	
	class HelmholtzSingleLayerKernel : public Kernel<std::complex<double>>
	{
		public:

		/// Constructor
		HelmholtzSingleLayerKernel()
		:
		Kernel(),
		mWavenumber(MaterialParam::Instance().getParam("k")) {}
		
		private:

		returnType evalImpl(const Element& el,
							const double uf,
							const double vf,
							const Point3d& spt) const
		{
			const Point3d xfield = el.eval(uf,vf);
			const double r = dist(el.eval(uf,vf), spt);
			return std::exp(std::complex<double>(0.0, -wavenumber() * r)) / (4.0 * PI * r);
		}

		/// Wavenumber getter
		double wavenumber() const { return mWavenumber; }

		/// Wavenumber 
		const double mWavenumber;
		
	};
	
	
	
}

#endif
