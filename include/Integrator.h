#ifndef BEMQUAD_INTEGRATOR_H
#define BEMQUAD_INTEGRATOR_H

#include "common.h"
#include "Point3d.h"
#include "Element.h"
#include "Kernel.h"

namespace bemquad 
{
	///
	/// A class to encapsulate adaptive quadrature routines for boundary element
	/// analysis.
	///
	/// It is paramterised by K which defines the relevant kernel.
	///

	template<typename K>
	class Integrator 
	{
		public:

		/// Typedef of kernel for this integrator
		typedef K KernelType;

		/// Typedef of the matrix type returned by the integrator. Depends on the kernelType
		typedef std::vector<std::vector<typename K::returnType>> Matrix;
		
		/// Evaluate element/source point combination for each basis function over this element
		Matrix eval(const Element& el, const Point3d& spt) const;

		/// Evaluate element/source point combination for source point defined by local parent coordinate
		Matrix eval(const Element& el,
					const double u,
					const double v) const;

		/// Evaluate the area of this element
		double evalArea(const Element& el) const;

		private:

		/// Kernel getter
		const KernelType& kernel() const {return mKernel;}
		
		/// The kernel instance
		const KernelType mKernel;
		
	};

	template<typename K>
		double Integrator<K>::evalArea(const Element& el) const
	{
		return 0.0;
	}
	
	template<typename K>
		typename Integrator<K>::Matrix Integrator<K>::eval(const Element& el,
														   const Point3d& spt) const
	{
		// The general concept is to perform adaptive quadrature until we
		// reach a converged result for the boundary integral.

		for(uint ibasis = 0; ibasis < el.basisFuncN(); ++ibasis) {
			
		}
		
		return {{0.0}}; // Return an empty matrix for now.
	}
		
	template<typename K>
		typename Integrator<K>::Matrix Integrator<K>::eval(const Element& el,
														   const double u,
														   const double v) const
	{
		return {{0.0}}; // Return an empty matrix for now.
	}

	/// Helmholtz integrator typedef
	typedef Integrator<HelmholtzSingleLayerKernel>  HelmholtzIntegrator;
	
	
}

#endif
