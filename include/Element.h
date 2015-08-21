#ifndef BEMQUAD_ELEMENT_H
#define BEMQUAD_ELEMENT_H

namespace bemquad 
{
	///
	/// An element interface for boundary element integration
	///
	
	class Element 
	{
		public:

		/// Number of non-zero basis functions over this element
		virtual uint basisFuncN() const = 0;

		/// Get normal at given local parent coordinate
		virtual Point3d normal(const double u, const double v) const = 0;

		/// Get jacobian determinant at given local parent coordinate
		virtual double jacDet(const double u, const double v) const = 0;

		/// Return physical coordinate for given local parent coordinate
		virtual Point3d eval(const double u, const double v) const = 0;

		/// Evaluate the basis function for local index and local parent coordinate
		virtual double evalBasis(const uint i, const double u, const double v) const = 0;

		/// Evaluate all basis functions over this element for given local parent coordinate
		virtual std::vector<double> evalBasis(const double u, const double v) const = 0;
		
		protected:

		private:
		
	};

	///
	/// A dummy element class for testing
	///

	class DummyElement : public Element 
	{
		public:

		/// Number of non-zero basis functions over this element
		uint basisFuncN() const{ return 0;}

		/// Get normal at given local parent coordinate
		Point3d normal(const double u, const double v) const 
		{
			return Point3d(0.0, 0.0, 0.0);
		}

		/// Get jacobian determinant at given local parent coordinate
		double jacDet(const double u, const double v) const 
		{
			return 0.0;
		}

		/// Return physical coordinate for given local parent coordinate
		Point3d eval(const double u, const double v) const
		{
			return Point3d(0.0, 0.0, 0.0);
		}

		/// Evaluate the basis function for local index and local parent coordinate
		double evalBasis(const uint i, const double u, const double v) const
		{
			return 0.0;
		}

		/// Evaluate all basis functions over this element for given local parent coordinate
		std::vector<double> evalBasis(const double u, const double v) const
		{
			return {0.0};	
		}
		
	};

	
	
	
}

#endif
