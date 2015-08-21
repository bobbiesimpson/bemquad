#ifndef ADAPTIVE_1D_TRAP_H
#define ADAPTIVE_1D_TRAP_H

#include <utility>
#include <stdexcept>
#include <map>
#include <algorithm>

#include "common.h"
#include "IQuadrature.h"

namespace bemquad 
{

	/// Evaluate the trapezoidal rule given a step size and vector of function values.
	/// We assume this vector is ordered in ascending order of evaluation points.
	
	template<typename T>
		T trapezoidalRule(const double h, const std::vector<T>& fvals)
	{
		T sum = 0.0;
        if(1 == fvals.size()) {
            sum += h * fvals[0];
        }
        else {
            for(uint i = 0; i < fvals.size(); ++i) {
                if(0 == i || fvals.size() - 1 == i)
                    sum += fvals[i] * h * 0.5;
                else
                    sum += fvals[i]* h;
            }
        }
		return sum;
	}
	
	/// A class that encapsulates the necessary data and functions to integrate a given
	/// function to a prescribed accuracy using adaptive techniques with the trapezoidal
	/// rule. It is parameterised by T, the return type and F, the function object which
	/// we are integrating.

	template< typename T, typename F >
	class AdaptiveTrap1DIntegrator 
	{
		private:

		/// A TCell represents a interval over which we apply the trapezoidal rule
		/// It is a simple struct that contains the step size, global function
		/// indices and the value of the trapezoidal rule applied over the cell.
		
		class TCell 
		{
			public:
			
			/// constructor
			TCell(const double hval = 0.0,
				  const std::vector<uint>& gvec = {},
				  const uint l = 0,
				  const double val = 0.0,
                  const double lower = 0.0,
				  const double upper = 0.0)
				: h(hval),
				indices(gvec),
				level(l),
				value(val),
				lowerL(lower),
				upperL(upper) {}

				/// Step size for this cell
				double h;

				/// Global function indices 
				std::vector<uint> indices;

				/// What level of trapezoidal rule is this cell. (0,1,2,...)
				uint level;

				/// approximation of integrand over this cell
				double value;

				/// Lower range of cell interval
				double lowerL;

				/// Upper range of cell interval
				double upperL;
				
		};
		
		public:

		/// Constructor
		AdaptiveTrap1DIntegrator(const double a,
								 const double b,
								 const F& fun,
								 const double tol,
								 const uint nmax,
                                 const uint nmin,
                                 const double sratio = 10.0)
        :
        mInterval(std::make_pair(a,b)),
        mFunctor(fun),
        mTol(tol),
        mNmax(nmax),
        mNmin(nmin),
        mSubDivideRatio(sratio),
        mCurrentGFuncIndex(0),
        mPreviousResult(0.0),
        mCurrentResult(0.0) {}

		/// Try and integrate until converged.
		T eval();

		/// Number of maximum levels allowed
		uint maxLevelN() const { return mNmax;}
        
        /// Minimum number of levels allowed
        uint minLevelN() const { return mNmin; }

        /// Subdivide ratio getter
        double subDivideRatio() const {return mSubDivideRatio;}
        
		/// tolerance getter
		double tolerance() const { return mTol; }

		void setFunctor(const F& f)
		{
			reset(); // if we change the functor, reset all data
			mFunctor = f;
		}
		
		private:

		/// Reset data
		void reset()
		{
			setPreviousResult(0.0);
			setCurrentResult(0.0);
			funcMap().clear();
			setCurrentGIndex(0);
		}

		/// Get the function value given the global index
		T funcEval(const uint gindex) const;

		/// evalute according to trapezoidal rule with given step size and
		/// vector of global function indices
		T eval(const double h, const std::vector<uint>& gvec) const;
		

		/// Evaluate the functino and return its global index and value
		std::pair<uint, T> evalFunc(const double x)
		{
			auto p = mFuncEvalMap.insert(std::make_pair(mCurrentGFuncIndex, mFunctor(x)));
			return std::make_pair(mCurrentGFuncIndex++, p.first->second); // remember to increment global index
		}

		/// Recursive refinement of cell until we 'globally' converged.
		void globalRecursiveEval(const TCell& cell);
        
        /// Recursively refine until 'locally' converged
        void localRecursiveEval(const TCell& cell);

        /// Returned the 'p-'refined cell with integral value computed.
        TCell refineCell(const TCell& cell);
        
		/// Subdivide given cell into subcells with equal length.
		/// does not compute integral value
		std::pair<TCell, TCell> subdivide(const TCell& cell) const;
		
		/// upper limit
		double upperL() const { return mInterval.second; }

		/// lower limit
		double lowerL() const { return mInterval.first; }

		/// Global function index setter
		void setCurrentGIndex(const uint i) { mCurrentGFuncIndex = i; }
		
		/// Function map getter
		std::map<uint, T>& funcMap() { return mFuncEvalMap;}

		/// Previous result setter
		void setPreviousResult(const T& r) { mPreviousResult = r; }

		/// Previous result getter
		double previousResult() const { return mPreviousResult; }
		
		/// Current result setter
		void setCurrentResult(const T& r) {mCurrentResult = r; }

		/// Increment current result
		void incrementCurrentResult(const T& r) { mCurrentResult += r;}

		/// Current result setter
		double currentResult() const { return mCurrentResult; }
		
		/// the interval we are integrating over
		const std::pair<double, double> mInterval;
		
		/// The functino we are integrating. Must define the ()(double) operator
		F mFunctor;

		/// The prescribed tolerance that decides convergence
		const double mTol;
		
		/// Max number of levels.
		const uint mNmax;
        
        /// Min number of levels
        const uint mNmin;
        
        /// Ratio of subcell errors that dictates when subdivision occurs.
        const double mSubDivideRatio;

		/// The index number assigned to the next function evaluation
		uint mCurrentGFuncIndex;
		
		/// Mapping from global function evaluation index to its value
		std::map<uint, T> mFuncEvalMap;

		/// The previous approximation to the integral
		T mPreviousResult;
		
		/// The current approximation to the integral
		T mCurrentResult;
		
	};

	template<typename T, typename F>
		T AdaptiveTrap1DIntegrator<T,F>::funcEval(const uint gindex) const
	{ return mFuncEvalMap.at(gindex); }
	
	template<typename T, typename F>
		T AdaptiveTrap1DIntegrator<T,F>::eval(const double h, const std::vector<uint>& gvec) const
	{
		std::vector<T> fvals;
		for(const auto& index : gvec)
			fvals.push_back(funcEval(index));
		return trapezoidalRule(h, fvals);
	}

	template<typename T, typename F>
		T AdaptiveTrap1DIntegrator<T,F>::eval() 
	{
		setCurrentResult(0.0);
		setPreviousResult(0.0);
        
        // start at level 1
        TCell cell;
        cell.level = 1;
        cell.h = upperL() - lowerL();
        cell.lowerL = lowerL(); cell.upperL = upperL();
        auto p1 = evalFunc(cell.lowerL);
        auto p2 = evalFunc(cell.upperL);
        cell.indices = {p1.first, p2.first};
        cell.value = 0.5 * cell.h * (p1.second + p2.second);
        setCurrentResult(cell.value);

		globalRecursiveEval(cell);
		return currentResult();
		
	}

	template<typename T, typename F>
		void AdaptiveTrap1DIntegrator<T,F>::globalRecursiveEval(const TCell& cell)
	{
        if(cell.level > maxLevelN()) {
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
            //std::cout << "Reached max. level. Returning result...\n";
            //return;
        }
			

		// this function computes the trapezoidal rule for the next level (i.e. higher order) of
		// the trapezoidal rule. If convergence is obtained, we return the converged result.
		// Otherwise we move to the next level.

		setPreviousResult(currentResult());
        TCell next_cell = refineCell(cell);
		setCurrentResult(currentResult() - cell.value + next_cell.value);
        
        //std::cout << "current result = " << currentResult() << "\n";
		
		// check for global convergence
       if(std::abs(currentResult() - previousResult()) / std::abs(previousResult()) < tolerance()
           &&
           cell.level > minLevelN()) {
            //std::cout << "Converged at cell level " << next_cell.level << "\n";
            return;
        }
		else if(next_cell.level < 3) { // If current cell level < 3, we cannot subdivide.
			globalRecursiveEval(next_cell);
			return;
		}
		else {
			auto coarse_cells = subdivide(cell);
			const auto& left_coarse = coarse_cells.first; const auto& right_coarse = coarse_cells.second;
			auto fine_cells = subdivide(next_cell);
			const auto& left_fine = fine_cells.first; const auto& right_fine = fine_cells.second;
            const double left_error = std::abs(left_coarse.value - left_fine.value) / std::abs(left_coarse.value);
            const double right_error = std::abs(right_coarse.value - right_fine.value) / std::abs(right_coarse.value);
            //std::cout << "left subcell error: " << left_error << "\n";
            //std::cout << "right subcell error: " << right_error << "\n";
			if(right_error > subDivideRatio() * left_error) {
                //std::cout << "subdividing right subcell at level: " << cell.level << "\n";
                localRecursiveEval(left_fine);
				globalRecursiveEval(right_fine);
				return;
			}
			else if(left_error > subDivideRatio() * right_error) {
                //std::cout << "subdividing left subcell at level: " << cell.level << "\n";
                localRecursiveEval(right_fine);
                globalRecursiveEval(left_fine);
				return;
			}
			else {
				globalRecursiveEval(next_cell);
				return;
			}
		}
	}
    
    template<typename T, typename F>
    void AdaptiveTrap1DIntegrator<T,F>::localRecursiveEval(const TCell& cell)
    {
        // the idea is that we recursively refine and only terminate when 'local' convergence
        // is achieved.
        if(cell.level > maxLevelN()) {
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
            //std::cout << "Reached max. no. of levels. REturning result...\n";
            //return;
        }
        
        
        setPreviousResult(currentResult());
        TCell next_cell = refineCell(cell);
        setCurrentResult(currentResult() - cell.value + next_cell.value);
        
        // now check for local convergence
        const double error_w = (cell.upperL - cell.lowerL) / (upperL() - lowerL());
        if(std::abs(next_cell.value - cell.value) / std::abs(cell.value) < error_w * tolerance()) {
            //std::cout << "local convergence achieved\n";
            return;
        }
        // otherwise keep raising the order of this subelement until locally converged.
        // Subdivision is NOT performed.
        else {
            localRecursiveEval(next_cell);
            return;
        }
    }
    
    template<typename T, typename F>
    typename AdaptiveTrap1DIntegrator<T,F>::TCell AdaptiveTrap1DIntegrator<T,F>::refineCell(const TCell& cell)
    {
        TCell next_cell;
        next_cell.upperL = cell.upperL; next_cell.lowerL = cell.lowerL;
        next_cell.level = cell.level + 1;
        const double h = 0.5 * cell.h;
        next_cell.h = h;
        std::vector<uint> new_indices; // index of new function evaluation indices
        double del = next_cell.lowerL + h;
        double temp_val = 0.0;
        while(del < next_cell.upperL) {
            auto p = evalFunc(del);
            new_indices.push_back(p.first);
            temp_val += p.second * h;
            del += 2.0 * h;
        }
        auto new_i = new_indices.begin();
        for(const auto& g : cell.indices) {
            next_cell.indices.push_back(g);
            if(new_i != new_indices.end())
                next_cell.indices.push_back(*new_i++);
        }
        next_cell.value = temp_val + 0.5 * cell.value;
        return next_cell;
    }

	template<typename T, typename F>
		std::pair< typename AdaptiveTrap1DIntegrator<T,F>::TCell, typename AdaptiveTrap1DIntegrator<T,F>::TCell >
		AdaptiveTrap1DIntegrator<T,F>::subdivide(const TCell& cell) const
	{
		TCell left, right; // split cell into equal sized subcells, left and right

		// set the cell step size and level
		left.h = right.h = cell.h;
		left.level = right.level = cell.level;

		// set the function indices
		const uint n = cell.indices.size();
		const uint mid_i = n / 2; // making use of floor integer division
		for(uint i = 0; i < mid_i; ++i) 
			left.indices.push_back(cell.indices[i]);
		left.indices.push_back(cell.indices[mid_i]);
		right.indices.push_back(cell.indices[mid_i]);
		for(uint i = mid_i + 1; i < n; ++i)
			right.indices.push_back(cell.indices[i]);

		// set the integral values
		left.value = eval(left.h, left.indices);
		right.value = eval(right.h, right.indices);

		// set limits
		const double mid = 0.5 * (cell.upperL + cell.lowerL);
		left.lowerL = cell.lowerL; left.upperL = mid;
		right.lowerL = mid; right.upperL = cell.upperL;
		
		return std::make_pair(left, right);
	}
	
	
	/// Wrapper function for adaptive trapezoidal rule
	template< typename F >
		auto adaptive1DTrapezoidal(const double a, const double b, F fun,
								   const double tol = 1.0e-3,
								   const uint nmax = 20,
                                   const uint nmin = 1) -> decltype(fun(a))
	{
		return AdaptiveTrap1DIntegrator<decltype(fun(a)), F>(a,b,fun,tol,nmax,nmin).eval();
	}
    
    /// Assume function lies in interval [-1,1]
    template<typename F>
    auto adaptiveGaussLegendre(F fun,
                               const double tol = 1.0e-3,
                               const uint nmax = 20,
                               const uint nmin = 1) ->decltype(fun(0.0))
    {
        typedef decltype(fun(0.0)) T;
        T previous = 0.0;
        T current = 0.0;
        
        uint currentlevel = 1;
        while(true) {
            current = 0.0;
            if(currentlevel > nmax) {
                throw std::runtime_error("Maximum number of levels reached");
            }
            auto rule = gaussLegendreRule(currentlevel);
            for(const auto& p : rule) {
                current += fun(p.first) * p.second;
            }
            if(std::abs(current - previous) / std::abs(previous) < 0.1 * tol) {
                break;
            }
            currentlevel += 1;
            previous = current;
        }
        return current;
    }
    
    
}

#endif
