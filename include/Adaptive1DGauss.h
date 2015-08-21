//
//  Adaptive1DGauss.h
//  subvibro
//
//  Created by Robert Simpson on 10/08/2015.
//
//

#ifndef subvibro_Adaptive1DGauss_h
#define subvibro_Adaptive1DGauss_h

#include "Adaptive1DClenshaw.h"
#include "IQuadrature.h"

namespace bemquad {
    
    /// A class that represents adaptive quadrature using Gauss-Legendre
    /// quadrature.  In practice this class would not be used since GL
    /// quadrature is not nested and therefore very inefficienct for
    /// adaptive algorithms. Its main purpose is for testing and comparision.
    
    template<typename T, typename F>
    class AdaptiveGauss1DIntegrator {
    
    public:
        
        /// Constructor
        AdaptiveGauss1DIntegrator(const double a,
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
        mPreviousResult(0.0),
        mCurrentResult(0.0)
        {
            if(a > b)
                throw std::runtime_error("Invalid interval: a must be less than b");
        }
        
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
            //reset(); // if we change the functor, reset all data
            mFunctor = f;
        }
        
    private:
        
        /// Static map of quadrature rules. Allows for efficient caching.
        static std::map<uint, DPairVec> sQRuleCache;
        
        /// Get the clenshaw curtis rule for a given order
        static DPairVec gaussRule(const uint n)
        {
            auto it = sQRuleCache.find(n);
            if(it != sQRuleCache.end())
                return it->second;
            else {
                //std::cout << "Generating CC rule with " << neworder + 1 << " points\n";
                auto i = sQRuleCache.insert(std::make_pair(n, gaussLegendreRule(n)));
                if(false == i.second)
                    throw std::runtime_error("Could not insert quadrature rule");
                return i.first->second;
            }
        }
        
        /// Same as above except scaled to interval [a,b]
        static DPairVec gaussRule(const uint n, const double a, const double b)
        { return scalePts(gaussRule(n), a, b); }
        
        /// Recursively evaluate the cell until global convergence is acheived.
        void globalRecursiveEval(const Cell& cell);
        
        /// Recursively evaluate the cell until local convergence is acheived.
        void localRecursiveEval(const Cell& cell);
        
        /// Get the function value at the given quadrature point. Uses cache.
        T evalFunc(const double x)
        {
            auto find = mFuncCacheMap.find(x);
            if(find != mFuncCacheMap.end()) {
                //std::cout << "found function value: " << find->second << " at " << x << "\n";
                return find->second;
            }
            else {
                //std::cout << "inserting function value at x = " << x << "\n";
                auto i = mFuncCacheMap.insert(std::make_pair(x,mFunctor(x)));
                return i.first->second;
            }
        }
        
        /// Evaluate the integral over the given interval with CC.
        /// Cached quadrature points are automatically used.
        T eval(const uint order, const double a, const double b);
        
        /// Evaluate the cell, filling in the appropriate integral value.
        void evalCell(Cell& cell)
        {
            cell.value = eval(cell.level, cell.lowerL, cell.upperL);
        }
        
        /// upper limit
        double upperL() const { return mInterval.second; }
        
        /// lower limit
        double lowerL() const { return mInterval.first; }
        
        /// Previous result setter
        void setPreviousResult(const T& r) { mPreviousResult = r; }
        
        /// Previous result getter
        double previousResult() const { return mPreviousResult; }
        
        /// Current result setter
        void setCurrentResult(const T& r) { mCurrentResult = r; }
        
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
        /// I.e. if relative error in left subcell > mSubDivideRatio * right subcell error
        /// we subdivide the left subcell and vice versa.
        const double mSubDivideRatio;
        
        /// Map of quadrature point in interval [a,b] and its associated function value
        std::map<double, T> mFuncCacheMap;
        
        /// The previous approximation to the integral
        T mPreviousResult;
        
        /// The current approximation to the integral
        T mCurrentResult;
        
    };
    
    /// Define the static cache map for storing quadrature rules.
    template<typename T, typename F>
    std::map<uint, DPairVec> AdaptiveGauss1DIntegrator<T,F>::sQRuleCache;
    
    template<typename T, typename F>
    T AdaptiveGauss1DIntegrator<T,F>::eval()
    {
        Cell base_cell(1,lowerL(), upperL());
        evalCell(base_cell);
        setCurrentResult(base_cell.value);
        globalRecursiveEval(base_cell);
        return currentResult();
    }
    
    template<typename T, typename F>
    T AdaptiveGauss1DIntegrator<T,F>::eval(const uint n,
                                           const double a,
                                           const double b)
    {
        auto rule = gaussRule(n,a,b);
        T sum  = 0.0;
        for(const auto& p : rule)
            sum += evalFunc(p.first) * p.second;
        return sum;
    }
    
    
    template<typename T, typename F>
    void AdaptiveGauss1DIntegrator<T,F>::globalRecursiveEval(const Cell& cell)
    {
        if(cell.level > maxLevelN())
        {
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
        }
        setPreviousResult(currentResult());
        Cell next_cell = cell; // first copy, then modify
        next_cell.level += 1;
        evalCell(next_cell);
        setCurrentResult(currentResult() - cell.value + next_cell.value);
        
        // check for global convergence
        if(std::abs(currentResult() - previousResult()) / std::abs(previousResult()) < tolerance()
           &&
           next_cell.level >= minLevelN()) {
            //std::cout << "Converged at cell level " << next_cell.level << "\n";
            return;
        }
        else if(next_cell.level < 3){ // If current cell level < 3, we cannot subdivide.
            globalRecursiveEval(next_cell);
            return;
        }
        else {
            // check errors in subintervals
            auto coarse_pair = subdivide(cell);
            auto fine_pair = subdivide(next_cell);
            auto& lcoarse = coarse_pair.first; auto& rcoarse = coarse_pair.second;
            auto& lfine = fine_pair.first; auto& rfine = fine_pair.second;
            evalCell(lcoarse);
            evalCell(rcoarse);
            evalCell(lfine);
            evalCell(rfine);
            const double lerror = std::abs(lfine.value - lcoarse.value) / std::abs(lcoarse.value);
            const double rerror = std::abs(rfine.value - rcoarse.value) / std::abs(rcoarse.value);
            if(rerror > subDivideRatio() * lerror) {
                //std::cout << "subdividing right subcell at level: " << cell.level << "\n";
                setCurrentResult(currentResult() - next_cell.value + lfine.value + rfine.value);
                localRecursiveEval(lfine);
                globalRecursiveEval(rfine);
                return;
            }
            else if(lerror > subDivideRatio() * rerror) {
                //std::cout << "subdividing left subcell at level: " << cell.level << "\n";
                setCurrentResult(currentResult() - next_cell.value + lfine.value + rfine.value);
                localRecursiveEval(rfine);
                globalRecursiveEval(lfine);
                return;
            }
            else {
                //std::cout << "no subdivision. continuing to next level\n";
                globalRecursiveEval(next_cell);
                return;
            }
        }
    }
    
    template<typename T, typename F>
    void AdaptiveGauss1DIntegrator<T,F>::localRecursiveEval(const Cell& cell)
    {
        // the idea is that we recursively refine and only terminate when 'local' convergence
        // is achieved.
        if(cell.level > maxLevelN())
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
        
        setPreviousResult(currentResult());
        Cell next_cell = cell;
        next_cell.level += 1;
        evalCell(next_cell);
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
    
    
    /// Wrapper function for adaptive Gauss rule
    template< typename F >
    auto adaptive1DGauss(const double a,
                         const double b,
                         F fun,
                         const double tol = 1.0e-3,
                         const uint nmax = 20,
                         const uint nmin = 1) -> decltype(fun(a))
    {
        return AdaptiveGauss1DIntegrator<decltype(fun(a)), F>(a,b,fun,tol,nmax,nmin).eval();
    }
    
}


#endif
