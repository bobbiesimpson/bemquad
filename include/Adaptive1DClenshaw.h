#ifndef ADAPTIVE_1D_CLENSHAW_H
#define ADAPTIVE_1D_CLENSHAW_H

#include "IQuadrature.h"

#include <mutex>

namespace bemquad {
    
    /// A struct which hold all the necessary information for a 'cell' that represents
    /// a range of our integration domain.
    
    struct Cell {
        
        /// Constructor
        Cell(const uint l = 0,
             const double lower = 0.0,
             const double upper = 0.0,
             const double val = 0.0)
        :
        level(l),
        lowerL(lower),
        upperL(upper),
        value(val){}

        /// The level were are on.
        uint level;
        
        /// Lower range of cell interval
        double lowerL;
        
        /// Upper range of cell interval
        double upperL;
        
        /// approximation of integrand over this cell
        double value;
        
    };
    
    std::pair<Cell, Cell> subdivide(const Cell& cell)
    {
        assert(cell.level > 0);
        Cell left, right;
        const double mid = 0.5 * (cell.lowerL + cell.upperL);
        left.lowerL = cell.lowerL; left.upperL = mid;
        right.lowerL = mid; right.upperL = cell.upperL;
        left.level = cell.level - 1; right.level = cell.level - 1;
        return std::make_pair(left, right);
    }

    /// A class to evaluate a 1D integral adaptively using Clenshaw Curtis (CC) quadrature.
    /// CC quadrature exhibits nested quadrature points that make it advantageous over
    /// G-L quadrature for adaptive methods.  CC quadrature also strikes a compromise
    /// between accuracy for integrating polynomials and an ability to reuse function
    /// evaluations for adaptivity. For oscillatory integrands, an adaptive trapezoidal rule
    /// should be preferred.
    
    template<typename T, typename F>
    class AdaptiveClenshaw1DIntegrator {
    
    public:
        
        /// Constructor
        AdaptiveClenshaw1DIntegrator(const double a,
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
        
        /// Static mutex member
        static std::mutex sMutex;
        
        /// Static map of quadrature rules. Allows for efficient caching.
        static std::map<uint, DPairVec> sQRuleCache;
        
        /// Get the clenshaw curtis rule for a given order
        static DPairVec clenshawRule(const uint order)
        {
            std::lock_guard<std::mutex> lock(sMutex);
            uint neworder = (order < 3) ? order : 3.0 * std::pow(2.0, order - 3);
            auto it = sQRuleCache.find(neworder);
            if(it != sQRuleCache.end())
                return it->second;
            else {
                //std::cout << "Generating CC rule with " << neworder + 1 << " points\n";

                auto i = sQRuleCache.insert(std::make_pair(neworder, clenshawCurtisRule(neworder + 1)));
                if(false == i.second)
                    throw std::runtime_error("Could not insert quadrature rule");
                return i.first->second;
            }
        }
        
        /// Same as above except scaled to interval [a,b]
        static DPairVec clenshawRule(const uint order, const double a, const double b)
        { return scalePts(clenshawRule(order), a, b); }
        
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
    std::map<uint, DPairVec> AdaptiveClenshaw1DIntegrator<T,F>::sQRuleCache;
    
    /// Define the static mutex
    template<typename T, typename F>
    std::mutex AdaptiveClenshaw1DIntegrator<T,F>::sMutex;
    
    
    template<typename T, typename F>
    T AdaptiveClenshaw1DIntegrator<T,F>::eval()
    {
        Cell base_cell(0,lowerL(), upperL());
        evalCell(base_cell);
        setCurrentResult(base_cell.value);
        globalRecursiveEval(base_cell);
        return currentResult();
    }
    
    template<typename T, typename F>
    void AdaptiveClenshaw1DIntegrator<T,F>::globalRecursiveEval(const Cell& cell)
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
           // std::cout << "Converged at cell level " << next_cell.level << "\n";
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
    void AdaptiveClenshaw1DIntegrator<T,F>::localRecursiveEval(const Cell& cell)
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
    
    template<typename T, typename F>
    T AdaptiveClenshaw1DIntegrator<T,F>::eval(const uint order, const double a, const double b)
    {
        auto rule = clenshawRule(order,a,b);
        T sum  = 0.0;
        for(const auto& p : rule)
            sum += evalFunc(p.first) * p.second;
        return sum;
    }
    
    
    /// Wrapper function for adaptive Clenshaw curtis rule
    template< typename F >
    auto adaptive1DClenshaw(const double a,
                            const double b,
                            F fun,
                            const double tol = 1.0e-3,
                            const uint nmax = 20,
                            const uint nmin = 1) -> decltype(fun(a))
    {
        return AdaptiveClenshaw1DIntegrator<decltype(fun(a)), F>(a,b,fun,tol,nmax,nmin).eval();
    }
    
}

#endif