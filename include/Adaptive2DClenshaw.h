//
//  Header.h
//  subvibro
//
//  Created by Robert Simpson on 14/08/2015.
//
//

#ifndef ADAPTIVE_2D_CLENSHAW_H
#define ADAPTIVE_2D_CLENSHAW_H

#include <map>

#include "common.h"
#include "IQuadrature.h"

namespace bemquad {
    
    /// A representation of a two-dimensional interval defined over
    /// the integration domain.
    
    struct Cell2D {
      
        /// Constructor
        Cell2D(const double ulower = 0.0,
               const double uupper = 0.0,
               const double vlower = 0.0,
               const double vupper = 0.0,
               const uint l = 0,
               const double v = 0.0)
        :
        lowerU(ulower),
        upperU(uupper),
        lowerV(vlower),
        upperV(vupper),
        level(l),
        value(v) {}
        
        /// Lower limit of u-coordinate
        double lowerU;
        
        /// Upper limit of u-coordinate
        double upperU;
        
        /// Lower limit of v-coordinate
        double lowerV;
        
        /// Upper limit of v-coordinate
        double upperV;
        
        /// The cell level in the hierarchy of quadrature levels
        uint level;
        
        /// The current approximation of the integrand over the cell.
        double value;
        
    };
    
    /// Subdivide a 2D cell uniformly into four cells of equal area
    std::vector<Cell2D> subdivide(const Cell2D& cell)
    {
        assert(cell.level > 0);
        const double midu = 0.5 * (cell.upperU + cell.lowerU);
        const double midv = 0.5 * (cell.upperV + cell.lowerV);
        Cell2D sw(cell.lowerU, midu, cell.lowerV, midv, cell.level - 1);
        Cell2D se(midu, cell.upperU, cell.lowerV, midv, cell.level -1);
        Cell2D nw(cell.lowerU, midu, midv, cell.upperV, cell.level - 1);
        Cell2D ne(midu, cell.upperU, midv, cell.upperV, cell.level - 1);
        return {sw, se, nw, ne};
    }
    
    /// Calculate the ratio of cell areas according to numer / denom
    double ratio(const Cell2D& numer, const Cell2D& denom)
    {
        return ((numer.upperU - numer.lowerU) * (numer.upperV - numer.lowerV)) / ((denom.upperU - denom.lowerU) * (denom.upperV - denom.lowerV));
    }
    
    /// This class is capable of integrating a given function to a prescribed
    /// accuracy using the nested Clenshaw-Curtis (CC) quadrature rule. It is
    /// based on a tensor product of 1D CC rules.
    ///
    /// It is designed primarily for non-singular and nearly singular
    /// BE integrals. It has particulary advantages when integrating
    /// oscillatory integrands over conventional a GL approach.
    ///
    /// For strongly singular integrands, the polar transformation version
    /// of this class should be used (AdaptivePolarClenshaw2DIntegrator).
    
    template<typename T, typename F>
    class AdaptiveClenshaw2DIntegrator {
        
    public:
        
        /// Constructor. Limits of integration are [a,b] x [c,d]
        AdaptiveClenshaw2DIntegrator(const double a,
                                     const double b,
                                     const double c,
                                     const double d,
                                     const F& fun,
                                     const double tol,
                                     const uint nmax,
                                     const uint nmin,
                                     const double sratio = 5.0)
        :
        mUInterval(std::make_pair(a,b)),
        mVInterval(std::make_pair(c,d)),
        mFunctor(fun),
        mTol(tol),
        mNmax(nmax),
        mNmin(nmin),
        mSubDivideRatio(sratio)
        {
            if(a > b)
                throw std::runtime_error("Invalid interval: a must be less than b");
            if(c > d)
                throw std::runtime_error("Invalid interval: c must be less than d");
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
        
    protected:
        
    private:
        
        /// Static map of quadrature rules. Allows for efficient caching.
        static std::map<uint, DPairVec> sQRuleCache;
        
        /// Get the clenshaw curtis rule for a given order
        static DPairVec clenshawRule(const uint order)
        {
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
        void globalRecursiveEval(const Cell2D& cell);
        
        /// Recursively evaluate the cell until local convergence is acheived.
        void localRecursiveEval(const Cell2D& cell);
        
        /// Get the function value at the given quadrature point. Uses cache.
        T evalFunc(const double u, const double v)
        {
            auto find = mFuncCacheMap.find(std::make_pair(u,v));
            if(find != mFuncCacheMap.end()) {
                //std::cout << "found function value: " << find->second << " at " << x << "\n";
                return find->second;
            }
            else {
                //std::cout << "inserting function value at x = " << x << "\n";
                auto i = mFuncCacheMap.insert(std::make_pair(std::make_pair(u,v), mFunctor(u,v)));
                return i.first->second;
            }
        }
        
        /// Evaluate the integral over the given interval with CC.
        /// Cached quadrature points are automatically used.
        T eval(const uint order,
               const double a,
               const double b,
               const double c,
               const double d);
        
        /// Evaluate the cell, filling in the appropriate integral value.
        void evalCell(Cell2D& cell)
        {
            cell.value = eval(cell.level, cell.lowerU, cell.upperU, cell.lowerV, cell.upperV);
        }
        
        /// Lower limit, u-direction
        double lowerU() { return mUInterval.first; }
        
        /// Lower limit, v-direction
        double lowerV() { return mVInterval.first; }
        
        /// Upper limit, u-direction
        double upperU() { return mUInterval.second; }
        
        /// Upper limit, v-direction
        double upperV() { return mVInterval.second; }
        
        /// Previous result setter
        void setPreviousResult(const T& r) { mPreviousResult = r; }
        
        /// Previous result getter
        double previousResult() const { return mPreviousResult; }
        
        /// Current result setter
        void setCurrentResult(const T& r) { mCurrentResult = r; }
        
        /// Current result setter
        double currentResult() const { return mCurrentResult; }
        
        /// u-coordinate interval
        std::pair<double, double> mUInterval;
        
        /// v-coordinate interval
        std::pair<double, double> mVInterval;
        
        /// The function we are integrating.
        F mFunctor;
        
        /// The tolerance we are aiming for with the quadrature scheme
        const double mTol;
        
        /// Maximum number of levels of quadrature
        const uint mNmax;
        
        /// Minimum number of levels of quadrature
        const uint mNmin;
        
        /// Ratio of quadrature errors that dictates when subdivision occurs
        const double mSubDivideRatio;
        
        /// Cache of function evaluations
        std::map<std::pair<double, double>, T> mFuncCacheMap;
        
        /// The previous integral approximation
        T mPreviousResult;
        
        /// The current integral approximation
        T mCurrentResult;
        
    };
    
    /// Define the static cache map for storing quadrature rules.
    template<typename T, typename F>
    std::map<uint, DPairVec> AdaptiveClenshaw2DIntegrator<T,F>::sQRuleCache;
    
    template<typename T, typename F>
    T AdaptiveClenshaw2DIntegrator<T,F>::eval()
    {
        const uint base_level = 0;
        Cell2D base_cell(lowerU(), upperU(), lowerV(), upperV(), base_level);
        evalCell(base_cell);
        setCurrentResult(base_cell.value);
        globalRecursiveEval(base_cell);
        return currentResult();
    }
    
    template<typename T, typename F>
    void AdaptiveClenshaw2DIntegrator<T,F>::globalRecursiveEval(const Cell2D& cell)
    {
        if(cell.level > maxLevelN())
        {
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
        }
        setPreviousResult(currentResult());
        Cell2D next_cell = cell; // first copy, then modify
        next_cell.level += 1;
        evalCell(next_cell);
        setCurrentResult(currentResult() - cell.value + next_cell.value);
        
        // check for global convergence
        if(std::abs(currentResult() - previousResult()) / std::abs(previousResult()) < tolerance()
           &&
           next_cell.level >= minLevelN()) {
             std::cout << "Converged at cell level " << next_cell.level << "\n";
            return;
        }
        else if(next_cell.level < 3){ // If current cell level < 3, we cannot subdivide.
            globalRecursiveEval(next_cell);
            return;
        }
        else {
            // check for local errors of subcells and subdivide
            // as appropriate.
            auto coarse_vec = subdivide(cell);
            for(auto& c : coarse_vec)
                evalCell(c);
            auto fine_vec = subdivide(next_cell);
            for(auto& c : fine_vec)
                evalCell(c);
            assert(coarse_vec.size() == fine_vec.size());
            
            // now check for error differences between subcells
            double min_error = std::numeric_limits<double>::max();
            std::vector<double> error_vec;
            double finecell_sum = 0.0;
            for(uint c = 0; c < 4; ++c) {
                const double error = std::abs(coarse_vec[c].value - fine_vec[c].value) / std::abs(coarse_vec[c].value);
                if(error < min_error)
                    min_error = error;
                error_vec.push_back(error);
                finecell_sum += fine_vec[c].value;
            }
            
            // generate vectors of cell for local/global refinement
            std::vector<uint> localrefine, globalrefine;
            for(uint c = 0; c < 4; ++c) {
                if(error_vec[c] > subDivideRatio() * min_error)
                    globalrefine.push_back(c);
                else
                    localrefine.push_back(c);
            }
            if(globalrefine.size() > 0) {
                setCurrentResult(currentResult() - next_cell.value + finecell_sum);
                for(const auto& i : localrefine)
                    localRecursiveEval(fine_vec[i]);
                for(const auto& i : globalrefine) {
                    std::cout << "local refinement of cell index: " << i << "\n";
                    globalRecursiveEval(fine_vec[i]);
                }
                return;
            }
            else {
                std::cout << "no local refinement applied. Carrying on...\n";
                globalRecursiveEval(next_cell);
                return;
            }
        }
    }
    
    template<typename T, typename F>
    void AdaptiveClenshaw2DIntegrator<T,F>::localRecursiveEval(const Cell2D& cell)
    {
        // the idea is that we recursively refine and only terminate when 'local' convergence
        // is achieved.
        if(cell.level > maxLevelN())
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
        
        setPreviousResult(currentResult());
        Cell2D next_cell = cell;
        next_cell.level += 1;
        evalCell(next_cell);
        setCurrentResult(currentResult() - cell.value + next_cell.value);
        
        // now check for local convergence
        const double error_w = ratio(next_cell, Cell2D(lowerU(), upperU(), lowerV(), upperV()));
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
    T AdaptiveClenshaw2DIntegrator<T,F>::eval(const uint order,
                                              const double a,
                                              const double b,
                                              const double c,
                                              const double d)
    {
        auto urule = clenshawRule(order,a,b);
        auto vrule = clenshawRule(order,c,d);
        T sum  = 0.0;
        for(const auto& upt : urule)
            for(const auto& vpt : vrule)
                sum += evalFunc(upt.first, vpt.first) * upt.second * vpt.second;
        return sum;
    }
    
    template< typename F >
    auto adaptive2DClenshaw(const double a,
                            const double b,
                            const double c,
                            const double d,
                            F fun,
                            const double tol = 1.0e-3,
                            const uint nmax = 20,
                            const uint nmin = 1) -> decltype(fun(a,c))
    {
        return AdaptiveClenshaw2DIntegrator<decltype(fun(a,c)), F>(a,b,c,d,fun,tol,nmax,nmin).eval();
    }
    
    
    
}

#endif
