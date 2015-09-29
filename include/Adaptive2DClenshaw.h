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
#include <mutex>
#include <limits>

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
        value(v),
        parent(nullptr) {}
        
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
        
        /// Reference to parent (nil if no parent)
        Cell2D* parent;
        
        /// References to children (nil if no children)
        std::vector<Cell2D*> children;
        
    };
    
    /// Are the ranges of the two cells equal?
    bool equalRange(const Cell2D& c1, const Cell2D& c2)
    {
        return c1.lowerU == c2.lowerU && c1.upperU == c2.upperU && c1.lowerV == c2.lowerV && c1.upperV == c2.upperV;
    }
    
    /// Subdivide a 2D cell uniformly into four cells of equal area
    std::vector<Cell2D> subdivide(Cell2D& cell)
    {
        assert(cell.level > 0);
        const double midu = 0.5 * (cell.upperU + cell.lowerU);
        const double midv = 0.5 * (cell.upperV + cell.lowerV);
        Cell2D sw(cell.lowerU, midu, cell.lowerV, midv, cell.level - 1);
        Cell2D se(midu, cell.upperU, cell.lowerV, midv, cell.level -1);
        Cell2D nw(cell.lowerU, midu, midv, cell.upperV, cell.level - 1);
        Cell2D ne(midu, cell.upperU, midv, cell.upperV, cell.level - 1);
        std::vector<Cell2D> subcells{sw, se, nw, ne};
        for(auto& c : subcells) {
            c.parent = &cell;
            cell.children.push_back(&c);
        }
        return subcells;
    }
    
    /// Generate a subcell that is formed from the 'centre' of the
    /// interval of the given cell. It is used to check for singularities
    /// that may lie near the midpoint of the cell.
    Cell2D centreSubcell(const Cell2D& cell)
    {
        assert(cell.level > 0);
        const double u_h = 0.25 * (cell.upperU - cell.lowerU);
        const double v_h = 0.25 * (cell.upperV - cell.lowerV);
        return Cell2D(cell.lowerU + u_h,
                      cell.upperU - u_h,
                      cell.lowerV + v_h,
                      cell.upperV - v_h,
                      cell.level - 1);
    }
    
    /// Return the set of cell that correspond to the set intersection of the given cell
    /// and the 'centre' sub cell defined by the function above
    std::vector<Cell2D> rimCells(const Cell2D& cell)
    {
        const uint ncell = 4;
        const double u_h = 1.0 / ncell * (cell.upperU - cell.lowerU);
        const double v_h = 1.0 / ncell * (cell.upperV - cell.lowerV);

        std::vector<Cell2D> rvec;
        for(uint i = 0; i < ncell; ++i) {
            const double lower_u = cell.lowerU + i * u_h;
            const double upper_u = lower_u + u_h;
            for(uint j = 0; j < ncell; ++j) {
                if((1 == i || 2 == i) && (1 == j || 2 == j))
                    continue;
                const double lower_v = cell.lowerV + j * v_h;
                const double upper_v = lower_v + v_h;
                rvec.push_back(Cell2D(lower_u, upper_u, lower_v, upper_v, cell.level -1));
            }
        }
        return rvec;
    }
    
    /// Calculate the ratio of cell areas according to numer / denom
    double ratio(const Cell2D& numer, const Cell2D& denom)
    {
        return ((numer.upperU - numer.lowerU) * (numer.upperV - numer.lowerV))
        / ((denom.upperU - denom.lowerU) * (denom.upperV - denom.lowerV));
    }
    
    /// This class is capable of integrating a given function to a prescribed
    /// accuracy using the nested Clenshaw-Curtis (CC) quadrature rule. It is
    /// based on a tensor product of 1D CC rules.
    ///
    /// It is designed primarily for non-singular and nearly singular
    /// BE integrals. It has particulary advantages when integrating
    /// oscillatory integrands over a conventional GL approach.
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
                                     F& fun,
                                     const double tol,
                                     const AdaptivityType adapt = AdaptivityType::LOCAL,
                                     const uint nmax = 20,
                                     const uint nmin = 1,
                                     const double sratio = 5.0)
        :
        mUInterval(std::make_pair(a,b)),
        mVInterval(std::make_pair(c,d)),
        mFunctor(fun),
        mTol(tol),
        mAdaptivity(adapt),
        mNmax(nmax),
        mNmin(nmin),
        mSubDivideRatio(sratio)
        {
            if(a > b)
                throw std::runtime_error("Invalid interval: a must be less than b");
            if(c > d)
                throw std::runtime_error("Invalid interval: c must be less than d");
        }
        
        /// Clear the cached data
        void clear()
        {
            mFuncCacheMap.clear();
            mPreviousResult = 0.0;
            mCurrentResult = 0.0;
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
        
        /// Adaptivity type getter
        const AdaptivityType adaptivityType() const { return mAdaptivity; }
        
        /// Adaptivity type setter
        void setAdaptivity(const AdaptivityType t) { mAdaptivity = t; }
        
        
    protected:
        
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
        void globalRecursiveEval(Cell2D& cell);
        
        /// Recursively evaluate the cell until local convergence is acheived.
        void localRecursiveEval(Cell2D& cell);
        
        /// Given two cell of differing levels of quadrature order, perform
        /// adaptive h-refinement until convergence is achieved.
        void adaptiveRecursiveEval(Cell2D& coarse,
                                   Cell2D& fine);
        
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
        
        T relative_error(T a, T b)
        {
            if((std::abs(a) < DEFAULT_TOLERANCE) && (std::abs(b) < DEFAULT_TOLERANCE))
                return std::max(std::abs((a - b) / (1.0 + a)), std::abs((a - b) / (1.0 + b)));
            else
                return std::max(std::abs((a - b) / (a)), std::abs((a - b) / (b)));
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
        
        /// Reference to the function we are integrating.
        F& mFunctor;
        
        /// The tolerance we are aiming for with the quadrature scheme
        const double mTol;
        
        /// Local or global adaptivity flag.
        AdaptivityType mAdaptivity;
        
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
    
    template<typename T, typename F>
    std::mutex AdaptiveClenshaw2DIntegrator<T,F>::sMutex;
    
    /// Define the static cache map for storing quadrature rules.
    template<typename T, typename F>
    std::map<uint, DPairVec> AdaptiveClenshaw2DIntegrator<T,F>::sQRuleCache;
    
    template<typename T, typename F>
    T AdaptiveClenshaw2DIntegrator<T,F>::eval()
    {
        // Perform one subdivision to deal with singularity at midpoint case
        const uint base_level = 1;
        Cell2D base(lowerU(), upperU(), lowerV(), upperV(), base_level);
        auto subcells = subdivide(base);
        double integral = 0.0;
        for(auto & c : subcells) {
            evalCell(c);
            integral += c.value;
        }
        setCurrentResult(integral);
        for(auto& c : subcells)
            globalRecursiveEval(c);
        
//        const uint base_level = 0;
//        Cell2D base_cell(lowerU(), upperU(), lowerV(), upperV(), base_level);
//        evalCell(base_cell);
//        setCurrentResult(base_cell.value);
//        globalRecursiveEval(base_cell);
        
        return currentResult();
    }
    
    template<typename T, typename F>
    void AdaptiveClenshaw2DIntegrator<T,F>::globalRecursiveEval(Cell2D& cell)
    {
        if(cell.level > maxLevelN())
        {
            /// We haven't converged, so let's try integrating with a polar transformation
            /// using an approximate nearest point.
            
            /// TODO: write a function for determining the nearest point (within a certain tolerance)
            /// and then pass this approximate point to the polar version of this integrator.
            
            /// Actually, I don't want to call it here - instead, I want to call the polar transformation when
            /// the number of local subdivisions haas reached a maximum.
            
            throw std::runtime_error("Reaced maximum number of levels in adaptive quadrature.");
        }
        setPreviousResult(currentResult());
        Cell2D next_cell = cell; // first copy, then modify
        next_cell.level += 1;
        evalCell(next_cell);
        setCurrentResult(currentResult() - cell.value + next_cell.value);
        
        // check for global convergence
        if(relative_error(currentResult(), previousResult())/*(std::abs(currentResult() - previousResult()) / std::abs(previousResult()))*/ < tolerance()
           &&
           next_cell.level >= minLevelN()) {
             //std::cout << "Converged at cell level " << next_cell.level << "\n";
            //std::cout << "converged with cell ratio: " << ratio(Cell2D(lowerU(), upperU(), lowerV(), upperV()), cell) << "\n";
            return;
        }
        else if(next_cell.level < 3 || AdaptivityType::GLOBAL == adaptivityType()){ // If current cell level < 3, we cannot subdivide.
            globalRecursiveEval(next_cell);
            return;
        }
        else {
            // check for local errors of subcells and subdivide
            // as appropriate.
            adaptiveRecursiveEval(cell, next_cell);
        }
    }
    
    template<typename T, typename F>
    void AdaptiveClenshaw2DIntegrator<T,F>::adaptiveRecursiveEval(Cell2D& coarse,
                                                                  Cell2D& fine)
    {
        assert(equalRange(coarse, fine)); // make sure the ranges are equal
        auto coarse_vec = subdivide(coarse);
        for(auto& c : coarse_vec)
            evalCell(c);
        auto fine_vec = subdivide(fine);
        for(auto& c : fine_vec) {
            evalCell(c);
        }
        assert(coarse_vec.size() == fine_vec.size());
        
        // now check for error differences between subcells
        double min_error = std::numeric_limits<double>::max();
        std::vector<double> error_vec;
        double finecell_sum = 0.0;
        double coarsecell_sum = 0.0;
        for(uint c = 0; c < 4; ++c) {
            const double error = std::abs(coarse_vec[c].value - fine_vec[c].value) / std::abs(coarse_vec[c].value);
            //const double error = relative_error(coarse_vec[c].value, fine_vec[c].value);
            if(error < min_error)
                min_error = error;
            error_vec.push_back(error);
            finecell_sum += fine_vec[c].value;
            coarsecell_sum += coarse_vec[c].value;
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
            setCurrentResult(currentResult() - fine.value + finecell_sum);
            for(const auto& i : localrefine)
                localRecursiveEval(fine_vec[i]);
            for(const auto& i : globalrefine) {
                //std::cout << "local refinement of cell index: " << i << "\n";
                globalRecursiveEval(fine_vec[i]);
            }
            return;
        }
        else {
            
            // let's see if the error on the current cell is much larger than
            // the error seen in its near neighbours
            Cell2D* parent = fine.parent;
            double cell_error = 0.0;
            double min_e = std::numeric_limits<double>::max();
            for(uint i = 0; i < 4; ++i) {
                Cell2D c = *parent->children[i];
                c.level = coarse.level;
                evalCell(c);
                Cell2D f = c;
                f.level = fine.level;
                evalCell(f);
                //const double e = relative_error(c.value, f.value);
                const double e = std::abs(f.value - c.value) / std::abs(f.value);
                if(equalRange(f, fine))
                    cell_error = e;
                else {
                    if(e < min_e)
                        min_e = e;
                }
            }
            if(cell_error > subDivideRatio() * min_e) {
                setCurrentResult(currentResult() - fine.value + finecell_sum);
                for(auto& c : fine_vec)
                    globalRecursiveEval(c);
                return;
            }
            else {
                globalRecursiveEval(fine);
                return;
            }
        
            // let's try subdivsion centred around the midpoint of the integration domain
            // to cope with singularities that might lie near to the midpoint.
           /*auto centre_coarse = centreSubcell(coarse);
            evalCell(centre_coarse);
            auto centre_fine = centreSubcell(fine);
            evalCell(centre_fine);
            const double coarse_rim_val = coarse.value - centre_coarse.value;
            const double fine_rim_val = fine.value - centre_fine.value;
            
//            auto rim_cells_coarse = rimCells(coarse);
//            double coarse_rim_val = 0.0;
//            for(auto& c : rim_cells_coarse) {
//                evalCell(c);
//                coarse_rim_val += c.value;
//            }
//            auto rim_cells_fine = rimCells(fine);
//            double fine_rim_val = 0.0;
//            for(auto& c : rim_cells_fine) {
//                evalCell(c);
//                fine_rim_val += c.value;
//            }
            
            const double centre_error = std::abs(centre_fine.value - centre_coarse.value) / std::abs(centre_fine.value);
            const double rim_error = std::abs(fine_rim_val - coarse_rim_val) / std::abs(fine_rim_val);
            
            if(centre_error > subDivideRatio() * rim_error) {
                std::cout << "Subdividing from midpoint\n";
                setCurrentResult(currentResult() - fine.value + finecell_sum);
                for(auto& c : fine_vec)
                    globalRecursiveEval(c);
                return;
//
//                globalRecursiveEval(centre_fine);
//                for(auto& c : rim_cells_fine)
//                    localRecursiveEval(c);
//                return;
            }
            else*/
//                globalRecursiveEval(fine);
//                return;
//            }
        }
    }
    
    template<typename T, typename F>
    void AdaptiveClenshaw2DIntegrator<T,F>::localRecursiveEval(Cell2D& cell)
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
                            const AdaptivityType adapt = AdaptivityType::LOCAL,
                            const uint nmax = 20,
                            const uint nmin = 1) -> decltype(fun(a,c))
    {
        return AdaptiveClenshaw2DIntegrator<decltype(fun(a,c)), F>(a,b,c,d,fun,tol,adapt,nmax,nmin).eval();
    }
    
    
    
}

#endif
