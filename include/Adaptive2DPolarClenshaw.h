//
//  Adapative2DPolarClenshaw.h
//  subvibro
//
//  Created by Robert Simpson on 18/08/2015.
//
//

#ifndef subvibro_Adaptive2DPolarClenshaw_h
#define subvibro_Adaptive2DPolarClenshaw_h

#include "Adaptive2DClenshaw.h"
#include "common.h"
#include "IQuadrature.h"

namespace bemquad {
    
    /// Enumeration of polar integration domains
    
    enum class PolarDir {
        N = 0, E, S, W
    };
    
    /// A functor class that applies the appropriate transformation
    /// for polar integration. It is assumed that the variables lie
    /// in the standard parent interval [-1,1].
    
    template<typename F>
    class PolarFunctor {
        
    private:
        
        /// Limits of integration domain, u-direction
        const std::pair<double, double> mULimits;
        
        /// Limits of integration domain, v-direction
        const std::pair<double, double> mVLimits;
        
        /// Source point coordinates (u,v) coordinate system
        const std::pair<double, double> mGlbSrc;
        
        /// Source points coordinates in local parent coordinate system
        const std::pair<double, double> mLocalSrc;
        
        /// Domain of theta coordinate
        std::pair<double, double> mThetaLimits;
        
        /// Non-const Functor getter
        F& functor() { return mFunc; }
        
        /// The original functor
        F mFunc;
        
        /// The polar integration domain enumeration
        const PolarDir mDir;
        
        /// Determinant of jacobian (xi,eta) -> (u,v)
        const double mDetJacob;
        
        /// Tolerance for shifting quadrature points that lie directly
        /// on the source point coordinate.
        const double mTol;
        
    public:
        
        /// Constructor
        PolarFunctor(const std::pair<double, double>& ulimits,
                     const std::pair<double, double>& vlimits,
                     const std::pair<double, double>& scoords,
                     const F& fun,
                     const PolarDir dir,
                     const double tol = 1.0e-7)
        :
        mULimits(ulimits),
        mVLimits(vlimits),
        mGlbSrc(scoords),
        mLocalSrc(scaleToParentInterval(scoords, ulimits, vlimits)),
        mFunc(fun),
        mDir(dir),
        mDetJacob(0.25 * (ulimits.second - ulimits.first) * (vlimits.second - vlimits.first)),
        mTol(tol)
        {
            init();
        }
        
        /// Polar direction getter
        const PolarDir& polarDir() const { return mDir; }
        
        /// Functor getter
        //const F& functor() const { return mFunc; }
        
        const DPair& ulimits() const { return mULimits; }
        
        const DPair& vlimits() const { return mVLimits; }
        
        const DPair& glbSrcPt() const { return mGlbSrc; }
        
        /// Given a parent coordinate from the integrator, return the function value
        /// (along with any jacobian determinants) of the transformed integrand.
        auto operator()(const double xi, const double eta) -> decltype(functor()(xi,eta))
        {
            double xi_n = xi;
            double eta_n = eta;
            
            // If we are on the edge eta = -1 we have to deal with rho = 0 in which
            // the integrand will be undefined.
            if(essentiallyEqual(eta_n, -1.0, tolerance()))
                eta_n = eta_n + tolerance();
            
            const double theta = 0.5 * (upperTheta() - lowerTheta()) * xi_n + 0.5 * (upperTheta() + lowerTheta());
            double rhohat;
            
            switch (polarDir()) {
                case PolarDir::E:
                    rhohat = (1.0 - localSrcPtXi()) / std::cos(theta);
                    break;
                case PolarDir::N:
                    rhohat = (1.0 - localSrcPtEta()) / std::sin(theta);
                    break;
                case PolarDir::W:
                    rhohat = -(1.0 + localSrcPtXi()) / std::cos(theta);
                    break;
                case PolarDir::S:
                    rhohat = -(1.0 + localSrcPtEta()) / std::sin(theta);
                    break;
            }
            const double rho = 0.5 * (eta_n + 1.0) * rhohat;
            const double polarjacob = 0.25 * (upperTheta() - lowerTheta()) * rhohat;
            auto parentpair = std::make_pair(localSrcPtXi() + rho * std::cos(theta),
                                                  localSrcPtEta() + rho * std::sin(theta));
            trimToParentInterval(parentpair.first, parentpair.second, tolerance());

            const auto glbpair = scalePt(parentpair, ulimits(), vlimits());
            return functor()(glbpair.first, glbpair.second) * rho * polarjacob * jacobDet();
        }
        
        /// Tolerance getter
        double tolerance() const {return mTol; }
        
    private:
        
        /// Initialise the transformation data
        void init()
        {
            // intialise the theta limits
            switch (polarDir()) {
                case PolarDir::E:
                    mThetaLimits = std::make_pair(-std::atan((1.0 + localSrcPtEta()) / (1.0 - localSrcPtXi())),
                                                  std::atan((1.0 - localSrcPtEta()) / (1.0 - localSrcPtXi())));
                    break;
                case PolarDir::N:
                    mThetaLimits = std::make_pair(PI / 2.0  - std::atan((1.0 - localSrcPtXi()) / (1.0 - localSrcPtEta())),
                                                  PI / 2.0 + std::atan((localSrcPtXi() + 1.0) / (1.0 - localSrcPtEta())));
                    break;
                case PolarDir::W:
                    mThetaLimits = std::make_pair(PI - std::atan((1.0 - localSrcPtEta()) / (localSrcPtXi() + 1.0)),
                                                  PI + std::atan((localSrcPtEta() + 1.0) / (localSrcPtXi() + 1.0)));

                    break;
                case PolarDir::S:
                    mThetaLimits = std::make_pair(3.0 * PI / 2.0 - std::atan((localSrcPtXi() + 1.0) / (localSrcPtEta() + 1.0)),
                                                  3.0 * PI / 2.0 + std::atan((1.0 - localSrcPtXi()) / (localSrcPtEta() + 1.0)));
                    break;
            }
        }
        
        /// Get the determinant of the jacobian that maps the integration domain to the parent
        /// interval [-1,1]x[-1,1]
        double jacobDet() const { return mDetJacob; }
        
        double lowerU() const { return mULimits.first; }
        
        double upperU() const { return mULimits.second; }
        
        double lowerV() const { return mVLimits.first; }
        
        double upperV() const { return mVLimits.second; }
        
        double localSrcPtXi() const { return mLocalSrc.first; }
        
        double localSrcPtEta() const { return mLocalSrc.second; }
        
        double lowerTheta() const { return mThetaLimits.first; }
        
        double upperTheta() const { return mThetaLimits.second; }
        
        /// Map a parent coordinate to a 'global' coordinate in the original integration domain
        std::pair<double, double> globalCoord(const double xi, const double eta) const
        {
            return std::make_pair(scalePt(xi, lowerU(), upperU()), scalePt(eta, lowerV(), upperV()));
        }
        

        
    };
    
    /// A class to integrate adaptively using Clenshaw-Curtis quadrature
    /// and a polar transformaiton applied at the specified
    /// local coordinate.
    
    template<typename T, typename F>
    class Adaptive2DPolarClenshawIntegrator {
      
    public:
        
        /// Constructor where an integrand over
        /// [a,b] x [c,d] is assumed. The singular point
        /// (p1,p2) is assumed to lie in this interval.
        
        Adaptive2DPolarClenshawIntegrator(const double a,
                                          const double b,
                                          const double c,
                                          const double d,
                                          const double p1,
                                          const double p2,
                                          F& fun,
                                          const double tol,
                                          const AdaptivityType adapt = AdaptivityType::LOCAL,
                                          const uint nmax = 20,
                                          const uint nmin = 1,
                                          const double sratio = 5.0)
        :
        mUInterval(std::make_pair(a,b)),
        mVInterval(std::make_pair(c,d)),
        mSrcCoords(std::make_pair(p1,p2)),
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
        
        /// Functor getter
        const F& functor() const { return mFunctor; }
        
        void setFunctor(const F& f)
        {
            //reset(); // if we change the functor, reset all data
            mFunctor = f;
        }
        
        /// Adaptivity type getter
        const AdaptivityType adaptivityType() const { return mAdaptivity; }
        
        /// Adaptivity type setter
        void setAdaptivity(const AdaptivityType t) { mAdaptivity = t; }
        
    private:
        
        const DPair& ulimits() const { return mUInterval; }
        
        const DPair& vlimits() const { return mVInterval; }
        
        const DPair& glbSrcPt() const { return mSrcCoords; }
        
        /// Lower limit, u-direction
        double lowerU() { return mUInterval.first; }
        
        /// Lower limit, v-direction
        double lowerV() { return mVInterval.first; }
        
        /// Upper limit, u-direction
        double upperU() { return mUInterval.second; }
        
        /// Upper limit, v-direction
        double upperV() { return mVInterval.second; }
        
        /// Source coordinate, u-direction
        double sourceU() { return mSrcCoords.first; }
        
        /// Source coordinate, v-direction
        double sourceV() { return mSrcCoords.second; }
        
        /// Non const functor getter
        F& functor() { return mFunctor; }
        
        /// u-coordinate interval
        std::pair<double, double> mUInterval;
        
        /// v-coordinate interval
        std::pair<double, double> mVInterval;
            
        /// Coordinates of source point
        std::pair<double, double> mSrcCoords;
        
        /// Reference to the function we are integrating.
        F& mFunctor;
        
        /// The tolerance we are aiming for with the quadrature scheme
        const double mTol;
        
        /// Flag for adaptivity type
        AdaptivityType mAdaptivity;
        
        /// Maximum number of levels of quadrature
        const uint mNmax;
        
        /// Minimum number of levels of quadrature
        const uint mNmin;
        
        /// Ratio of quadrature errors that dictates when subdivision occurs
        const double mSubDivideRatio;
        
    };
    
    template<typename T, typename F>
    T Adaptive2DPolarClenshawIntegrator<T,F>::eval()
    {
        const double tol = DEFAULT_TOLERANCE;
        T sum = 0.0;

        // E polar interval
        if(!essentiallyEqual(upperU(), sourceU(), tol)) {
            PolarFunctor<F> t1(ulimits(), vlimits(), glbSrcPt(), functor(), PolarDir::E);
            sum += adaptive2DClenshaw(-1.0, 1.0, -1.0, 1.0, t1, tolerance(), adaptivityType(), maxLevelN(), minLevelN());
        }

        //std::cout << sum << "\n";
        
        // N polar interval
        if(!essentiallyEqual(upperV(), sourceV(), tol)) {
            PolarFunctor<F> t1( ulimits(), vlimits(), glbSrcPt(),functor(), PolarDir::N);
            sum += adaptive2DClenshaw(-1.0, 1.0, -1.0, 1.0, t1, tolerance(), adaptivityType(), maxLevelN(), minLevelN());
        }
        
        //std::cout << sum << "\n";
        // W polar interval
        if(!essentiallyEqual(lowerU(), sourceU(), tol)) {
            PolarFunctor<F> t1( ulimits(), vlimits(), glbSrcPt(), functor(), PolarDir::W);
            sum += adaptive2DClenshaw(-1.0, 1.0, -1.0, 1.0, t1, tolerance(), adaptivityType(), maxLevelN(), minLevelN());
        }
        //std::cout << sum << "\n";
        
        // S polar interval
        if(!essentiallyEqual(lowerV(), sourceV(), tol)) {
            PolarFunctor<F> t1(ulimits(), vlimits(), glbSrcPt(), functor(), PolarDir::S);
            sum += adaptive2DClenshaw(-1.0, 1.0, -1.0, 1.0, t1, tolerance(), adaptivityType(), maxLevelN(), minLevelN());
        }
        //std::cout << sum << "\n";
        return sum;
    }
    
    template< typename F >
    auto adaptive2DPolarClenshaw(const double a,
                                 const double b,
                                 const double c,
                                 const double d,
                                 const double p1,
                                 const double p2,
                                 F fun,
                                 const double tol = 1.0e-3,
                                 const AdaptivityType adapt = AdaptivityType::LOCAL,
                                 const uint nmax = 20,
                                 const uint nmin = 1) -> decltype(fun(a,c))
    {
        return Adaptive2DPolarClenshawIntegrator<decltype(fun(a,c)), F>(a,b,c,d,p1,p2,fun,tol,adapt,nmax,nmin).eval();
    }
    
    
    
}

#endif
