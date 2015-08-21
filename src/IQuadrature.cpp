#include <cassert>
#include <stdexcept>

#include "IQuadrature.h"

namespace bemquad 
{
	DPairVec gaussLegendreRule(const uint npts) 
	{
		DPairVec rvec;
		if( npts == 0 )
			rvec = {};
		else if(npts == 1)
			rvec.push_back(std::make_pair(0.0,2.0));
		else if(npts == 2) {
			rvec.push_back(std::make_pair(-1 / sqrt( 3 ), 1.0));
			rvec.push_back(std::make_pair(1 / sqrt( 3 ), 1.0));
		}
		else if(npts == 3) {
			rvec.push_back(std::make_pair(-sqrt( 3.0 / 5.0 ), 5.0 / 9.0));
			rvec.push_back(std::make_pair(0.0, 8.0 / 9.0));
			rvec.push_back(std::make_pair(sqrt( 3.0 / 5.0 ), 5.0 / 9.0));
		}
		else if( npts == 4 ) {
			const double w1 =  ( 18.0 - sqrt( 30.0 ) ) / 36.0;
			const double w2 =  ( 18.0 + sqrt( 30.0 ) ) / 36.0;
                
			double inner = sqrt( ( 3.0 - 2.0 * sqrt( 6.0 / 5.0 ) ) / 7.0 );
			double outer = sqrt( ( 3.0 + 2.0 * sqrt( 6.0 / 5.0 ) ) / 7.0 );
			
			rvec.push_back(std::make_pair(-outer, w1));
			rvec.push_back(std::make_pair(-inner, w2));
			rvec.push_back(std::make_pair(inner, w2));
			rvec.push_back(std::make_pair(outer, w1));
		}
		else if( npts == 5 ) {
			const double w1 = (322.0 - 13.0*sqrt(70.0))/900.0;
			const double w2 = (322.0 + 13.0*sqrt(70.0))/900.0;
			const double w3 = 128.0/225.0;
                
			double inner = 1.0/3.0 * sqrt( 5.0 - 2.0*sqrt(10.0/7.0));
			double outer = 1.0/3.0 * sqrt( 5.0 + 2.0*sqrt(10.0/7.0));
                
			rvec.push_back(std::make_pair(-outer, w1));
			rvec.push_back(std::make_pair(-inner, w2));
			rvec.push_back(std::make_pair(0.0, w3));
			rvec.push_back(std::make_pair(inner, w2));
			rvec.push_back(std::make_pair(outer, w1));
		}
		else if( npts == 6 ) {
			const double w1 = 0.171324492379170;
			const double w2 = 0.360761573048139;
			const double w3 = 0.467913934572691;
                
			double inner = 0.238619186083197;
			double mid = 0.661209386466265;
			double outer = 0.932469514203152;

			rvec.push_back(std::make_pair(-outer, w1));
			rvec.push_back(std::make_pair(-mid, w2));
			rvec.push_back(std::make_pair(-inner, w3));
			rvec.push_back(std::make_pair(inner, w3));
			rvec.push_back(std::make_pair(mid, w2));
			rvec.push_back(std::make_pair(outer, w1));
		}
		else if(npts == 7) {
			const double w1 = 0.129484966168870;
			const double w2 = 0.279705391489277;
			const double w3 = 0.381830050505119;
			const double w4 = 0.417959183673469;
                
			double inner = 0.405845151377397;
			double mid = 0.741531185599394;
			double outer = 0.949107912342759;

			rvec.push_back(std::make_pair(-outer, w1));
			rvec.push_back(std::make_pair(-mid, w2));
			rvec.push_back(std::make_pair(-inner, w3));
			rvec.push_back(std::make_pair(0.0, w4));
			rvec.push_back(std::make_pair(inner, w3));
			rvec.push_back(std::make_pair(mid, w2));
			rvec.push_back(std::make_pair(outer, w1));
		}
		else if( npts == 8 ) {
			const double w1 = 0.101228536290376;
			const double w2 = 0.222381034453374;
			const double w3 = 0.313706645877887;
			const double w4 = 0.362683783378362;
                
			const double p1 = 0.960289856497536; 
			const double p2 = 0.796666477413627; 
			const double p3 = 0.525532409916329; 
			const double p4 = 0.183434642495650; 

			rvec.push_back(std::make_pair(-p1, w1));
			rvec.push_back(std::make_pair(-p2, w2));
			rvec.push_back(std::make_pair(-p3, w3));
			rvec.push_back(std::make_pair(-p4, w4));
			rvec.push_back(std::make_pair(p4, w4));
			rvec.push_back(std::make_pair(p3, w3));
			rvec.push_back(std::make_pair(p2, w2));
			rvec.push_back(std::make_pair(p1, w1));
		}
		else if( npts == 9 ) {
			const double w1 = 0.081274388361574;
			const double w2 = 0.180648160694857;
			const double w3 = 0.260610696402935;
			const double w4 = 0.312347077040003;
			const double w5 = 0.330239355001260;
                
			const double p1 =  0.968160239507626; 
			const double p2 = 0.836031107326636; 
			const double p3 = 0.613371432700590; 
			const double p4 = 0.324253423403809; 

			rvec.push_back(std::make_pair(-p1, w1));
			rvec.push_back(std::make_pair(-p2, w2));
			rvec.push_back(std::make_pair(-p3, w3));
			rvec.push_back(std::make_pair(-p4, w4));
			rvec.push_back(std::make_pair(0.0, w5));
			rvec.push_back(std::make_pair(p4, w4));
			rvec.push_back(std::make_pair(p3, w3));
			rvec.push_back(std::make_pair(p2, w2));
			rvec.push_back(std::make_pair(p1, w1));
		}
		else if( npts == 10 )
		{       
			const double w1 = 0.066671344308688;
			const double w2 = 0.149451349150581;
			const double w3 = 0.219086362515982;
			const double w4 = 0.269266719309996;
			const double w5 = 0.295524224714753;
                
			const double p1 = 0.973906528517172; 
			const double p2 = 0.865063366688985; 
			const double p3 = 0.679409568299024; 
			const double p4 = 0.433395394129247; 
			const double p5 = 0.148874338981631;

			rvec.push_back(std::make_pair(-p1, w1));
			rvec.push_back(std::make_pair(-p2, w2));
			rvec.push_back(std::make_pair(-p3, w3));
			rvec.push_back(std::make_pair(-p4, w4));
			rvec.push_back(std::make_pair(-p5, w5));
			rvec.push_back(std::make_pair(p5, w5));
			rvec.push_back(std::make_pair(p4, w4));
			rvec.push_back(std::make_pair(p3, w3));
			rvec.push_back(std::make_pair(p2, w2));
			rvec.push_back(std::make_pair(p1, w1));
		}
		else
		{
			// numerically determine gauss points and weights
			// slow but nice to have for higher order gauss rules.
			std::vector<double> pts, wts;
			pts.resize( npts, 0.0 );
			wts.resize( npts, 0.0 );
                
			const double x1 = -1.0;
			const double x2 = 1.0;
			const double EPS = 1.0e-14;
			int m = ( npts + 1 ) / 2;
			double xm = 0.5 * ( x2 + x1 );
			double xl = 0.5 * ( x2 - x1 );
			for( int i = 0; i < m; ++i )
			{
				double z = cos( PI * ( i + 0.75 ) / ( npts + 0.5 ) );
				double z1, pp, p3;
				do
				{
					double p1( 1.0 );
					double p2( 0.0 );
					for( uint j = 0; j < npts; ++j )
					{
						p3 = p2;
						p2 = p1;
						p1 = ( ( 2.0 * j + 1.0 ) * z * p2 - j * p3 ) / ( j + 1 );
					}
					pp = npts * ( z * p1 - p2 ) / ( z * z - 1.0 );
					z1 = z;
					z = z1 - p1 / pp;
				} while ( fabs( z - z1 ) > EPS );
				pts[ i ] = xm - xl * z;
				pts[ npts - 1 - i ] = xm + xl * z;
				wts[ i ] = 2.0 * xl / ( ( 1.0 - z * z ) * pp * pp );
				wts[ npts - 1 - i ] = wts[ i ];
			}
			for(uint i = 0; i < npts; ++i) 
				rvec.push_back(std::make_pair(pts[i], wts[i]));
			
		}
		return rvec;
	}

	DPairVec clenshawCurtisRule(const uint npts) 
	{
		assert(npts >=1);

		std::vector<double> x(npts);
		std::vector<double> w(npts);
		double b;

		if(npts == 1) {
			x[0] = 0.0;
			w[0] = 2.0;
		}
		else {
			for(uint i = 0; i < npts; i++) {
				x[i] =  std::cos ( ( double ) ( npts - 1 - i ) * PI
								   / ( double ) ( npts - 1 ) );
			}
			x[0] = -1.0;
			if ( ( npts % 2 ) == 1 ) 
				x[(npts-1)/2] = 0.0;
			x[npts-1] = +1.0;

			for(uint  i = 0; i < npts; i++) {
				const double theta = ( double ) ( i ) * PI / ( double ) ( npts - 1 );
				w[i] = 1.0;
				for (uint j = 1; j <= ( npts - 1 ) / 2; j++) {
					if ( 2 * j == ( npts - 1 ) )
						b = 1.0;
					else
						b = 2.0;
					w[i] = w[i] - b *  std::cos ( 2.0 * ( double ) ( j ) * theta )
						/ ( double ) ( 4 * j * j - 1 );
				}
			}
			w[0] = w[0] / ( double ) ( npts - 1 );
			for(uint i = 1; i < npts - 1; i++) 
				w[i] = 2.0 * w[i] / ( double ) ( npts - 1 );
			w[npts-1] = w[npts-1] / ( double ) ( npts - 1 );
		}
		DPairVec rvec;
		for(uint i = 0; i < npts; ++i) 
			rvec.push_back(std::make_pair(x[i], w[i]));
		return rvec;
	}

	DPairVec trapezoidalRule(const uint npts) 
	{
		uint n;
		if(npts < 2) {
			throw std::runtime_error("Cannot use trapezoidal rule with npts < 2. Attempting to use npts = 2");
			n = 2;
		}
		else
			n = npts;

		const double h = 2.0 / (n - 1.0);
		DPairVec rvec;
		for(uint i = 0; i < n; ++i) 
			rvec.push_back(std::make_pair(h * i - 1.0, h));
		rvec.front().second *= 0.5; // modify weights of first and last points
		rvec.back().second *= 0.5;
		return rvec;
	}
    
    DPairVec scalePts(const DPairVec& pvec, const double a, const double b)
    {
        assert(b > a);
        DPairVec newpvec = pvec;
        const double c1 = 0.5 * (b - a);
        const double c2 = 0.5 * (a + b);
        for(auto& p : newpvec) {
            p.first = c1 * p.first + c2;
            p.second = p.second * c1;
        }
        return newpvec;
    }
	
    double scalePt(const double xi, const double a, const double b)
    {
        return xi * (b - a) * 0.5 + 0.5 * (a + b);
    }
    
    DPair scalePt(const std::pair<double, double>& xi,
                  const std::pair<double, double>& ulimits,
                  const std::pair<double, double>& vlimits)
    {
        return std::make_pair(scalePt(xi.first, ulimits.first, ulimits.second), scalePt(xi.second, vlimits.first, vlimits.second));
    }
	
    std::pair<double, double> scaleToParentInterval(const std::pair<double, double>& p,
                                                    const std::pair<double, double>& ulimits,
                                                    const std::pair<double, double>& vlimits)
    {
        const double urange = ulimits.second - ulimits.first;
        const double umid = 0.5 * (ulimits.second + ulimits.first);
        const double vrange = vlimits.second - vlimits.first;
        const double vmid = 0.5 * (vlimits.second + vlimits.first);
        return std::make_pair(2.0 / urange * (p.first - umid),
                              2.0 / vrange * (p.second - vmid));
    }
    
	QPt IQuadrature::get() const
	{
		auto lpair = currentLocalIndices();
		return QPt(currentPt(ParentDir::U, lpair.first),
				   currentPt(ParentDir::V, lpair.second));
	}
	
	double IQuadrature::getWeight() const
	{
		auto lpair = currentLocalIndices();
		return currentWeight(ParentDir::U, lpair.first) *
			currentWeight(ParentDir::V, lpair.second);
	}
}
