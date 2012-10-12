#include "kernels.h"

namespace math
{
	double pow(double x, short n)
	{
		switch(n)
		{
		case 0: return 1;
		case 1: return x;
		case 2: return x*x;
		case 3: return x*x*x;
		case 4: return (x*x)*(x*x);
		case 5: return (x*x)*(x*x)*x;
		case 6: return (x*x*x)*(x*x*x);
		case 7: return (x*x*x)*(x*x*x)*x;
		case 8: return ((x*x)*(x*x))*((x*x)*(x*x));
		case 9: return (x*x*x)*(x*x*x)*(x*x*x);
		default:
			//__debugbreak();
			double r=1;
			for(int i=0;i<n;i++)
				r*=x;
			return r;
		}
	}

	double 
		w_poly6_factor(0),
		w_poly6_gradient_factor(0),
		w_poly6_laplacian_factor(0),
		w_viscosity_laplacian_factor(0),
		w_spiky_gradient_factor(0);
	double h(0), h_squared(0);

	void calculate_common_factors(double _h)
	{
		h=_h;
		h_squared=h*h;
		w_poly6_factor=315.0/(64.0*PI*pow(h, 9));
		w_poly6_laplacian_factor=945.0/(8*PI*pow(h, 9));
		w_poly6_gradient_factor=-945.0/(32*PI*pow(h, 9));
		w_spiky_gradient_factor=-45.0/(PI*pow(h, 6));
		w_viscosity_laplacian_factor=45.0/(PI*pow(h, 5));
	}

	double w_poly6( double r)
	{
		return w_poly6_factor*pow(h_squared-r*r, 3);
	}

	math::vec w_poly6_gradient( math::vec& r, double length_hint) 
	{
		return w_poly6_gradient_factor*r*pow(h_squared-length_hint*length_hint, 2);
	}

	double w_poly6_laplacian( double r) 
	{
		return w_poly6_laplacian_factor*(h_squared-r*r)*(1.75*r*r-0.75*h_squared);
	}

	math::vec w_spiky_gradient( math::vec& r, double length_hint) 
	{
		return (w_spiky_gradient_factor/length_hint)*pow(h-length_hint, 2)*r;
	}

	double w_viscosity_laplacian( double r) 
	{
		return w_viscosity_laplacian_factor*(1.0-r/h);
	}
}