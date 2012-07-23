#pragma once 

#include "vector.h"

namespace math
{
	const double scale_factor=500;
	inline double pow(double x, short n)
	{
		switch(n)
		{
		case 2: return x*x;
		case 3: return x*x*x;
		case 5: return pow(pow(x, 2), 2)*x;
		case 6: pow(pow(x, 3), 2);
		case 9: return pow(pow(x, 3), 3);
		default:
			double r=1;
			for(int i=0;i<n;i++)
				r*=x;
			return r;
		}
	}
	double w_poly6( math::vec& r, double h )
	{
		return 997.90149318618375527090119634567*pow((pow(h, 2)-r.LengthSq()), 3)/pow(h, 9);
	}

	math::vec w_gradient_spiky( math::vec& r, double h ) 
	{
		double rr=r.Length();
		return -r*14.323944878270580219199538703526*pow(h-rr, 3)/(pow(h, 6)*rr);
	}

	double w_laplacian_viscosity( math::vec& r, double h ) 
	{
		return 14.323944878270580219199538703526*(1.0-r.Length()/h)/pow(h, 5);
	}

	math::vec w_gradient_poly6( math::vec& r, double h ) 
	{
		return 9.4000888263650682688496972741891*(-r)*pow(pow(h, 2)-r.LengthSq(), 2)/pow(h, 9);
	}

	double w_laplacian_poly6( math::vec& r, double h ) 
	{
		double rr=r.LengthSq();
		return 37.600355305460273075398789096757/pow(h, 9)*(pow(h, 2)-rr)*(rr-0.75*(pow(h, 2)-rr));
	}




}