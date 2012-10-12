#pragma once 

#include "vector.h"

namespace math
{
	const double PI=3.1415926535897932384626433832795;
	double pow(double x, short n);

	extern double 
		w_poly6_factor,
		w_poly6_gradient_factor,
		w_poly6_laplacian_factor,
		w_viscosity_laplacian_factor,
		w_spiky_gradient_factor;
	extern double h, h_squared;

	void calculate_common_factors(double _h);
	double w_poly6( double r);
	math::vec w_poly6_gradient( math::vec& r, double length_hint);
	double w_poly6_laplacian( double r);
	math::vec w_spiky_gradient( math::vec& r, double length_hint);
	double w_viscosity_laplacian( double r);

}