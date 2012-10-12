#pragma once

#include "graphics.h"

namespace sph
{
	struct particle
	{
		math::vec pos, v, forces, color_field_gradient;
		double density, color_field_laplacian;
		graphics::color col;
		const static graphics::color default_color;
		void draw();
		void dump();
	};

	struct particle_pair
	{
		particle *p1, *p2;
		double dist;
		particle_pair(particle* _p1=nullptr, particle* _p2=nullptr, double _dist=0): p1(_p1), p2(_p2), dist(_dist) {}
	};
}