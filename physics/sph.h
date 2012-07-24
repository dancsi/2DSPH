#pragma once

#include <vector>
#include "math/vector.h"

namespace sph
{
	struct fluid
	{
		struct particle
		{
			math::vec pos, v, f_pressure, f_viscosity, colour_field_gradient;
			double density, colour_field_laplacian;
			uint16_t hash;
			unsigned char flags;
			enum FLAGS {TRACK = 1};
			inline void draw();
			inline void track();
			inline void rehash();
		} *particles;
		std::vector<int> *hash_table;
		double eta, sigma, mass, smoothing_length, density, k, gradient_length_treshold, damping, rest_density;
		int n;
		fluid(int particles_x, int particles_y, double spacing_x, double spacing_y, math::vec pos );
		~fluid();
		void prepare_step();
		void step(double dt);
		void draw();
	};
	extern std::vector<fluid> fluids;
	void init();
	int add_fluid(fluid& f);
	void step(double dt);
	void draw();
};