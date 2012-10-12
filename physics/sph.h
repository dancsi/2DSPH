#pragma once

#include <vector>
#include "math/vector.h"
#include <boost/function.hpp>

#include "particle.h"
#include "particle_grid.h"

namespace sph
{
	struct fluid
	{
		particle *particles;
		particle_grid* grid;
		double eta, sigma, mass, smoothing_length, density, k, gradient_length_treshold, damping, rest_density, wall_damping;
		double glass_pushback_distance, glass_pushback_minimum_pressure, glass_density, glass_viscosity;
		int north_particle_idx, south_particle_idx, east_particle_idx, west_particle_idx;
		static double max_force;
		int n;
		fluid(int particles_x, int particles_y, double spacing_x, double spacing_y, math::vec pos );
		~fluid();
		void prepare_step();
		void step(double dt);
		void glass_common_update(boost::function<void (fluid*, particle*, math::vec)> f);
		//template <typename T>
		//void glass_common_update(T f);
		void calculate_densities();
		void calculate_densities(int x, int y);
		void calculate_glass_fluid_densities(particle* p, math::vec r);
		void calculate_forces();
		void calculate_forces(particle_pair p);
		void calculate_glass_fluid_forces(particle* p, math::vec r);
		void update_positions(double dt);
		void enforce_walls(particle* p, double& dt, math::vec& newpos);
		void enforce_glass(particle* p, double& dt, math::vec& newpos);
		void set_bounding_particle_indices();
		void draw();
	};
	extern std::vector<fluid> fluids;
	void init();
	int add_fluid(fluid& f);
	void step(double dt);
	void draw();
};