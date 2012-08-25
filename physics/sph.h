#pragma once

#include <vector>
#include "math/vector.h"
#include <boost/function.hpp>
#include "graphics.h"

namespace sph
{
	struct particle
	{
		math::vec pos, v, forces, color_field_gradient;
		double density, color_field_laplacian;
		graphics::color col;
		const static graphics::color default_color;
		inline void draw();
		inline void dump();
	};

	struct particle_pair
	{
		particle *p1, *p2;
		double dist;
		particle_pair(particle* _p1=nullptr, particle* _p2=nullptr, double _dist=0): p1(_p1), p2(_p2), dist(_dist) {}
	};

	typedef std::vector<particle*> cell;

	struct particle_grid
	{
		const static size_t initial_side_length=50,  
			initial_capacity=initial_side_length*initial_side_length;
		size_t n_cells, n_cells_x, n_cells_y, n_particles;
		double cell_width;
		math::vec min_corner, max_corner;
		std::vector<cell*> cells;
		particle_grid(double _cell_size);
		void fill(particle* particles, size_t n_particles);
		void set_bounds(particle* particles, size_t n_particles);
		size_t get_cell_index(math::vec particle_pos);
		size_t get_cell_index(size_t idx_x, size_t idx_y);
		static const size_t invalid_index=-1u;

	};

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