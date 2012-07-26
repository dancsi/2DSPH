#pragma once

#include <vector>
#include "math/vector.h"

namespace sph
{
	struct particle
	{
		math::vec pos, v, forces, color_field_gradient;
		double density, color_field_laplacian;
		inline void draw();
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
		int n;
		fluid(int particles_x, int particles_y, double spacing_x, double spacing_y, math::vec pos );
		~fluid();
		void prepare_step();
		void step(double dt);
		void enforce_walls(particle* p, double dt, math::vec& newpos);
		void enforce_glass(particle* p, double dt, math::vec& newpos);
		void draw();
	};
	extern std::vector<fluid> fluids;
	void init();
	int add_fluid(fluid& f);
	void step(double dt);
	void draw();
};