#pragma once

#include "particle.h"

namespace sph
{
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
}