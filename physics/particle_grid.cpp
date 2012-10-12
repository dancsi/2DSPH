#include "particle_grid.h"

namespace sph
{
	particle_grid::particle_grid( double _cell_size ) : 
		cell_width(_cell_size), 
		min_corner(math::vec(DBL_MAX, DBL_MAX)), 
		max_corner(math::vec(DBL_MIN, DBL_MIN)),
		n_cells(0), n_cells_x(0), n_cells_y(0)
	{
		cells.reserve(initial_capacity);
	}

	void particle_grid::fill( particle* particles, size_t n_particles )
	{
		if(simulator::detailed_logging) 
			logger::log("filling with %d particles", n_particles);
		for(size_t i=0;i<n_cells;i++)
		{
			if(cells[i]!=nullptr)
			{
				delete cells[i];
				cells[i]=nullptr;
			}
		}
		particle_grid::n_particles=n_particles;
		set_bounds(particles, n_particles);
		n_cells_x=size_t((max_corner.x-min_corner.x+cell_width)/cell_width);
		n_cells_y=size_t((max_corner.y-min_corner.y+cell_width)/cell_width);
		n_cells=n_cells_x*n_cells_y;
		if(simulator::detailed_logging)
			logger::log("calculated bounds -> %dx%d grid.\n\tcorners: (%lf, %lf), (%lf, %lf), cell width: %.2lf", n_cells_x, n_cells_y, min_corner.x, min_corner.y, max_corner.x, max_corner.y, cell_width);
		cells.resize(n_cells);

		for(int i=0;i<n_particles;i++)
		{
			size_t idx=get_cell_index(particles[i].pos);
			if(idx==invalid_index) continue;
			if(cells[idx]==nullptr)
			{
				cells[idx]=new cell();
				cells[idx]->reserve(initial_capacity);
			}
			cells[idx]->push_back(&particles[i]);
		}
	}

	void particle_grid::set_bounds( particle* particles, size_t n_particles )
	{
		min_corner.Set(DBL_MAX, DBL_MAX);
		max_corner.Set(DBL_MIN, DBL_MIN);

		for(int i=0;i<n_particles;i++)
		{
			math::vec pos=particles[i].pos;
			if(pos.x<min_corner.x) min_corner.x=pos.x;
			if(pos.y<min_corner.y) min_corner.y=pos.y;
			if(pos.x>max_corner.x) max_corner.x=pos.x;
			if(pos.y>max_corner.y) max_corner.y=pos.y;
		}
	}

	size_t particle_grid::get_cell_index( math::vec particle_pos )
	{
		size_t idx_x=size_t((particle_pos.x-min_corner.x)/cell_width);
		size_t idx_y=size_t((particle_pos.y-min_corner.y)/cell_width);
		return get_cell_index(idx_x, idx_y);
	}

	size_t particle_grid::get_cell_index( size_t idx_x, size_t idx_y )
	{
		if(idx_x>=n_cells_x||idx_y>=n_cells_y)
			return invalid_index;
		return n_cells_y*idx_x+idx_y;
	}
}