#include "sph.h"

#include "simulator.h"
#include "math/vector.h"
#include "math/kernels.h"
#include "math/line.h"
#include "graphics.h"
#include "config.h"
#include "boost/foreach.hpp"
#include "physics.h"
#include <omp.h>
#include <windows.h>


namespace sph
{
	std::vector<fluid> fluids;
	double fluid::max_force;
	const graphics::color particle::default_color = graphics::colors::blue;
	void init()
	{
		BOOST_FOREACH(config::ptree::value_type& pt , config::cfg.get_child("fluids"))
		{
			math::vec pos(pt.second.get<double>("start_x"), pt.second.get<double>("start_y"));
			fluid f(pt.second.get<int>("particles_x"), pt.second.get<int>("particles_y"), pt.second.get<double>("spacing_x"), pt.second.get<double>("spacing_y"), pos);
			f.mass=pt.second.get<double>("mass");
			f.density=pt.second.get<double>("density");
			f.rest_density=pt.second.get<double>("rest_density");
			f.damping=pt.second.get<double>("damping");
			f.eta=pt.second.get<double>("eta");
			f.sigma=pt.second.get<double>("sigma");
			f.k=pt.second.get<double>("k");
			f.smoothing_length=pt.second.get<double>("smoothing_length");
			f.gradient_length_treshold=pt.second.get<double>("gradient_length_treshold");
			f.glass_pushback_distance=pt.second.get<double>("glass_pushback_distance");
			f.glass_pushback_minimum_pressure=pt.second.get<double>("glass_pushback_minimum_pressure");
			f.glass_density=pt.second.get<double>("glass_density");
			f.glass_viscosity=pt.second.get<double>("glass_viscosity");
			f.wall_damping=0.25*f.damping;
			math::calculate_common_factors(f.smoothing_length);
			f.grid=new particle_grid(f.smoothing_length);
			add_fluid(f);
		}
	}

	int add_fluid( fluid& f )
	{
		fluids.push_back(f);
		return fluids.size()-1;
	}

	void step( double dt )
	{
		for(fluid& f : fluids)
		{
			f.step(dt);
		}
	}


	fluid::~fluid()
	{

	}

	fluid::fluid(int particles_x, int particles_y, double spacing_x, double spacing_y, math::vec pos)
	{
		n=particles_x*particles_y;
		particles=new particle[n];
		int k=0;
		for(int i=0;i<particles_x;i++)
		{
			for(int j=0;j<particles_y;j++)
			{
				particles[k].pos.Set(pos.x+spacing_x*i, pos.y+spacing_y*j);
				particles[k].v=math::vec(0, 0);
				k++;
			}
		}
	}

	std::vector<particle_pair> pairs_cache;

	void fluid::step( double dt )
	{
		prepare_step();

		/*STEP 1: calculate densities*/
		calculate_densities();
		glass_common_update(&fluid::calculate_glass_fluid_densities);

		/*if(simulator::detailed_logging)
		{
		printf("densities: ");
		for(int i=0;i<n;i++)
		{
		printf("%lf ", particles[i].density);
		if(particles[i].density!=density) system("pause");
		}
		printf("\n");
		}*/

		/*STEP 2: calculate forces*/
		calculate_forces();
		glass_common_update(&fluid::calculate_glass_fluid_forces);
		/*STEP 3: move particles*/
		update_positions(dt);
		set_bounding_particle_indices();
	}

	void fluid::calculate_densities( int x, int y)
	{
		const int16_t neighbor_dx[]={/*3x3*/ /*0,*/ 1, 1, 0, -1, -1, -1, 0, 1, /*5x5 2, 2, 2, 2, 2, 1, 0, -1, -2, -2, -2, -2, -2, -1, 0, 1*/};
		const int16_t neighbor_dy[]={/*3x3*/ /*0,*/ 0, 1, 1, 1, 0, -1, -1, -1, /*5x5 -2, -1, 0, 1, 2, 2, 2, 2, 2, 1, 0, -1, -2, -2, -2, -2*/};
		const int n_neighbor_cells=sizeof(neighbor_dx)/sizeof(int16_t);

		size_t idx = grid->get_cell_index(x, y);
		if(idx==particle_grid::invalid_index)
			return;
		cell* c=grid->cells[idx];
		if(c==nullptr) return;

		for(int i=0;i<c->size();i++)
		{
			//same cell
			particle *a=(*c)[i], *b;
			for (int j = i+1; j < c->size(); j++)
			{
				b=(*c)[j];
				math::vec r=(a->pos)-(b->pos);
				double r_length=r.Length();
				if(r_length<smoothing_length)
				{
					pairs_cache.push_back(particle_pair(a, b, r_length));
					double additional_density=mass*math::w_poly6(r_length);
					/*if(r.x!=0)	
					{
					logger::log("(%lf, %lf) -> %lf", r.x, r.y, additional_density);
					system("pause");
					//__debugbreak();
					}*/
					//#pragma omp atomic
					a->density+=additional_density;
					//#pragma omp atomic
					b->density+=additional_density;
				}

			}

			//neighbor cells
			for(int ni=0;ni<n_neighbor_cells;ni++)
			{
				size_t idx = grid->get_cell_index(x+neighbor_dx[ni], y+neighbor_dy[ni]);
				if(idx==particle_grid::invalid_index) continue;
				cell* c2=grid->cells[idx];
				if(c2==nullptr) continue;
				for(int j=0;j<c2->size();j++)
				{
					b=(*c2)[j];
					math::vec r=(a->pos)-(b->pos);
					double r_length=r.Length();
					if(r_length<=smoothing_length)
					{
						pairs_cache.push_back(particle_pair(a, b, r_length));
						double additional_density=mass*math::w_poly6(r_length);
						//#pragma omp atomic
						a->density+=additional_density;
						//#pragma omp atomic
						b->density+=additional_density;
					}
				}
			}
		}
	}

	void fluid::calculate_densities()
	{
		//#pragma omp parallel for schedule(dynamic)
		for(int x=0;x<grid->n_cells_x;x++)
		{
			for(int y=0;y<grid->n_cells_y;y++)
			{
				calculate_densities(x, y);
			}
		}
	}

	void fluid::prepare_step()
	{
		max_force=DBL_MIN;
		grid->fill(particles, n);
		pairs_cache.clear();
		for(int i=0;i<n;i++)
		{
			particles[i].col=particle::default_color;
			particles[i].density=density;
			particles[i].forces.Set(0,0);
			particles[i].color_field_gradient=math::vec(0, 0);
			particles[i].color_field_laplacian=0;
		}
	}

	void fluid::draw()
	{
		for(int i=0;i<n;i++)
		{
			particles[i].draw();
		}
	}

	void fluid::calculate_forces()
	{
		size_t pairs_cache_size=pairs_cache.size();
		//#pragma omp parallel for schedule(dynamic) 
		for(int i=0;i<pairs_cache_size;i++)
		{
			calculate_forces(pairs_cache[i]);
		}
	}

	void fluid::calculate_forces( particle_pair p )
	{
		particle *a=p.p1, *b=p.p2;
		double r_length=p.dist;
		math::vec r=a->pos-b->pos;

		//pressure
		double 
			pressure_a=k*(a->density-rest_density),
			pressure_b=k*(b->density-rest_density);
		math::vec common_pressure_term(-mass*(pressure_a+pressure_b)*0.5*math::w_spiky_gradient(r, r_length));
		math::vec 
			pressure_force_a=common_pressure_term/b->density,
			pressure_force_b=-common_pressure_term/a->density;

		//viscosity
		math::vec common_viscosity_term=mass*eta*((b->v)-(a->v))*math::w_viscosity_laplacian(r_length);
		math::vec
			viscosity_force_a=common_viscosity_term / b->density,
			viscosity_force_b=common_viscosity_term / a->density;

		//surface tension
		math::vec 
			common_color_field_gradient_term=mass*math::w_poly6_gradient(r, r_length);

		a->color_field_gradient.atomic_increment(common_color_field_gradient_term/b->density);
		b->color_field_gradient.atomic_increment(common_color_field_gradient_term/a->density);

		double 
			common_color_field_laplacian_term = mass * math::w_poly6_laplacian(r_length);
		//#pragma omp atomic
		a->color_field_laplacian+=common_color_field_laplacian_term / b->density;
		//#pragma omp atomic
		b->color_field_laplacian+=common_color_field_laplacian_term / a->density;
		a->forces.atomic_increment(viscosity_force_a+pressure_force_a);
		b->forces.atomic_increment(viscosity_force_b+pressure_force_b);
	}

	void fluid::update_positions( double dt )
	{
		//#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<n;i++)
		{
			//surface tension
			math::vec f_surface_tension(0, 0);
			double color_field_gradient_length=particles[i].color_field_gradient.Length();
			if(color_field_gradient_length>gradient_length_treshold)
			{
				f_surface_tension=-sigma*particles[i].color_field_laplacian*particles[i].color_field_gradient/color_field_gradient_length;
			}
			math::vec force=particles[i].forces;
			math::vec acceleration =  force / particles[i].density+physics::gravity;
			max_force=std::max(max_force, (force+physics::gravity*particles[i].density).Length());
			math::vec old_v=particles[i].v;
			particles[i].v=pow(damping, dt)*particles[i].v+acceleration*dt;

			math::vec newpos = particles[i].pos+0.5*(old_v+particles[i].v)*dt;

			double deltat=dt;

			enforce_glass(&particles[i], deltat, newpos);
			enforce_walls(&particles[i], deltat, newpos);

			//logger::log("(%.2lf, %.2lf) ", particles[i].pos.x, particles[i].pos.y);
		}
	}


	void fluid::set_bounding_particle_indices()
	{
		north_particle_idx=0, south_particle_idx=0, east_particle_idx=0, west_particle_idx=0;
		//#pragma omp parallel for schedule(dynamic)
		for(int i=1;i<n;i++)
		{
			if(particles[i].pos.x<particles[west_particle_idx].pos.x)
			{
				west_particle_idx=i;
			}
			if(particles[i].pos.y<particles[south_particle_idx].pos.y)
			{
				south_particle_idx=i;
			}
			if(particles[i].pos.x>particles[east_particle_idx].pos.x)
			{
				east_particle_idx=i;
			}
			if(particles[i].pos.y>particles[north_particle_idx].pos.y)
			{
				north_particle_idx=i;
			}
		}
		particles[south_particle_idx].col=particles[north_particle_idx].col=particles[east_particle_idx].col=particles[west_particle_idx].col=graphics::colors::red;
	}

	void draw()
	{
		for(fluid& f:fluids)
		{
			glColor3f(0, 0, 1);
			f.draw();
		}
	}

};