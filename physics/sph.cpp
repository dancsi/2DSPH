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

		logger::log("northmost particle: %d", north_particle_idx);
		//particles[north_particle_idx].dump();
		//particles[south_particle_idx].dump();
		//particles[west_particle_idx].dump();
		//particles[east_particle_idx].dump();
		//if(particles[east_particle_idx])
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

	void fluid::enforce_walls( particle* p, double dt, math::vec& newpos )
	{
		if(newpos.x<0)
		{
			double passed_time=(p->pos.x-0)/p->v.x;
			p->v.x*=-wall_damping; //reversing direction
			newpos.x=(dt-passed_time)*p->v.x;
		}
		else if(newpos.x>simulator::width)
		{
			double passed_time=(simulator::width-p->pos.x)/p->v.x;
			p->v.x*=-wall_damping; //reversing direction
			newpos.x=simulator::width-(dt-passed_time)*p->v.x;
		}
		if(newpos.y<0)
		{
			double passed_time=(p->pos.y-0)/p->v.y;
			p->v.y*=-wall_damping; //reversing direction
			newpos.y=(dt-passed_time)*p->v.y;
		}
		else if(newpos.y>simulator::height)
		{
			double passed_time=(simulator::height-p->pos.y)/p->v.y;
			p->v.y*=-1; //reversing direction
			newpos.y=simulator::height-(dt-passed_time)*p->v.y;
		}
	}

	void fluid::enforce_glass( particle* p, double dt, math::vec& newpos )
	{
		/*
		if(newpos.x<190)
		{
			double passed_time=(p->pos.x-190)/p->v.x;
			p->v.x*=-wall_damping; //reversing direction
			newpos.x=190+(dt-passed_time)*p->v.x;
		}
		else if(newpos.x>210)
		{
			double passed_time=(210-p->pos.x)/p->v.x;
			p->v.x*=-wall_damping; //reversing direction
			newpos.x=210-(dt-passed_time)*p->v.x;
		}
		*/
		for(int j=0;j<physics::scenery.size();j++)
		{
			math::line l1=physics::scenery[j], l2=math::line(p->pos, newpos);
			math::vec inter(0, 0);
			bool r=math::intersection(l1, l2, inter);
			if(r)
			{
				double passed_time=sqrt((p->pos-inter).LengthSq()/p->v.LengthSq());
				p->v=p->v.deflect(l1)*wall_damping;
				newpos=inter+p->v*(dt-passed_time);
			}
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
#if 0
			if(physics::blackholes_enabled)
			{
				for(physics::black_hole& b: physics::black_holes)
				{
					force+=b.get_force_field(particles[i].pos);
				}
			}
#endif
			math::vec acceleration =  force / particles[i].density+physics::gravity;
			max_force=std::max(max_force, (force+physics::gravity*particles[i].density).Length());
			/*if(simulator::detailed_logging)
			{
			logger::log("%d: surface_tension: %s, pressure: %s, viscosity: %s, acceleration: %s", i, std::string(f_surface_tension).c_str(), std::string(particles[i].f_pressure).c_str(),	std::string(particles[i].f_viscosity).c_str(), std::string(acceleration).c_str());
			}*/
			math::vec old_v=particles[i].v;
			particles[i].v=pow(damping, dt)*particles[i].v+acceleration*dt;

			math::vec newpos = particles[i].pos+0.5*(old_v+particles[i].v)*dt;

			enforce_walls(&particles[i], dt, newpos);
			enforce_glass(&particles[i], dt, newpos);

			particles[i].pos=newpos;
			//logger::log("(%.2lf, %.2lf) ", particles[i].pos.x, particles[i].pos.y);
		}
	}

	void fluid::glass_common_update(boost::function<void (fluid*, particle*, math::vec)> f)
	{
		static const double pb_dist=std::min(glass_pushback_distance, smoothing_length);
		//logger::log("glass pushback: %lf", pb_dist);
		//#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<physics::scenery.size();j++)
			{
				math::line l=physics::scenery[j];
				math::vec proj=math::point_projection(l, particles[i].pos);
				if(point_is_on_line_segment(l, proj))
				{
					double wall_dist=(particles[i].pos-proj).Length();
					if(wall_dist<=pb_dist)
					{
						//logger::log("collision waiting to happen between (%.2lf, %.2lf) and ((%.2lf, %.2lf) -> (%.2lf, %.2lf)", particles[i].pos.x, particles[i].pos.y, l.a.x, l.a.y, l.b.x, l.b.y);
						f(this, &particles[i], -(particles[i].pos-proj));
					}
				}
			}
			/*
			//glass
			//vertical
			double wall_dist=pb_dist-particles[i].pos.y;

			if(wall_dist>0)
			{
			//if(simulator::detailed_logging) logger::log("glass bottom:");
			f(this, &particles[i], math::vec(0, wall_dist));
			}
			//horizontal

			wall_dist=190.0+pb_dist-particles[i].pos.x;
			if(wall_dist>0)
			{
			//if(simulator::detailed_logging) logger::log("glass sides:");
			f(this, &particles[i], math::vec(wall_dist, 0));
			}
			wall_dist=particles[i].pos.x-(210.0-pb_dist);
			if(wall_dist>0)
			{
			//if(simulator::detailed_logging) logger::log("glass sides:");
			f(this, &particles[i], math::vec(-wall_dist, 0));
			}
			*/
		}

	}

	void fluid::calculate_glass_fluid_forces( particle* p, math::vec r )
	{
		double pressure=k*(p->density-rest_density);
		if(pressure<glass_pushback_minimum_pressure)
		{
			//logger::log("pressure too small (%lf), resetting to (%lf)", pressure, glass_pushback_minimum_pressure);
			pressure=glass_pushback_minimum_pressure;
		}
		math::vec pressure_force=-mass*pressure*math::w_spiky_gradient(r, r.Length())/p->density*glass_density;
		math::vec v_difference=-p->v;
		math::vec viscosity_force = mass*eta*v_difference*math::w_viscosity_laplacian(r.Length())/p->density*glass_viscosity;
		//if(simulator::detailed_logging) logger::log("pressure: %s, viscosity: %s", ((std::string)(pressure_force)).c_str(), ((std::string)viscosity_force).c_str());
		math::vec forces=viscosity_force+pressure_force;
		p->forces+=forces;
		if(simulator::detailed_logging) logger::log("\tforces delta: (%lf, %lf)",forces.x, forces.y);

	}

	void fluid::calculate_glass_fluid_densities( particle* p, math::vec r )
	{
		double additional_density=mass*math::w_poly6(r.Length());
		p->density+=additional_density;
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


	void particle::draw()
	{
		const float force_scale=30;
		col.set_current();
		if(col!=particle::default_color)
		{
			col=particle::default_color;
			glTranslatef(0, 0, 0.5);
			graphics::circle(pos, 1);
			glTranslatef(0, 0, -0.5);
		}
		else
		{
			graphics::circle(pos, 1);
		}
		/*glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		glVertex2f(pos.x, pos.y);
		glVertex2f(pos.x+(forces.x/fluid::max_force)*force_scale, pos.y+(forces.y/fluid::max_force)*force_scale);
		glEnd();*/
	}

	void particle::dump()
	{
		logger::log("\t{pos: (%.2lf, %.2lf), v: (%.2lf, %.2lf), force: (%.2lf, %.2lf)}", pos.x, pos.y, v.x, v.y, forces.x, forces.y);
	}


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


};