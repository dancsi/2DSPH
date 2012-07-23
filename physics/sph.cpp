#include "sph.h"

#include "simulator.h"
#include "math/vector.h"
#include "math/kernels.h"
#include "graphics.h"
#include "config.h"
#include "boost/foreach.hpp"
#include "physics.h"

namespace sph
{
	std::vector<fluid> fluids;

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
				particles[k].flags=0;
				k++;
			}
		}
		int idx=0;
		particles[idx].track();
		logger::log("following particle %d", idx);
	}

	void fluid::step( double dt )
	{
		static unsigned long long tickcount=0;
		tickcount++;
		if(tickcount%100==0) logger::log("%llu ticks", tickcount);
		double h=smoothing_length;
		prepare_step();

		/*STEP 1: calculate densities*/
		for(int i=0;i<n;i++)
		{
			//logger::log("(%.2lf, %.2lf)", particles[i].pos.x, particles[i].pos.y);
			for(int j=0;j<n;j++)
			{
				if(j!=i)
				{
					math::vec r= particles[i].pos-particles[j].pos;
					if(r.Length()<h)
					{
						particles[i].density+=math::w_poly6(r, h);
					}
				}
			}
			particles[i].density*=mass;
			//logger::log("density: %lf", particles[i].density);
		}

		/*STEP 2: calculate forces*/
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				if(j!=i)
				{
					math::vec r= particles[i].pos-particles[j].pos;
					if(i==0&&j<50) logger::log("%d: %s", j, std::string(r).c_str());
					if(r.Length()<=h)
					{
						double pressure_i=k*(particles[i].density-rest_density);
						double pressure_j=k*(particles[j].density-rest_density);
						//logger::log("pressure_i: %.2lf, pressure_j: %.2lf", pressure_i, pressure_j);

						if(i==0&&j==30) 
							logger::log("FOUND %d", j);

						particles[i].f_pressure+=mass*(pressure_i+pressure_j)/(2*particles[j].density)*math::w_gradient_spiky(r, h);
						particles[i].f_viscosity+=eta*mass*(particles[j].v-particles[i].v)/particles[j].density*math::w_laplacian_viscosity(r, h)*1000;
						particles[i].colour_field_gradient+=mass/particles[j].density*math::w_gradient_poly6(r, h)*1e7;
						particles[i].colour_field_laplacian+=mass/particles[j].density*math::w_laplacian_poly6(r, h)*1e7;
					}
				}
			}
		}

		/*STEP 3: move particles*/
		for(int i=0;i<n;i++)
		{
			//surface tension
			math::vec f_surface_tension(0, 0);
#if 1
			if(particles[i].colour_field_gradient.Length()>gradient_length_treshold)
			{
				f_surface_tension=-sigma*particles[i].colour_field_laplacian*particles[i].colour_field_gradient/particles[i].colour_field_gradient.Length();
			}
#endif
			math::vec force=math::vec(0, 0);

			if(physics::blackholes_enabled)
			{
				for(physics::black_hole& b: physics::black_holes)
				{
					force+=b.get_force_field(particles[i].pos);
				}
			}

			force+=f_surface_tension+particles[i].f_pressure+particles[i].f_viscosity;
			math::vec acceleration =  force / particles[i].density+physics::gravity;
			if(i==0)
			{
				logger::log("surface_tension: %s, pressure: %s, viscosity: %s, acceleration: %s", 
					std::string(f_surface_tension).c_str(), 
					std::string(particles[i].f_pressure).c_str(),
					std::string(particles[i].f_viscosity).c_str(),
					std::string(acceleration).c_str());
			}
			particles[i].v+=acceleration*dt;

			math::vec newpos = particles[i].pos+particles[i].v*dt;
			if(newpos.x<0)
			{
				double passed_time=(particles[i].pos.x-0)/particles[i].v.x;
				particles[i].v.x*=-damping; //reversing direction
				newpos.x=(dt-passed_time)*particles[i].v.x;
			}
			else if(newpos.x>simulator::width)
			{
				double passed_time=(simulator::width-particles[i].pos.x)/particles[i].v.x;
				particles[i].v.x*=-damping; //reversing direction
				newpos.x=simulator::width-(dt-passed_time)*particles[i].v.x;
			}
			if(newpos.y<0)
			{
				double passed_time=(particles[i].pos.y-0)/particles[i].v.y;
				particles[i].v.y*=-damping; //reversing direction
				newpos.y=(dt-passed_time)*particles[i].v.y;
			}
			else if(newpos.y>simulator::height)
			{
				double passed_time=(simulator::height-particles[i].pos.y)/particles[i].v.y;
				particles[i].v.y*=-damping; //reversing direction
				newpos.y=simulator::height-(dt-passed_time)*particles[i].v.y;
			}
			particles[i].pos=newpos;
			//logger::log("(%.2lf, %.2lf) ", particles[i].pos.x, particles[i].pos.y);
		}
		/*char buf[100];
		sprintf(buf, "dowwin (x=%.2lf, y=%.2lf)\n", particles[0].v.x, particles[0].v.y);
		SDL_WM_SetCaption(buf, nullptr);*/
	}

	void fluid::prepare_step()
	{
		for(int i=0;i<n;i++)
		{
			particles[i].density=density;
			particles[i].f_pressure=0;
			particles[i].f_viscosity=0;
			particles[i].colour_field_gradient=0;
			particles[i].colour_field_laplacian=0;
		}
	}

	void fluid::draw()
	{
		for(int i=0;i<n;i++)
		{
			particles[i].draw();
		}
	}

	void draw()
	{
		for(fluid& f:fluids)
		{
			glColor3f(0, 0, 1);
			f.draw();
		}
	}


	void fluid::particle::draw()
	{
		if(flags&TRACK)
		{
			glColor3f(1, 0, 0);
			glTranslatef(0, 0, .5);
		}
		else
		{
			glColor3f(0, 0, 1);
		}
		graphics::circle(pos, 2);
		if(flags&TRACK)
		{
			glTranslatef(0, 0, -.5);
		}
	}

	void fluid::particle::track()
	{
		flags|=TRACK;
	}

};