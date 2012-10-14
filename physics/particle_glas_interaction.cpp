#pragma once

#include "sph.h"
#include "math/kernels.h"

namespace sph
{
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
	void fluid::enforce_walls( particle* p, double& dt, math::vec& newpos )
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

		p->pos=newpos;
	}
	void fluid::enforce_glass( particle* p, double& dt, math::vec& newpos )
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

		bool critical_particle=false;

		if(newpos.x<190.0) critical_particle=true;
		if(newpos.x>210.0) critical_particle=true;
		critical_particle=false;
		if(critical_particle)
		{
			logger::log("critical particle %p\n", p);
			logger::log("pos=%s, v=%s, newpos=%s", STR(p->pos), STR(p->v), STR(newpos));
		}

		//#define DEBUGGING
		bool r=true;
		while(dt>0 && r)
		{
			math:: line l2=math::line(p->pos, newpos);
			auto it =min_element(physics::scenery.begin(), physics::scenery.end(), [&] (math::line& a, math::line& b) -> bool {
				math::vec inter_a(0, 0), inter_b(0, 0);
				bool r_a=math::intersection(a, l2, inter_a), r_b=math::intersection(b, l2, inter_b);
				if(!r_a && !r_b)
				{
					return (p->pos-inter_a).LengthSq()<(p->pos-inter_b).LengthSq();;
				}
				else
				{
					if(!r_b) return true;
					if(!r_a) return false;
					return (p->pos-inter_a).LengthSq()<(p->pos-inter_b).LengthSq();
				}
			});
			math::line l1=*it; 
			math::vec inter(0, 0);
			r=math::intersection(l1, l2, inter);
			if(r)
			{

				double passed_time=sqrt((p->pos-inter).LengthSq()/p->v.LengthSq());
				math::vec oldv=p->v;
				p->v=p->v.deflect(l1)*wall_damping;
				dt-=passed_time;
				math::vec oldpos=p->pos;
				p->pos=inter;
				newpos=inter+(p->v)*(dt);
				if(critical_particle)
				{
					logger::log("deflected to pos=%s, v=%s, newpos=%s", STR(p->pos), STR(p->v), STR(newpos));
					logger::log("force: %s", STR(mass*(1.0/passed_time)*(p->v-oldv)));
				}
			}
		}
		p->pos=newpos;
		if(critical_particle)
		{
			if(newpos.x>210.0+1e-3)
				__debugbreak();
			if(newpos.x<190.0-1e-3)
				__debugbreak();
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
		//forces*=100;
		p->forces+=forces;
		//if(simulator::detailed_logging) logger::log("\tforces delta: (%lf, %lf)",forces.x, forces.y);

	}
	void fluid::calculate_glass_fluid_densities( particle* p, math::vec r )
	{
		double additional_density=mass*math::w_poly6(r.Length());
		p->density+=additional_density;
	}
}