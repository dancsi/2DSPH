#include "physics.h"

#include "config.h"
#include "graphics.h"
#include "line.h"
#include "sph.h"
#include <boost/foreach.hpp>
#include <cmath>
#include <vector>

namespace physics
{
	math::vec gravity;
	std::vector<line> scenery;
	std::vector<black_hole> black_holes;
	bool blackholes_enabled=false;
	double cell_size=4;
	bool init()
	{
		using boost::property_tree::ptree;
		physics::gravity=math::vec(0, config::cfg.get<double>("constants.gravity"));
		BOOST_FOREACH(ptree::value_type& obj, config::cfg.get_child("objects.fixed"))
		{
			logger::log("%s", obj.second.get<std::string>("type").c_str());
			if(obj.second.get<std::string>("type")=="line")
			{
				std::vector<math::vec> vertices;
				BOOST_FOREACH(ptree::value_type& vert, obj.second.get_child("vertices"))
				{
					vertices.push_back(math::vec(vert.second.get<float>("x"), vert.second.get<float>("y")));
				}
				scenery.push_back(line(vertices[0], vertices[1]));
			}
		}

		sph::init();
		return true;
	}

	void step(double dt)
	{
		sph::step(dt);
	}

	bool destroy()
	{
		return 1;
	}

	void draw()
	{
		for(line& l: scenery)
		{
			l.draw();
		}
		sph::draw();
	}

	int add_black_hole( black_hole& bh )
	{
		black_holes.push_back(bh);
		return black_holes.size()-1;
	}

	uint16_t hash( math::vec& pos )
	{
		return ((uint16_t(pos.x)&0xff)<<8)+(uint16_t(pos.y)&0xff);
	}

	uint16_t get_neighbor_cell( uint16_t hash, uint16_t dx, uint16_t dy )
	{
		//logger::log("first component: 0x%x, second: 0x%x, retval: 0x%x", ((hash>>8+dx)&0xff)<<8, ((hash&0xff)+dy)&0xff, ((((hash>>8)+dx)&0xff)<<8) + (((hash&0xff)+dy)&0xff));
		return ((((hash>>8)+dx)&0xff)<<8)+(((hash&0xff)+dy)&0xff);
	}

}