#include "particle.h"

namespace sph
{
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
}