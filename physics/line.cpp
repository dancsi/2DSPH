#include "line.h"
#include "graphics.h"

namespace physics
{
	line::line() {}
	line::line(math::vec& _a, math::vec& _b): a(_a), b(_b) {}

	void line::draw()
	{
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);
		{
			glVertex2d(a.x, a.y);
			glVertex2d(b.x, b.y);
		}
		glEnd();
		//logger::log("line drawn from (%lf, %lf) to (%lf, %lf)", a.x, a.y, b.x, b.y);
	}

}