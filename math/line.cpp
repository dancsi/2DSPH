#include "line.h"
#include "graphics.h"

namespace math
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

	line::operator math::vec()
	{
		return b-a;
	}

	math::vec point_projection( vec line, vec point )
	{
		return ((line*point)/(line.LengthSq()))*line;
	}

	bool point_is_on_line_segment( line& l, vec point )
	{
		if(((l.a-l.b)*point)*((l.b-l.a)*point)>=0)
			return true;
	}

	bool intersection(const line &_ls1,const line &_ls2, math::vec &_ip /*= 0*/ )
	{
		float s1_x, s1_y, s2_x, s2_y;
		s1_x = _ls1.b.x - _ls1.a.x;     s1_y = _ls1.b.y - _ls1.a.y;
		s2_x = _ls2.b.x - _ls2.a.x;     s2_y = _ls2.b.y - _ls2.a.y;

		float s, t;
		s = (-s1_y * (_ls1.a.x -  _ls2.a.x) + s1_x * (_ls1.a.y - _ls2.a.y)) / (-s2_x * s1_y + s1_x * s2_y);
		t = ( s2_x * (_ls1.a.y -  _ls2.a.y) - s2_y * (_ls1.a.x  -  _ls2.a.x)) / (-s2_x * s1_y + s1_x * s2_y);

		if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
		{
			_ip.x = _ls1.a.x + (t * s1_x);
			_ip.y = _ls1.a.y + (t * s1_y);
			return 1;
		}

		return 0;
	}
}