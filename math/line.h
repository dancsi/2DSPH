#pragma once

#include "math/vector.h"

namespace math
{
	class line
	{
	public:
		line();
		line(math::vec& a, math::vec& b);
		operator math::vec();
		bool collides(line& l);
		void draw();
		math::vec a, b;
	};
	math::vec point_projection(math::vec line, math::vec point);
	bool intersection(const line &_ls1,const line &_ls2, math::vec &_ip /*= 0*/ );
	bool point_is_on_line_segment(line& l, math::vec point);
}