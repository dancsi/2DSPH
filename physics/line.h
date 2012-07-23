#pragma once

#include "obj.h"
#include "math/vector.h"

namespace physics
{
	class line: public obj
	{
	public:
		line();
		line(math::vec& a, math::vec& b);
		bool collides(line& l);
		void draw();
		math::vec a, b;
	};

}