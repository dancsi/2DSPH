#pragma once

#include "math/vector.h"
namespace physics
{
	struct black_hole
	{
		math::vec pos;
		double force;
		black_hole();
		black_hole(math::vec pos, double force=2);
		math::vec get_force_field(math::vec point);
	};
}