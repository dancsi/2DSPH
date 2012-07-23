#include "black_hole.h"

namespace physics 
{
	black_hole::black_hole()
	{
		force=0;
	}

	black_hole::black_hole( math::vec pos, double force/*=100*/ )
	{
		black_hole::pos=pos;
		black_hole::force=force;
	}

	math::vec black_hole::get_force_field( math::vec point )
	{
		return force*(pos-point);
	}
}