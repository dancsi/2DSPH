#pragma once

#include "math/vector.h"

namespace physics
{
	class obj
	{
	public:
		obj(void);
		~obj(void);
		virtual bool collides(obj& other);
	};
}