#pragma once

#include "simulator.h"
#include "logger.h"
#include "math/vector.h"

namespace graphics
{
	const double PI = 3.1415926535897932384626433832795;
	
	struct color
	{
		float r, g, b;
		void set_current();
		const bool operator!=(const color& rhs) const;
	};

	namespace colors
	{
		const color red={1, 0, 0};
		const color blue={0, 0, 1};
		const color green={0, 1, 0};
	}

	bool init();
	bool destroy();
	void proseri();
	void circle(math::vec p, double r=10);
}