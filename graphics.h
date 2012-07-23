#pragma once

#include "simulator.h"
#include "logger.h"
#include "math/vector.h"

namespace graphics
{
	const double PI = 3.1415926535897932384626433832795;
	bool init();
	bool destroy();
	void proseri();
	void circle(math::vec p, double r=10);
}