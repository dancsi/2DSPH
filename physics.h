#pragma once

#include "simulator.h"
#include "math/vector.h"
#include <vector>

#include "physics/line.h"
#include "physics/black_hole.h"

namespace physics
{
	extern math::vec gravity;
	extern std::vector<line> scenery;
	extern std::vector<black_hole> black_holes;
	extern bool blackholes_enabled;
	extern double cell_size;
	bool init();
	bool destroy();
	void draw();
	void step(double dt);
	int add_black_hole(black_hole& bh);
	uint16_t hash(math::vec& pos);
	uint16_t get_neighbor_cell(uint16_t hash, uint16_t dx, uint16_t dy);
}