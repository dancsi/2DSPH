#pragma once

/* STD HEADERS */

#include <iostream>
#include <cstdio>

/* SDL */
#include <SDL.h>
#include <SDL_opengl.h>

/* INTERNAL INCLUDES */

#include "graphics.h"
#include "physics.h"

namespace simulator
{

	const double very_long_time=0.5;
	const double expected_time=0.0016;

	extern int height, width;

	extern int mousex, mousey;
	extern bool detailed_logging;

	inline int64_t SDL_GetPerformanceCounter()
	{
		LARGE_INTEGER r;
		QueryPerformanceCounter(&r);
		return r.QuadPart;
	}

	inline int64_t SDL_GetPerformanceFrequency()
	{
		LARGE_INTEGER r;
		QueryPerformanceFrequency(&r);
		return r.QuadPart;
	}


}