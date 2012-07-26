#include "simulator.h"
#include "logger.h"
#include "graphics.h"
#include "config.h"
#include <omp.h>

#define SHITY_SDL_DEFINE

using namespace std;
using namespace simulator;

static bool running = true;
bool	simulator::detailed_logging=false;
int simulator::mousex=0, simulator::mousey=0, simulator::width, simulator::height;

void init()
{
	config::read("world.json");
	simulator::width=config::cfg.get<int>("dimensions.w");
	simulator::height=config::cfg.get<int>("dimensions.h");
	if(!graphics::init())
	{	
		cerr<<"graphics init failed\n";
		exit(1);
	}
	if(!physics::init())
	{
		cerr<<"physics init failed\n";
		exit(1);
	}

	omp_set_num_threads(8);

	int x, y;
	SDL_PumpEvents();
	SDL_GetMouseState(&x, &y);
	physics::add_black_hole(physics::black_hole(math::vec(x, y), 10000000));

}

void destroy()
{
	graphics::destroy();
	physics::destroy();
}


void handle_event(SDL_Event &evt, float dt)
{
	if (evt.type == SDL_QUIT)
		running = false;
	if(evt.type == SDL_KEYDOWN)
	{
		if(evt.key.keysym.sym==SDLK_ESCAPE)
		{
			running=false;
		}
	}
	if(evt.type == SDL_MOUSEBUTTONUP)
	{
		mousex=evt.button.x; mousey=simulator::height-evt.button.y;
		logger::log("clicked at %d, %d", mousex, mousey);
		if(evt.button.button==SDL_BUTTON_LEFT)
		{
			physics::blackholes_enabled=!physics::blackholes_enabled;
			logger::log("blackholes enabled: %s", physics::blackholes_enabled?"YES":"NO");
		}
		if(evt.button.button==SDL_BUTTON_RIGHT)
		{
			detailed_logging=!detailed_logging;
			logger::log("logging enabled: %s", detailed_logging?"YES":"NO");
		}
	}
	if(evt.type==SDL_MOUSEMOTION)
	{
		mousex=evt.button.x; mousey=simulator::height-evt.button.y;
		physics::black_holes[0].pos.Set(mousex/2, mousey/2);
	}
}

void step(double dt)
{
	physics::step(dt);
}

void render()
{
	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	physics::draw();
}

int main( int argc, char *argv[] )
{
	init();

	SDL_Event evt;
	uint64_t oldtime=SDL_GetPerformanceCounter();

	while(running)
	{
		uint64_t newtime=SDL_GetPerformanceCounter();
		double timestep=double(newtime-oldtime)/SDL_GetPerformanceFrequency();
		oldtime=newtime;
		if(timestep>very_long_time)
		{
			logger::log("last frame too slow (%lfs)", timestep);
			timestep=expected_time;
		}

		while (SDL_PollEvent(&evt))
			handle_event(evt, timestep);

		if(timestep>0) //nesto se desilo, jer se za dt=0 ne moze nista desiti...
		{
			step(timestep);
		}
		render();

		SDL_GL_SwapBuffers();
		//SDL_Delay(1);
	}

	destroy();
	return 0;
}