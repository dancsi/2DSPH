#include "graphics.h"

namespace graphics
{
	bool init()
	{
		SDL_Init( SDL_INIT_EVERYTHING );
		SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
		logger::log("init sdl (%d, %d)", simulator::width, simulator::height);
		SDL_SetVideoMode( simulator::width, simulator::height, 32, SDL_OPENGL );

		glClearColor( 0, 0, 0, 0 );

		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		glOrtho( 0, simulator::width, 0, simulator::height, -1, 1 );

		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();

		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glEnable(GL_LINE_SMOOTH);
		//glEnable(GL_POLYGON_SMOOTH);
		//glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);

		SDL_WM_SetCaption( "dowwin", NULL );

		return 1;
	}

	bool destroy()
	{
		return 1;
	}

	void circle(float x, float y, float r=10, int n=10)
	{
		//logger::log("circle (%f, %f)", x, y);
		//glColor3f(1, 0, 0);
		glBegin(GL_TRIANGLE_FAN);
		double step=2*PI/n, angle=0;
		glVertex2f(x, y);
		for(int i=0;i<n;i++, angle+=step)
		{
			glVertex2d(x+cos(angle)*r, y+sin(angle)*r);
		}
		glEnd();
	}

	void circle( math::vec p, double r/*=10*/ )
	{
		circle(p.x, p.y, r);
	}

	void proseri()
	{
		logger::log("proseravam");
		glBegin(GL_TRIANGLES);
		glColor3f(1, 0, 0);
		glVertex2f(0, 0);
		glVertex2f(100, 100);
		glVertex2f(100, 000);
		glEnd();
		circle(simulator::mousex, simulator::mousey );
	}
}