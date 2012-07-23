#include "logger.h"

namespace logger
{
	void log(const char* format, ...)
	{
		va_list args;
		va_start(args, format);
		fprintf(stderr, "%lld: ", simulator::SDL_GetPerformanceCounter());
		vfprintf(stderr, format, args);
		fprintf(stderr, "\n");
		va_end(args);
	}
};