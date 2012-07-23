#include "config.h"
#include "logger.h"


namespace config
{
	ptree cfg;
	void read(const char* fname)
	{
		read_json(fname, cfg);
	}
}