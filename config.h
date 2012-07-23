#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace config
{
	using boost::property_tree::ptree;
	extern ptree cfg;
	void read(const char* fname);
}