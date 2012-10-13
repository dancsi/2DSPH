#pragma once

#include "vec2.h"
#include <DirectXMath.h>
#include <DirectXColors.h>

#define STR(x) std::string(x).c_str()

namespace math
{
	using namespace DirectX;
	typedef vec2 vec;
	typedef DirectX::XMFLOAT3 color_t;
}